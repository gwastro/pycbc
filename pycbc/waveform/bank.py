# Copyright (C) 2012  Alex Nitz, Josh Willis, Andrew Miller
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module provides classes that describe banks of waveforms
"""
import types
import logging
import os.path
import h5py
from copy import copy
import numpy as np
from ligo.lw import lsctables, utils as ligolw_utils
import pycbc.waveform
import pycbc.pnutils
import pycbc.waveform.compress
from pycbc import DYN_RANGE_FAC
from pycbc.types import FrequencySeries, zeros
import pycbc.io
from pycbc.io.ligolw import LIGOLWContentHandler
import hashlib
import warnings
import lal
import lalsimulation as lalsim


def sigma_cached(self, psd):
    """ Cache sigma calculate for use in tandem with the FilterBank class
    """
    if not hasattr(self, '_sigmasq'):
        from pycbc.opt import LimitedSizeDict
        self._sigmasq = LimitedSizeDict(size_limit=2**5)

    key = id(psd)
    if not hasattr(psd, '_sigma_cached_key'):
        psd._sigma_cached_key = {}

    if key not in self._sigmasq or id(self) not in psd._sigma_cached_key:
        psd._sigma_cached_key[id(self)] = True
        # If possible, we precalculate the sigmasq vector for all possible waveforms
        if pycbc.waveform.waveform_norm_exists(self.approximant):
            if not hasattr(psd, 'sigmasq_vec'):
                psd.sigmasq_vec = {}

            if self.approximant not in psd.sigmasq_vec:
                psd.sigmasq_vec[self.approximant] = \
                    pycbc.waveform.get_waveform_filter_norm(
                        self.approximant,
                        psd,
                        len(psd),
                        psd.delta_f,
                        self.min_f_lower
                    )

            if not hasattr(self, 'sigma_scale'):
                # Get an amplitude normalization (mass dependant constant norm)
                amp_norm = pycbc.waveform.get_template_amplitude_norm(
                                     self.params, approximant=self.approximant)
                amp_norm = 1 if amp_norm is None else amp_norm
                self.sigma_scale = (DYN_RANGE_FAC * amp_norm) ** 2.0

            curr_sigmasq = psd.sigmasq_vec[self.approximant]

            kmin = int(self.f_lower / psd.delta_f)
            self._sigmasq[key] = self.sigma_scale * \
                (curr_sigmasq[self.end_idx-1] - curr_sigmasq[kmin])

        else:
            if not hasattr(self, 'sigma_view'):
                from pycbc.filter.matchedfilter import get_cutoff_indices
                N = (len(self) -1) * 2
                kmin, kmax = get_cutoff_indices(
                        self.min_f_lower or self.f_lower, self.end_frequency,
                        self.delta_f, N)
                self.sslice = slice(kmin, kmax)
                self.sigma_view = self[self.sslice].squared_norm() * 4.0 * self.delta_f

            if not hasattr(psd, 'invsqrt'):
                psd.invsqrt = 1.0 / psd

            self._sigmasq[key] = self.sigma_view.inner(psd.invsqrt[self.sslice])
    return self._sigmasq[key]


# helper function for parsing approximant strings
def boolargs_from_apprxstr(approximant_strs):
    """Parses a list of strings specifying an approximant and where that
    approximant should be used into a list that can be understood by
    FieldArray.parse_boolargs.

    Parameters
    ----------
    apprxstr : (list of) string(s)
        The strings to parse. Each string should be formatted `APPRX:COND`,
        where `APPRX` is the approximant and `COND` is a string specifying
        where it should be applied (see `FieldArgs.parse_boolargs` for examples
        of conditional strings). The last string in the list may exclude a
        conditional argument, which is the same as specifying ':else'.

    Returns
    -------
    boolargs : list
        A list of tuples giving the approximant and where to apply them. This
        can be passed directly to `FieldArray.parse_boolargs`.
    """
    if not isinstance(approximant_strs, list):
        approximant_strs = [approximant_strs]
    return [tuple(arg.split(':')) for arg in approximant_strs]


def add_approximant_arg(parser, default=None, help=None):
    """Adds an approximant argument to the given parser.

    Parameters
    ----------
    parser : ArgumentParser
        The argument parser to add the argument to.
    default : {None, str}
        Specify a default for the approximant argument. Defaults to None.
    help : {None, str}
        Provide a custom help message. If None, will use a descriptive message
        on how to specify the approximant.
    """
    if help is None:
        help=str("The approximant(s) to use. Multiple approximants to use "
             "in different regions may be provided. If multiple "
             "approximants are provided, every one but the last must be "
             "be followed by a conditional statement defining where that "
             "approximant should be used. Conditionals can be any boolean "
             "test understood by numpy. For example, 'Apprx:(mtotal > 4) & "
             "(mchirp <= 5)' would use approximant 'Apprx' where total mass "
             "is > 4 and chirp mass is <= 5. "
             "Conditionals are applied in order, with each successive one "
             "only applied to regions not covered by previous arguments. "
             "For example, `'TaylorF2:mtotal < 4' 'IMRPhenomD:mchirp < 3'` "
             "would result in IMRPhenomD being used where chirp mass is < 3 "
             "and total mass is >= 4. The last approximant given may use "
             "'else' as the conditional or include no conditional. In either "
             "case, this will cause the last approximant to be used in any "
             "remaning regions after all the previous conditionals have been "
             "applied. For the full list of possible parameters to apply "
             "conditionals to, see WaveformArray.default_fields(). Math "
             "operations may also be used on parameters; syntax is python, "
             "with any operation recognized by numpy.")
    parser.add_argument("--approximant", nargs='+', type=str, default=default,
                        metavar='APPRX[:COND]',
                        help=help)


def parse_approximant_arg(approximant_arg, warray):
    """Given an approximant arg (see add_approximant_arg) and a field
    array, figures out what approximant to use for each template in the array.

    Parameters
    ----------
    approximant_arg : list
        The approximant argument to parse. Should be the thing returned by
        ArgumentParser when parsing the argument added by add_approximant_arg.
    warray : FieldArray
        The array to parse. Must be an instance of a FieldArray, or a class
        that inherits from FieldArray.

    Returns
    -------
    array
        A numpy array listing the approximants to use for each element in
        the warray.
    """
    return warray.parse_boolargs(boolargs_from_apprxstr(approximant_arg))[0]

def tuple_to_hash(tuple_to_be_hashed):
    """
    Return a hash for a numpy array, avoids native (unsafe) python3 hash function

    Parameters
    ----------
    tuple_to_be_hashed: tuple
        The tuple which is being hashed
        Must be convertible to a numpy array

    Returns
    -------
    int
        an integer representation of the hashed array
    """
    h = hashlib.blake2b(np.array(tuple_to_be_hashed).tobytes('C'),
                        digest_size=8)
    return np.fromstring(h.digest(), dtype=int)[0]


class TemplateBank(object):
    """Class to provide some basic helper functions and information
    about elements of a template bank.

    Parameters
    ----------
    filename : string
        The name of the file to load. Must end in '.xml[.gz]' or '.hdf'. If an
        hdf file, it should have a 'parameters' in its `attrs` which gives a
        list of the names of fields to load from the file. If no 'parameters'
        are found, all of the top-level groups in the file will assumed to be
        parameters (a warning will be printed to stdout in this case). If an
        xml file, it must have a `SnglInspiral` table.
    approximant : {None, (list of) string(s)}
        Specify the approximant(s) for each template in the bank. If None
        provided, will try to load the approximant from the file. The
        approximant may either be a single string (in which case the same
        approximant will be used for all templates) or a list of strings and
        conditionals specifying where to use the approximant. See
        `boolargs_from_apprxstr` for syntax.
    parameters : {None, (list of) sting(s)}
        Specify what parameters to load from the file. If None, all of the
        parameters in the file (if an xml file, this is all of the columns in
        the SnglInspiral table, if an hdf file, this is given by the
        parameters attribute in the file). The list may include parameters that
        are derived from the file's parameters, or functions thereof. For a
        full list of possible parameters, see `WaveformArray.default_fields`.
        If a derived parameter is specified, only the parameters needed to
        compute that parameter will be loaded from the file. For example, if
        `parameters='mchirp'`, then only `mass1, mass2` will be loaded from
        the file. Note that derived parameters can only be used if the
        needed parameters are in the file; e.g., you cannot use `chi_eff` if
        `spin1z`, `spin2z`, `mass1`, and `mass2` are in the input file.
    \**kwds :
        Any additional keyword arguments are stored to the `extra_args`
        attribute.

    Attributes
    ----------
    table : WaveformArray
        An instance of a WaveformArray containing all of the information about
        the parameters of the bank.
    has_compressed_waveforms : {False, bool}
        True if compressed waveforms are present in the the (hdf) file; False
        otherwise.
    indoc : {None, xmldoc}
        If an xml file was provided, an in-memory representation of the xml.
        Otherwise, None.
    filehandler : {None, pycbc.io.HFile}
        If an hdf file was provided, the file handler pointing to the hdf file
        (left open after initialization). Otherwise, None.
    extra_args : {None, dict}
        Any extra keyword arguments that were provided on initialization.
    """
    def __init__(self, filename, approximant=None, parameters=None,
                 **kwds):
        self.has_compressed_waveforms = False
        ext = os.path.basename(filename)
        if ext.endswith(('.xml', '.xml.gz', '.xmlgz')):
            self.filehandler = None
            self.indoc = ligolw_utils.load_filename(
                filename, False, contenthandler=LIGOLWContentHandler)
            self.table = lsctables.SnglInspiralTable.get_table(self.indoc)
            self.table = pycbc.io.WaveformArray.from_ligolw_table(self.table,
                columns=parameters)

            # inclination stored in xml alpha3 column
            names = list(self.table.dtype.names)
            names = tuple([n if n != 'alpha3' else 'inclination' for n in names])

            # low frequency cutoff in xml alpha6 column
            names = tuple([n if n!= 'alpha6' else 'f_lower' for n in names])
            self.table.dtype.names = names

        elif ext.endswith(('hdf', '.h5', '.hdf5')):
            self.indoc = None
            f = pycbc.io.HFile(filename, 'r')
            self.filehandler = f
            try:
                fileparams = list(f.attrs['parameters'])
            except KeyError:
                # just assume all of the top-level groups are the parameters
                fileparams = list(f.keys())
                logging.info("WARNING: no parameters attribute found. "
                    "Assuming that %s " %(', '.join(fileparams)) +
                    "are the parameters.")
            tmp_params = []
            # At this point fileparams might be bytes. Fix if it is
            for param in fileparams:
                try:
                    param = param.decode()
                    tmp_params.append(param)
                except AttributeError:
                    tmp_params.append(param)
            fileparams = tmp_params

            # use WaveformArray's syntax parser to figure out what fields
            # need to be loaded
            if parameters is None:
                parameters = fileparams
            common_fields = list(pycbc.io.WaveformArray(1,
                names=parameters).fieldnames)
            add_fields = list(set(parameters) &
                (set(fileparams) - set(common_fields)))
            # load
            dtype = []
            data = {}
            for key in common_fields+add_fields:
                data[key] = f[key][:]
                dtype.append((key, data[key].dtype))
            num = f[fileparams[0]].size
            self.table = pycbc.io.WaveformArray(num, dtype=dtype)
            for key in data:
                self.table[key] = data[key]
            # add the compressed waveforms, if they exist
            self.has_compressed_waveforms = 'compressed_waveforms' in f
        else:
            raise ValueError("Unsupported template bank file extension %s" %(
                ext))

        # if approximant is specified, override whatever was in the file
        # (if anything was in the file)
        if approximant is not None:
            # get the approximant for each template
            dtype = h5py.string_dtype(encoding='utf-8')
            apprxs = np.array(self.parse_approximant(approximant),
                              dtype=dtype)
            if 'approximant' not in self.table.fieldnames:
                self.table = self.table.add_fields(apprxs, 'approximant')
            else:
                self.table['approximant'] = apprxs
        self.extra_args = kwds
        self.ensure_hash()

    @property
    def parameters(self):
        """tuple: The parameters loaded from the input file.
        Same as `table.fieldnames`.
        """
        return self.table.fieldnames

    def ensure_hash(self):
        """Ensure that there is a correctly populated template_hash.

        Check for a correctly populated template_hash and create if it doesn't
        already exist.
        """
        fields = self.table.fieldnames
        if 'template_hash' in fields:
            return

        # The fields to use in making a template hash
        hash_fields = ['mass1', 'mass2', 'inclination',
                       'spin1x', 'spin1y', 'spin1z',
                       'spin2x', 'spin2y', 'spin2z',]

        fields = [f for f in hash_fields if f in fields]
        template_hash = np.array([tuple_to_hash(v) for v in zip(*[self.table[p]
                                  for p in fields])])
        if not np.unique(template_hash).size == template_hash.size:
            raise RuntimeError("Some template hashes clash. This should not "
                               "happen.")
        self.table = self.table.add_fields(template_hash, 'template_hash')

    def write_to_hdf(self, filename, start_index=None, stop_index=None,
                     force=False, skip_fields=None,
                     write_compressed_waveforms=True):
        """Writes self to the given hdf file.

        Parameters
        ----------
        filename : str
            The name of the file to write to. Must be a recognised HDF5
            file extension
        start_index : If a specific slice of the template bank is to be
            written to the hdf file, this would specify the index of the
            first template in the slice
        stop_index : If a specific slice of the template bank is to be
            written to the hdf file, this would specify the index of the
            last template in the slice
        force : {False, bool}
            If the file already exists, it will be overwritten if True.
            Otherwise, an OSError is raised if the file exists.
        skip_fields : {None, (list of) strings}
            Do not write the given fields to the hdf file. Default is None,
            in which case all fields in self.table.fieldnames are written.
        write_compressed_waveforms : {True, bool}
            Write compressed waveforms to the output (hdf) file if this is
            True, which is the default setting. If False, do not write the
            compressed waveforms group, but only the template parameters to
            the output file.

        Returns
        -------
        pycbc.io.HFile
            The file handler to the output hdf file (left open).
        """
        if not filename.endswith(('.hdf', '.h5', '.hdf5')):
            raise ValueError("Unrecoginized file extension")
        if os.path.exists(filename) and not force:
            raise IOError("File %s already exists" %(filename))
        f = pycbc.io.HFile(filename, 'w')
        parameters = self.parameters
        if skip_fields is not None:
            if not isinstance(skip_fields, list):
                skip_fields = [skip_fields]
            parameters = [p for p in parameters if p not in skip_fields]
        # save the parameters
        f.attrs['parameters'] = parameters
        write_tbl = self.table[start_index:stop_index]
        for p in parameters:
            f[p] = write_tbl[p]
        if write_compressed_waveforms and self.has_compressed_waveforms:
            for tmplt_hash in write_tbl.template_hash:
                compressed_waveform = pycbc.waveform.compress.CompressedWaveform.from_hdf(
                                        self.filehandler, tmplt_hash,
                                        load_now=True)
                compressed_waveform.write_to_hdf(f, tmplt_hash)
        return f

    def end_frequency(self, index):
        """ Return the end frequency of the waveform at the given index value
        """
        if hasattr(self.table[index], 'f_final'):
            return self.table[index].f_final

        return pycbc.waveform.get_waveform_end_frequency(
                                self.table[index],
                                approximant=self.approximant(index),
                                **self.extra_args)

    def parse_approximant(self, approximant):
        """Parses the given approximant argument, returning the approximant to
        use for each template in self. This is done by calling
        `parse_approximant_arg` using self's table as the array; see that
        function for more details."""
        return parse_approximant_arg(approximant, self.table)

    def approximant(self, index):
        """ Return the name of the approximant ot use at the given index
        """
        if 'approximant' not in self.table.fieldnames:
            raise ValueError("approximant not found in input file and no "
                "approximant was specified on initialization")
        apx = self.table["approximant"][index]
        if hasattr(apx, 'decode'):
            apx = apx.decode()
        return apx

    def __len__(self):
        return len(self.table)

    def template_thinning(self, inj_filter_rejector):
        """Remove templates from bank that are far from all injections."""
        if not inj_filter_rejector.enabled or \
                inj_filter_rejector.chirp_time_window is None:
            # Do nothing!
            return

        injection_parameters = inj_filter_rejector.injection_params.table
        fref = inj_filter_rejector.f_lower
        threshold = inj_filter_rejector.chirp_time_window
        m1= self.table['mass1']
        m2= self.table['mass2']
        tau0_temp, _ = pycbc.pnutils.mass1_mass2_to_tau0_tau3(m1, m2, fref)
        indices = []

        sort = tau0_temp.argsort()
        tau0_temp = tau0_temp[sort]

        for inj in injection_parameters:
            tau0_inj, _ = \
                pycbc.pnutils.mass1_mass2_to_tau0_tau3(inj.mass1, inj.mass2,
                                                       fref)
            lid = np.searchsorted(tau0_temp, tau0_inj - threshold)
            rid = np.searchsorted(tau0_temp, tau0_inj + threshold)
            inj_indices = sort[lid:rid]
            indices.append(inj_indices)

        indices_combined = np.concatenate(indices)
        indices_unique= np.unique(indices_combined)
        self.table = self.table[indices_unique]

    def ensure_standard_filter_columns(self, low_frequency_cutoff=None):
        """ Initialize FilterBank common fields

        Parameters
        ----------
        low_frequency_cutoff: {float, None}, Optional
            A low frequency cutoff which overrides any given within the
            template bank file.
        """

        # Make sure we have a template duration field
        if not hasattr(self.table, 'template_duration'):
            self.table = self.table.add_fields(np.zeros(len(self.table),
                                     dtype=np.float32), 'template_duration')

        # Make sure we have a f_lower field
        if low_frequency_cutoff is not None:
            if not hasattr(self.table, 'f_lower'):
                vec = np.zeros(len(self.table), dtype=np.float32)
                self.table = self.table.add_fields(vec, 'f_lower')
            self.table['f_lower'][:] = low_frequency_cutoff

        self.min_f_lower = min(self.table['f_lower'])
        if self.f_lower is None and self.min_f_lower == 0.:
            raise ValueError('Invalid low-frequency cutoff settings')


class LiveFilterBank(TemplateBank):
    def __init__(self, filename, sample_rate, minimum_buffer,
                       approximant=None, increment=8, parameters=None,
                       low_frequency_cutoff=None,
                       **kwds):

        self.increment = increment
        self.filename = filename
        self.sample_rate = sample_rate
        self.minimum_buffer = minimum_buffer
        self.f_lower = low_frequency_cutoff

        super(LiveFilterBank, self).__init__(filename, approximant=approximant,
                parameters=parameters, **kwds)
        self.ensure_standard_filter_columns(low_frequency_cutoff=low_frequency_cutoff)
        self.param_lookup = {}
        for i, p in enumerate(self.table):
            key =  (p.mass1, p.mass2, p.spin1z, p.spin2z)
            assert(key not in self.param_lookup) # Uh, oh, template confusion!
            self.param_lookup[key] = i

    def round_up(self, num):
        """Determine the length to use for this waveform by rounding.

        Parameters
        ----------
        num : int
            Proposed size of waveform in samples.

        Returns
        -------
        size: int
            The rounded size to use for the waveform buffer in samples.
            This is calculated using an internal `increment` attribute, which
            determines the discreteness of the rounding.
        """
        inc = self.increment
        size = np.ceil(num / self.sample_rate / inc) * self.sample_rate * inc
        return size

    def getslice(self, sindex):
        instance = copy(self)
        instance.table = self.table[sindex]
        return instance

    def id_from_param(self, param_tuple):
        """Get the index of this template based on its param tuple

        Parameters
        ----------
        param_tuple : tuple
            Tuple of the parameters which uniquely identify this template

        Returns
        --------
        index : int
            The ordered index that this template has in the template bank.
        """
        return self.param_lookup[param_tuple]

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.getslice(index)

        return self.get_template(index)

    def freq_resolution_for_template(self, index):
        """Compute the correct resolution for a frequency series that contains
        a given template in the bank.
        """
        from pycbc.waveform.waveform import props

        time_duration = self.minimum_buffer
        time_duration += 0.5
        params = props(self.table[index])
        params.pop('approximant')
        approximant = self.approximant(index)
        waveform_duration = pycbc.waveform.get_waveform_filter_length_in_time(
            approximant, **params
        )
        if waveform_duration is None:
            raise RuntimeError(
                'Template waveform {approximant} not recognized!'
            )
        time_duration += waveform_duration
        td_samples = self.round_up(time_duration * self.sample_rate)
        return self.sample_rate / float(td_samples)

    def get_template(self, index, delta_f=None):
        """Calculate and return the frequency-domain waveform for the template
        with the given index. The frequency resolution can optionally be given.

        Parameters
        ----------
        index: int
            Index of the template in the bank.
        delta_f: float, optional
            Resolution of the resulting frequency series. If not given, it is
            calculated from the time duration of the template.

        Returns
        -------
        htilde: FrequencySeries
            Template waveform in the frequency domain.
        """
        approximant = self.approximant(index)
        f_end = self.end_frequency(index)
        flow = self.table[index].f_lower

        if delta_f is None:
            delta_f = self.freq_resolution_for_template(index)

        flen = int(self.sample_rate / (2 * delta_f) + 1)

        if f_end is None or f_end >= (flen * delta_f):
            f_end = (flen - 1) * delta_f

        logging.info(
            "Generating %s, duration %s s, index %i, starting from %s Hz",
            approximant,
            1.0 / delta_f,
            index,
            flow
        )

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        htilde = pycbc.waveform.get_waveform_filter(
            zeros(flen, dtype=np.complex64), self.table[index],
            approximant=approximant, f_lower=flow, f_final=f_end,
            delta_f=delta_f, delta_t=1.0 / self.sample_rate, distance=distance,
            **self.extra_args)

        # If available, record the total duration (which may
        # include ringdown) and the duration up to merger since they will be
        # erased by the type conversion below.
        ttotal = template_duration = -1
        time_offset = None
        if hasattr(htilde, 'length_in_time'):
            ttotal = htilde.length_in_time
        if hasattr(htilde, 'chirp_length'):
            template_duration = htilde.chirp_length
        if hasattr(htilde, 'time_offset'):
            time_offset = htilde.time_offset

        self.table[index].template_duration = template_duration

        htilde = htilde.astype(np.complex64)
        htilde.f_lower = flow
        htilde.min_f_lower = self.min_f_lower
        htilde.end_idx = int(f_end / htilde.delta_f)
        htilde.params = self.table[index]
        htilde.chirp_length = template_duration
        htilde.length_in_time = ttotal
        htilde.approximant = approximant
        htilde.end_frequency = f_end

        if time_offset:
            htilde.time_offset = time_offset

        # Add sigmasq as a method of this instance
        htilde.sigmasq = types.MethodType(sigma_cached, htilde)

        htilde.id = self.id_from_param((htilde.params.mass1,
                                        htilde.params.mass2,
                                        htilde.params.spin1z,
                                        htilde.params.spin2z))
        return htilde


class FilterBank(TemplateBank):
    def __init__(self, filename, filter_length, delta_f, dtype,
                 out=None, max_template_length=None,
                 approximant=None, parameters=None,
                 enable_compressed_waveforms=True,
                 low_frequency_cutoff=None,
                 waveform_decompression_method=None,
                 **kwds):
        self.out = out
        self.dtype = dtype
        self.f_lower = low_frequency_cutoff
        self.filename = filename
        self.delta_f = delta_f
        self.N = (filter_length - 1 ) * 2
        self.delta_t = 1.0 / (self.N * self.delta_f)
        self.filter_length = filter_length
        self.max_template_length = max_template_length
        self.enable_compressed_waveforms = enable_compressed_waveforms
        self.waveform_decompression_method = waveform_decompression_method

        super(FilterBank, self).__init__(filename, approximant=approximant,
            parameters=parameters, **kwds)
        self.ensure_standard_filter_columns(low_frequency_cutoff=low_frequency_cutoff)

    def get_decompressed_waveform(self, tempout, index, f_lower=None,
                                  approximant=None, df=None):
        """Returns a frequency domain decompressed waveform for the template
        in the bank corresponding to the index taken in as an argument. The
        decompressed waveform is obtained by interpolating in frequency space,
        the amplitude and phase points for the compressed template that are
        read in from the bank."""

        from pycbc.waveform.waveform import props
        from pycbc.waveform import get_waveform_filter_length_in_time

        # Get the template hash corresponding to the template index taken in as argument
        tmplt_hash = self.table.template_hash[index]

        # Read the compressed waveform from the bank file
        compressed_waveform = pycbc.waveform.compress.CompressedWaveform.from_hdf(
                                self.filehandler, tmplt_hash,
                                load_now=True)

        # Get the interpolation method to be used to decompress the waveform
        if self.waveform_decompression_method is not None :
            decompression_method = self.waveform_decompression_method
        else :
            decompression_method = compressed_waveform.interpolation
        logging.info("Decompressing waveform using %s", decompression_method)

        if df is not None :
            delta_f = df
        else :
            delta_f = self.delta_f

        # Create memory space for writing the decompressed waveform
        decomp_scratch = FrequencySeries(tempout[0:self.filter_length], delta_f=delta_f, copy=False)

        # Get the decompressed waveform
        hdecomp = compressed_waveform.decompress(out=decomp_scratch, f_lower=f_lower, interpolation=decompression_method)
        p = props(self.table[index])
        p.pop('approximant')
        try:
            tmpltdur = self.table[index].template_duration
        except AttributeError:
            tmpltdur = None
        if tmpltdur is None or tmpltdur==0.0 :
            tmpltdur = get_waveform_filter_length_in_time(approximant, **p)
        hdecomp.chirp_length = tmpltdur
        hdecomp.length_in_time = hdecomp.chirp_length
        return hdecomp

    def generate_with_delta_f_and_max_freq(self, t_num, max_freq, delta_f,
                                           low_frequency_cutoff=None,
                                           cached_mem=None):
        """Generate the template with index t_num using custom length."""
        approximant = self.approximant(t_num)
        # Don't want to use INTERP waveforms in here
        if approximant.endswith('_INTERP'):
            approximant = approximant.replace('_INTERP', '')
        # Using SPAtmplt here is bad as the stored cbrt and logv get
        # recalculated as we change delta_f values. Fall back to TaylorF2
        # in lalsimulation.
        if approximant == 'SPAtmplt':
            approximant = 'TaylorF2'
        if cached_mem is None:
            wav_len = int(max_freq / delta_f) + 1
            cached_mem = zeros(wav_len, dtype=np.complex64)

        full_calculate_waveform = True
        if (self.has_compressed_waveforms and self.enable_compressed_waveforms):
            try:
                htilde = self.get_decompressed_waveform(
                    tempout,
                    index,
                    f_lower=f_low,
                    approximant=approximant,
                    df=None
                )
                full_calculate_waveform = False
            except KeyError:
                # This is the error caused when the compressed waveform is
                # not in the bank.
                warnings.warn(
                    "self.get_decompressed_waveform has raised a KeyError. "
                    "This may be as the compressed waveform has not been "
                    "generated for this approximant, but it could indicate "
                    "a more serious issue. Approximant: %s" % approximant
                )

        if full_calculate_waveform:
            htilde = pycbc.waveform.get_waveform_filter(
                cached_mem, self.table[t_num], approximant=approximant,
                f_lower=low_frequency_cutoff, f_final=max_freq, delta_f=delta_f,
                distance=1./DYN_RANGE_FAC, delta_t=1./(2.*max_freq)
            )

        return htilde

    def __getitem__(self, index):
        # Make new memory for templates if we aren't given output memory
        if self.out is None:
            tempout = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempout = self.out

        approximant = self.approximant(index)
        f_end = self.end_frequency(index)
        if f_end is None or f_end >= (self.filter_length * self.delta_f):
            f_end = (self.filter_length-1) * self.delta_f

        # Find the start frequency, if variable
        f_low = find_variable_start_frequency(approximant,
                                              self.table[index],
                                              self.f_lower,
                                              self.max_template_length)
        logging.info('%s: generating %s from %s Hz' % (index, approximant, f_low))

        # Clear the storage memory
        poke  = tempout.data # pylint:disable=unused-variable
        tempout.clear()

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        full_calculate_waveform = True
        if (self.has_compressed_waveforms and self.enable_compressed_waveforms):
            try:
                htilde = self.get_decompressed_waveform(
                    tempout,
                    index,
                    f_lower=f_low,
                    approximant=approximant,
                    df=None
                )
                full_calculate_waveform = False
            except KeyError:
                # This is the error caused when the compressed waveform is
                # not in the bank.
                warnings.warn(
                    "self.get_decompressed_waveform has raised a KeyError. "
                    "This may be as the compressed waveform has not been "
                    "generated for this approximant, but it could indicate "
                    "a more serious issue. Approximant: %s" % approximant
                )

        if full_calculate_waveform:
            htilde = pycbc.waveform.get_waveform_filter(
                tempout[0:self.filter_length], self.table[index],
                approximant=approximant, f_lower=f_low, f_final=f_end,
                delta_f=self.delta_f, delta_t=self.delta_t, distance=distance,
                **self.extra_args,
            )

        # If available, record the total duration (which may
        # include ringdown) and the duration up to merger since they will be
        # erased by the type conversion below.
        ttotal = template_duration = None
        if hasattr(htilde, 'length_in_time'):
            ttotal = htilde.length_in_time
        if hasattr(htilde, 'chirp_length'):
            template_duration = htilde.chirp_length

        self.table[index].template_duration = template_duration

        htilde = htilde.astype(self.dtype)
        htilde.f_lower = f_low
        htilde.min_f_lower = self.min_f_lower
        htilde.end_idx = int(f_end / htilde.delta_f)
        htilde.params = self.table[index]
        htilde.chirp_length = template_duration
        htilde.length_in_time = ttotal
        htilde.approximant = approximant
        htilde.end_frequency = f_end

        # Add sigmasq as a method of this instance
        htilde.sigmasq = types.MethodType(sigma_cached, htilde)
        htilde._sigmasq = {}
        return htilde


def find_variable_start_frequency(approximant, parameters, f_start, max_length,
                                  delta_f = 1):
    """ Find a frequency value above the starting frequency that results in a
    waveform shorter than max_length.
    """
    if (f_start is None):
        f = parameters.f_lower
    elif (max_length is not None):
        l = max_length + 1
        f = f_start - delta_f
        while l > max_length:
            f += delta_f
            l = pycbc.waveform.get_waveform_filter_length_in_time(approximant,
                                                          parameters, f_lower=f)
    else :
        f = f_start
    return f


class FilterBankSkyMax(TemplateBank):
    def __init__(self, filename, filter_length, delta_f,
                 dtype, out_plus=None, out_cross=None,
                 max_template_length=None, parameters=None,
                 low_frequency_cutoff=None, **kwds):
        self.out_plus = out_plus
        self.out_cross = out_cross
        self.dtype = dtype
        self.f_lower = low_frequency_cutoff
        self.filename = filename
        self.delta_f = delta_f
        self.N = (filter_length - 1 ) * 2
        self.delta_t = 1.0 / (self.N * self.delta_f)
        self.filter_length = filter_length
        self.max_template_length = max_template_length

        super(FilterBankSkyMax, self).__init__(filename, parameters=parameters,
              **kwds)

        self.ensure_standard_filter_columns(low_frequency_cutoff=low_frequency_cutoff)

    def __getitem__(self, index):
        # Make new memory for templates if we aren't given output memory
        if self.out_plus is None:
            tempoutplus = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutplus = self.out_plus
        if self.out_cross is None:
            tempoutcross = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutcross = self.out_cross

        approximant = self.approximant(index)

        # Get the end of the waveform if applicable (only for SPAtmplt atm)
        f_end = self.end_frequency(index)
        if f_end is None or f_end >= (self.filter_length * self.delta_f):
            f_end = (self.filter_length-1) * self.delta_f

        # Find the start frequency, if variable
        f_low = find_variable_start_frequency(approximant,
                                              self.table[index],
                                              self.f_lower,
                                              self.max_template_length)
        logging.info('%s: generating %s from %s Hz', index, approximant, f_low)

        # What does this do???
        poke1 = tempoutplus.data # pylint:disable=unused-variable
        poke2 = tempoutcross.data # pylint:disable=unused-variable

        # Clear the storage memory
        tempoutplus.clear()
        tempoutcross.clear()

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        hplus, hcross = pycbc.waveform.get_two_pol_waveform_filter(
            tempoutplus[0:self.filter_length],
            tempoutcross[0:self.filter_length], self.table[index],
            approximant=approximant, f_lower=f_low,
            f_final=f_end, delta_f=self.delta_f, delta_t=self.delta_t,
            distance=distance, **self.extra_args)

        if hasattr(hplus, 'chirp_length') and hplus.chirp_length is not None:
            self.table[index].template_duration = hplus.chirp_length

        hplus = hplus.astype(self.dtype)
        hcross = hcross.astype(self.dtype)
        hplus.f_lower = f_low
        hcross.f_lower = f_low
        hplus.min_f_lower = self.min_f_lower
        hcross.min_f_lower = self.min_f_lower
        hplus.end_frequency = f_end
        hcross.end_frequency = f_end
        hplus.end_idx = int(hplus.end_frequency / hplus.delta_f)
        hcross.end_idx = int(hplus.end_frequency / hplus.delta_f)
        hplus.params = self.table[index]
        hcross.params = self.table[index]
        hplus.approximant = approximant
        hcross.approximant = approximant

        # Add sigmasq as a method of this instance
        hplus.sigmasq = types.MethodType(sigma_cached, hplus)
        hplus._sigmasq = {}
        hcross.sigmasq = types.MethodType(sigma_cached, hcross)
        hcross._sigmasq = {}

        return hplus, hcross


class FilterBankTHA(TemplateBank):
    def __init__(self, filename, filter_length, delta_f,
                 dtype, out_comp1=None, out_comp2=None,
                 out_comp3=None, out_comp4=None, out_comp5=None,
                 max_template_length=None, parameters=None,
                 low_frequency_cutoff=None, max_num_comps=None, **kwds):
        self.out_comp1 = out_comp1
        self.out_comp2 = out_comp2
        self.out_comp3 = out_comp3
        self.out_comp4 = out_comp4
        self.out_comp5 = out_comp5
        self.dtype = dtype
        self.f_lower = low_frequency_cutoff
        self.filename = filename
        self.delta_f = delta_f
        self.N = (filter_length - 1 ) * 2
        self.delta_t = 1.0 / (self.N * self.delta_f)
        self.filter_length = filter_length
        self.max_template_length = max_template_length
        if max_num_comps is None:
            max_num_comps = 5
        self.max_num_comps = max_num_comps

        super(FilterBankTHA, self).__init__(filename, parameters=parameters,
              **kwds)

        self.ensure_standard_filter_columns(low_frequency_cutoff=low_frequency_cutoff)

        self.last_index = None

    def get(self, index, psd):
        if index == self.last_index:
             need_new_comps = False
        else:
             need_new_comps = True
        self.last_index = index
        # Make new memory for templates if we aren't given output memory
        if self.out_comp1 is None:
            tempoutcomp1 = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutcomp1 = self.out_comp1
        if self.out_comp2 is None:
            tempoutcomp2 = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutcomp2 = self.out_comp2
        if self.out_comp3 is None:
            tempoutcomp3 = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutcomp3 = self.out_comp3
        if self.out_comp4 is None:
            tempoutcomp4 = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutcomp4 = self.out_comp4
        if self.out_comp5 is None:
            tempoutcomp5 = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutcomp5 = self.out_comp5

        approximant = self.approximant(index)
        assert(approximant in ["IMRPhenomPv2", "IMRPhenomXP"])
        if approximant == "IMRPhenomPv2":
            _approx_cls = PhenomPv2Template
        else:
            _approx_cls = PhenomXPTemplate

        # Get the end of the waveform if applicable (only for SPAtmplt atm)
        f_end = self.end_frequency(index)
        if f_end is None or f_end >= (self.filter_length * self.delta_f):
            f_end = (self.filter_length-1) * self.delta_f

        # Find the start frequency, if variable
        f_low = self.f_lower
        logging.info('%s: generating %s from %s Hz', index, approximant, f_low)

        # What does this do???
        poke1 = tempoutcomp1.data # pylint:disable=unused-variable
        poke2 = tempoutcomp2.data # pylint:disable=unused-variable
        poke3 = tempoutcomp3.data # pylint:disable=unused-variable
        poke4 = tempoutcomp4.data # pylint:disable=unused-variable
        poke5 = tempoutcomp5.data # pylint:disable=unused-variable

        # Clear the storage memory
        tempoutcomp1.clear()
        tempoutcomp2.clear()
        tempoutcomp3.clear()
        tempoutcomp4.clear()
        tempoutcomp5.clear()

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        if need_new_comps:
            curr_tmp = _approx_cls(self.table[index], 1./self.delta_t,
                                      self.f_lower)
            self.curr_tmp = curr_tmp
        else:
            curr_tmp = self.curr_tmp
        num_comps = min(self.table["num_comps"][index], self.max_num_comps)
        hcomps = curr_tmp.get_whitened_normalized_comps(self.delta_f, psd,
                                                        num_comps=num_comps)

        tempoutcomp1.data[:] = hcomps[0].data[:]
        if hasattr(hcomps[0], 'chirp_length') and \
                hcomps[0].chirp_length is not None:
            self.table[index].template_duration = hcomps[0].chirp_length
        hcomps[0] = hcomps[0].astype(self.dtype)
        hcomps[0].f_lower = f_low
        hcomps[0].min_f_lower = f_low
        hcomps[0].end_frequency = f_end
        hcomps[0].end_idx = int(hcomps[0].end_frequency / hcomps[0].delta_f)
        hcomps[0].params = self.table[index]
        hcomps[0].approximant = approximant

        if hcomps[1] is not None:
            tempoutcomp2.data[:] = hcomps[1].data[:]
            hcomps[1] = hcomps[1].astype(self.dtype)
            hcomps[1].f_lower = f_low
            hcomps[1].min_f_lower = f_low
            hcomps[1].end_frequency = f_end
            hcomps[1].end_idx = int(hcomps[0].end_frequency / hcomps[0].delta_f)
            hcomps[1].params = self.table[index]
            hcomps[1].approximant = approximant

        if hcomps[2] is not None:
            tempoutcomp3.data[:] = hcomps[2].data[:]
            hcomps[2] = hcomps[2].astype(self.dtype)
            hcomps[2].f_lower = f_low
            hcomps[2].min_f_lower = f_low
            hcomps[2].end_frequency = f_end
            hcomps[2].end_idx = int(hcomps[0].end_frequency / hcomps[0].delta_f)
            hcomps[2].params = self.table[index]
            hcomps[2].approximant = approximant

        if hcomps[3] is not None:
            tempoutcomp4.data[:] = hcomps[3].data[:]
            hcomps[3] = hcomps[3].astype(self.dtype)
            hcomps[3].f_lower = f_low
            hcomps[3].min_f_lower = f_low
            hcomps[3].end_frequency = f_end
            hcomps[3].end_idx = int(hcomps[0].end_frequency / hcomps[0].delta_f)
            hcomps[3].params = self.table[index]
            hcomps[3].approximant = approximant

        if hcomps[4] is not None:
            tempoutcomp5.data[:] = hcomps[4].data[:]
            hcomps[4] = hcomps[4].astype(self.dtype)
            hcomps[4].f_lower = f_low
            hcomps[4].min_f_lower = f_low
            hcomps[4].end_frequency = f_end
            hcomps[4].end_idx = int(hcomps[0].end_frequency / hcomps[0].delta_f)
            hcomps[4].params = self.table[index]
            hcomps[4].approximant = approximant

        return hcomps

def _dpsi(theta_jn, phi_jl, beta):
    """Calculate the difference between the polarization with respect to the
    total angular momentum and the polarization with respect to the orbital
    angular momentum using code from
    https://git.ligo.org/lscsoft/pesummary/-/blob/master/pesummary/gw/conversions/angles.py#L13
    """
    if theta_jn == 0:
        return -1. * phi_jl
    n = np.array([np.sin(theta_jn), 0, np.cos(theta_jn)])
    j = np.array([0, 0, 1])
    l = np.array([
        np.sin(beta) * np.sin(phi_jl), np.sin(beta) * np.cos(phi_jl), np.cos(beta)
    ])
    p_j = np.cross(n, j)
    p_j /= np.linalg.norm(p_j)
    p_l = np.cross(n, l)
    p_l /= np.linalg.norm(p_l)
    cosine = np.inner(p_j, p_l)
    sine = np.inner(n, np.cross(p_j, p_l))
    dpsi = np.pi / 2 + np.sign(sine) * np.arccos(cosine)
    return dpsi

def _dphi(theta_jn, phi_jl, beta):
    """Calculate the difference in the phase angle between J-aligned
    and L-aligned frames using code from
    https://git.ligo.org/lscsoft/pesummary/-/blob/master/pesummary/gw/conversions/angles.py#L36

    Parameters
    ----------
    theta_jn: np.ndarray
        the angle between J and line of sight
    phi_jl: np.ndarray
        the precession phase
    beta: np.ndarray
        the opening angle (angle between J and L)
    """
    theta_jn = np.array([theta_jn])
    phi_jl = np.array([phi_jl])
    beta = np.array([beta])
    n = np.column_stack(
        [np.repeat([0], len(theta_jn)), np.sin(theta_jn), np.cos(theta_jn)]
    )
    l = np.column_stack(
        [
            np.sin(beta) * np.cos(phi_jl), np.sin(beta) * np.sin(phi_jl),
            np.cos(beta)
        ]
    )
    cosi = [np.inner(nn, ll) for nn, ll in zip(n, l)]
    inc = np.arccos(cosi)
    sign = np.sign(np.cos(theta_jn) - (np.cos(beta) * np.cos(inc)))
    cos_d = np.cos(phi_jl) * np.sin(theta_jn) / np.sin(inc)
    inds = np.logical_or(cos_d < -1, cos_d > 1)
    cos_d[inds] = np.sign(cos_d[inds]) * 1.
    dphi = -1. * sign * np.arccos(cos_d)
    return dphi[0]

def compute_beta(tmplt):
    """ Calculate beta (thetaJL) using code from
    https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_inspiral_8c_source.html#l06105
    """
    m1 = tmplt.mass1
    m2 = tmplt.mass2
    s1x = tmplt.spin1x
    s1y = tmplt.spin1y
    s1z = tmplt.spin1z
    s2x = tmplt.spin2x
    s2y = tmplt.spin2y
    s2z = tmplt.spin2z
    flow = tmplt.flow

    eta = m1 * m2 / (m1 + m2) / (m1 + m2);
    v0 = ((m1 + m2) * lal.MTSUN_SI * lal.PI * flow) ** (1. / 3.)

    lmag = (m1 + m2) * (m1 + m2) * eta / v0
    lmag *= (1.0 + v0 * v0 * (1.5 + eta / 6.))

    s1x = m1 * m1 * s1x
    s1y = m1 * m1 * s1y
    s1z = m1 * m1 * s1z
    s2x = m2 * m2 * s2x
    s2y = m2 * m2 * s2y
    s2z = m2 * m2 * s2z
    jx = s1x + s2x
    jy = s1y + s2y
    jz = lmag + s1z + s2z

    jnorm = (jx * jx + jy * jy + jz * jz) ** (1. / 2.)
    jhatz = jz / jnorm

    return np.arccos(jhatz)

class _PhenomTemplate():

    def __init__(self, template_params, sample_rate, f_lower):
        from pycbc.pnutils import get_imr_duration
        self.flow = float(f_lower)
        self.f_final = float(sample_rate / 2.)

        self.mass1 = float(template_params.mass1)
        self.mass2 = float(template_params.mass2)
        self.spin1x = float(template_params.spin1x)
        self.spin1y = float(template_params.spin1y)
        self.spin1z = float(template_params.spin1z)
        self.spin2x = float(template_params.spin2x)
        self.spin2y = float(template_params.spin2y)
        self.spin2z = float(template_params.spin2z)

        self.theta = float(template_params.latitude)
        self.phi = float(template_params.longitude)
        self.iota = float(template_params.inclination)
        self.psi = float(template_params.polarization)
        self.orb_phase = float(template_params.orbital_phase)
        self.fref = self.flow

        self.duration = get_imr_duration(self.mass1, self.mass2, self.spin1z,
                                         self.spin2z, self.flow, "IMRPhenomD")

        outs = self._model_parameters_from_source_frame(
            self.mass1*lal.MSUN_SI,
            self.mass2*lal.MSUN_SI,
            self.flow,
            self.orb_phase,
            self.iota,
            self.spin1x,
            self.spin1y,
            self.spin1z,
            self.spin2x,
            self.spin2y,
            self.spin2z,
        )
        chi1_l, chi2_l, chip, thetaJN, alpha0, phi_aligned, zeta_polariz = outs

        self.chi1_l = float(chi1_l)
        self.chi2_l = float(chi2_l)
        self.chip = float(chip)
        self.thetaJN = float(thetaJN)
        self.alpha0 = float(alpha0)
        self.phi0 = float(phi_aligned)
        # This is a correction on psi, currently unused
        self.psi_corr = zeta_polariz
        self.beta = compute_beta(self)

        self.comps = {}
        self.has_comps = False

    def gen_hp_hc(self, thetaJN, alpha0, phi0, df, f_final):
        # calculate cartesian spins for waveform generator
        a1 = np.sqrt(
            np.sum(np.square([self.spin1x, self.spin1y, self.spin1z]))
        )
        a2 = np.sqrt(
            np.sum(np.square([self.spin2x, self.spin2y, self.spin2z]))
        )
        phi1 = np.fmod(
            2 * np.pi + np.arctan2(self.spin1y, self.spin1x),
            2 * np.pi
        )
        phi2 = np.fmod(
            2 * np.pi + np.arctan2(self.spin2y, self.spin2x),
            2 * np.pi
        )
        phi12 = phi2 - phi1
        if phi12 < 0:
            phi12 += 2 * np.pi
        # prevent failure when s1z = a1 = 0
        if not self.spin1z and not a1:
            tilt1 = 0.
        else:
            tilt1 = np.arccos(self.spin1z / a1)
        if not self.spin2z and not a2:
            tilt2 = 0.
        else:
            tilt2 = np.arccos(self.spin2z / a2)
        iota, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = \
            lalsim.SimInspiralTransformPrecessingNewInitialConditions(
                thetaJN, alpha0, tilt1, tilt2, phi12, a1, a2,
                self.mass1*lal.MSUN_SI, self.mass2*lal.MSUN_SI, self.fref, phi0
            )
        # generate hp, hc for given orientation with lalsimulation
        return lalsim.SimInspiralChooseFDWaveform(
            self.mass1*lal.MSUN_SI, self.mass2*lal.MSUN_SI, spin1x, spin1y,
            spin1z, spin2x, spin2y, spin2z, 1.e6*lal.PC_SI, iota, phi0,
            0, 0, 0, df, self.flow, f_final, self.fref, lal.CreateDict(),
            lalsim.GetApproximantFromString(self.approximant)
        )

    def gen_harmonics_comp(self, thetaJN, alpha0, phi0, psi, df, f_final):
        # generate hp, hc
        hp, hc = self.gen_hp_hc(thetaJN, alpha0, phi0, df, f_final)
        # 1908.05707 defines psi in J-aligned frame. Need to rotate to
        # L-aligned frame and multiply by w+, wx
        dpsi = _dpsi(thetaJN, alpha0, self.beta)
        fp = np.cos(2 * (psi - dpsi))
        fc = -1. * np.sin(2 * (psi - dpsi))
        h = (fp * hp.data.data[:] + fc * hc.data.data[:])
        # 1908.05707 defines phi in J-aligned frame. Need to rotate to
        # L-aligned frame
        h *= np.exp(2j * _dphi(thetaJN, alpha0, self.beta))
        # create LAL frequency array and return precessing harmonic
        new = lal.CreateCOMPLEX16FrequencySeries(
            "", lal.LIGOTimeGPS(hp.epoch), 0, df, lal.SecondUnit, len(h)
        )
        new.data.data[:] = h[:]
        return new

    def get_interpolated_harmonic_comp(
        self, thetaJN, alpha0, phi0, psi, df, f_final
    ):
        # Waveform generation is a problem.
        # Compressed waveforms would be an option, but the file size will be
        # challenging. I'm trying out the idea of "INTERP" waveforms instead,
        # where we generate waveforms at much larger frequency spacing and then
        # upsample. (Technically these don't actually interpolate, but it's a
        # good way to explain the idea of what it's doing ..).
        from pycbc.filter import interpolate_complex_frequency

        def rulog2(val):
            return 2.0 ** np.ceil(np.log2(float(val)))


        # FIXME: THIS ONE IS IMPORTANT!!
        extra_padding = 5
        df_min = 1.0 / rulog2(self.duration + extra_padding)

        small = self.gen_harmonics_comp(
            thetaJN, alpha0, phi0, psi, df_min, f_final
        )

        offset = int(extra_padding * (small.data.length-1)*2 * small.deltaF)

        small = FrequencySeries(small.data.data[:], delta_f=small.deltaF,
                                epoch=small.epoch)

        large = interpolate_complex_frequency(small, df, zeros_offset=offset,
                                              side='left')
        return large

    def compute_waveform_five_comps(self, df, f_final, num_comps=5, interp=True):
        if interp:
            def _func(*args, **kwargs):
                return self.get_interpolated_harmonic_comp(*args, **kwargs)
        else:
            def _func(*args, **kwargs):
                h = self.gen_harmonics_comp(*args, **kwargs)
                return FrequencySeries(h.data.data[:], delta_f=h.deltaF, epoch=h.epoch)
        # calculate 5 harmonic decomposition as defined in 1908.05707
        hgen1a = _func(0., 0., 0., 0., df, f_final)
        hgen1b = _func(0., 0., np.pi/4., np.pi/4, df, f_final)
        tmp = hgen1a - hgen1b
        hgen1b = (hgen1a + hgen1b) / 2.
        hgen1a = tmp / 2.
        h1 = hgen1a

        if num_comps == 1:
            return h1, None, None, None, None


        # These are both negative w.r.t 1908.05707
        hgen2a = _func(np.pi/2., 0., np.pi/4., np.pi/4, df, f_final)
        hgen2b = _func(np.pi/2., np.pi/2., 0., np.pi/4, df, f_final)
        tmp = hgen2a + hgen2b
        hgen2b = -0.25 * (hgen2a - hgen2b)
        hgen2a = -0.25 * tmp
        h2 = hgen2a

        if num_comps == 2:
            return h1, h2, None, None, None

        hgen3a = _func(np.pi/2., 0., 0., 0., df, f_final)
        hgen3b = _func(np.pi/2., np.pi/2., 0., 0., df, f_final)
        hgen3a = 1./6. * (hgen3a + hgen3b)
        h3 = hgen3a

        if num_comps == 3:
            return h1, h2, h3, None, None

        h4 = hgen2b

        if num_comps == 4:
            return h1, h2, h3, h4, None

        h5 = hgen1b

        return h1, h2, h3, h4, h5

    def wn_cython(self, hs, ASD, flen, df, kmin, kmax):
        from pycbc.types.array_cpu import whiten_and_normalize_three
        from pycbc.types.array_cpu import whiten_and_normalize_two
        from pycbc.types.array_cpu import whiten_and_normalize_one
        if len(hs) == 3:
            whiten_and_normalize_three(hs[0].data, hs[1].data, hs[2].data,
                                       ASD.data, flen, df, kmin, kmax)
        elif len(hs) == 2:
            whiten_and_normalize_two(hs[0].data, hs[1].data,
                                     ASD.data, flen, df, kmin, kmax)
        elif len(hs) == 1:
            whiten_and_normalize_one(hs[0].data, ASD.data, flen, df, kmin, kmax)


    def whiten_and_normalize(self, hs, ASD, flen, df, kmin, kmax):
        if len(hs) in [3,2,1]:
            self.wn_cython(hs, ASD, flen, df, kmin, kmax)
            return
        raise NotImplementedError()

    def orthogonalize(self, hs, df):
        from pycbc.types.array_cpu import tha_orthogonalize_vecs
        for i in range(len(hs)):
            for j in range(i + 1, len(hs)):
                tha_orthogonalize_vecs(hs[i].data, hs[j].data, df, len(hs[i]))

    def get_whitened_normalized_comps(self, df, psd, num_comps=5):
        """
        Return a FrequencySeries of h+ and hx, whitened by the
        given ASD and normalized. The waveform is not zero-padded to
        match the length of the ASD, so its normalization depends on
        its own length.
        """
        # Generate a new wf
        if not self.has_comps:
            h1, h2, h3, h4, h5 = self.compute_waveform_five_comps(df, self.f_final, num_comps=num_comps)
            self.h1 = h1
            self.h2 = h2
            self.h3 = h3
            self.h4 = h4
            self.h5 = h5
            self.has_comps = True
        else:
            h1 = self.h1
            h2 = self.h2
            h3 = self.h3
            h4 = self.h4
            h5 = self.h5
        hs = [h1.copy()]
        if h2 is not None:
            hs += [h2.copy()]
        if h3 is not None:
            hs += [h3.copy()]
        if h4 is not None:
            hs += [h4.copy()]
        if h5 is not None:
            hs += [h5.copy()]

        flen = len(h1)
        ASD = psd ** 0.5
        kmin = int(self.flow / df)
        kmax = int(self.f_final / df)

        if len(h1) > len(ASD):
            err_msg = "waveform has length greater than ASD; cannot whiten"
            raise ValueError(err_msg)

        orthogonal = []

        # Whiten and Normalize
        self.whiten_and_normalize(hs, ASD, flen, df, kmin, kmax)

        # Orthogonalize
        self.orthogonalize(hs, df)

        for i in range(len(hs)):
            orthogonal += [hs[i]]

        orthogonal += [None] * (5 - len(hs))

        return orthogonal

class PhenomPv2Template(_PhenomTemplate):
    approximant = "IMRPhenomPv2"

    def _model_parameters_from_source_frame(self, *args):
        return lalsim.SimIMRPhenomPCalculateModelParametersFromSourceFrame(
            *args, lalsim.IMRPhenomPv2_V
        )

class PhenomXPTemplate(_PhenomTemplate):
    approximant = "IMRPhenomXP"

    def _model_parameters_from_source_frame(self, *args):
        return lalsim.SimIMRPhenomXPCalculateModelParametersFromSourceFrame(
            *args, None
        )

__all__ = ('sigma_cached', 'boolargs_from_apprxstr', 'add_approximant_arg',
           'parse_approximant_arg', 'tuple_to_hash', 'TemplateBank',
           'LiveFilterBank', 'FilterBank', 'find_variable_start_frequency',
           'FilterBankSkyMax', 'FilterBankTHA')

