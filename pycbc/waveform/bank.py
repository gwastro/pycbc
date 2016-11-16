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
import numpy
import logging
import os.path
import h5py
from copy import copy
import numpy as np
from glue.ligolw import ligolw, table, lsctables, utils as ligolw_utils
import pycbc.waveform
import pycbc.pnutils
import pycbc.waveform.compress
from pycbc import DYN_RANGE_FAC
from pycbc.types import zeros
import pycbc.io

def sigma_cached(self, psd):
    """ Cache sigma calculate for use in tandem with the FilterBank class
    """
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
                psd.sigmasq_vec[self.approximant] = pycbc.waveform.get_waveform_filter_norm(
                     self.approximant, psd, len(psd), psd.delta_f, self.f_lower)

            if not hasattr(self, 'sigma_scale'):
                # Get an amplitude normalization (mass dependant constant norm)
                amp_norm = pycbc.waveform.get_template_amplitude_norm(
                                     self.params, approximant=self.approximant)
                amp_norm = 1 if amp_norm is None else amp_norm
                self.sigma_scale = (DYN_RANGE_FAC * amp_norm) ** 2.0


            self._sigmasq[key] = psd.sigmasq_vec[self.approximant][self.end_idx] * self.sigma_scale

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
                psd.invsqrt = 1.0 / psd[self.sslice]

            return self.sigma_view.inner(psd.invsqrt)
    return self._sigmasq[key]

# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(LIGOLWContentHandler)

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
    load_compressed : {True, bool}
        If compressed waveforms are present in the file, load them. Only works
        for hdf files.
    load_compressed_now : {False, bool}
        If compressed waveforms are present in the file, load all of them into
        memory now. Otherwise, they will be loaded into memory on their first
        use. Note: if not `load_compressed_now`, the filehandler to the hdf
        file must be kept open.
    \**kwds :
        Any additional keyword arguments are stored to the `extra_args`
        attribute.


    Attributes
    ----------
    table : WaveformArray
        An instance of a WaveformArray containing all of the information about
        the parameters of the bank.
    compressed_waveforms : {None, dict}
        If compressed waveforms are present in the (hdf) file, a dictionary
        of those waveforms. Keys are the template hashes, values are
        `CompressedWaveform` instances.
    parameters : tuple
        The parameters loaded from the input file. Same as `table.fieldnames`.
    indoc : {None, xmldoc}
        If an xml file was provided, an in-memory representation of the xml.
        Otherwise, None.
    filehandler : {None, h5py.File}
        If an hdf file was provided, the file handler pointing to the hdf file
        (left open after initialization). Otherwise, None.
    extra_args : {None, dict}
        Any extra keyword arguments that were provided on initialization.
    """
    def __init__(self, filename, approximant=None, parameters=None,
            load_compressed=True, load_compressed_now=False,
            **kwds):
        ext = os.path.basename(filename)
        self.compressed_waveforms = None
        if ext.endswith('.xml') or ext.endswith('.xml.gz') or ext.endswith('.xmlgz'):
            self.filehandler = None
            self.indoc = ligolw_utils.load_filename(
                filename, False, contenthandler=LIGOLWContentHandler)
            self.table = table.get_table(
                self.indoc, lsctables.SnglInspiralTable.tableName)
            self.table = pycbc.io.WaveformArray.from_ligolw_table(self.table,
                columns=parameters)

            # inclination stored in xml alpha3 column
            names = list(self.table.dtype.names)
            names = tuple([n if n != 'alpha3' else 'inclination' for n in names]) 

            # low frequency cutoff in xml alpha6 column
            names = tuple([n if n!= 'alpha6' else 'f_lower' for n in names])
            self.table.dtype.names = names

        elif ext.endswith('hdf'):
            self.indoc = None
            f = h5py.File(filename, 'r')
            self.filehandler = f
            try:
                fileparams = map(str, f.attrs['parameters'])
            except KeyError:
                # just assume all of the top-level groups are the parameters
                fileparams = map(str, f.keys())
                logging.info("WARNING: no parameters attribute found. "
                    "Assuming that %s " %(', '.join(fileparams)) +
                    "are the parameters.")
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
                data[str(key)] = f[key][:]
                dtype.append((str(key), data[key].dtype))
            num = f[fileparams[0]].size
            self.table = pycbc.io.WaveformArray(num, dtype=dtype)
            for key in data:
                self.table[key] = data[key]
            # add the compressed waveforms, if they exist
            if load_compressed and 'compressed_waveforms' in f:
                self.compressed_waveforms = {}
                for tmplt_hash in self.table['template_hash']:
                    self.compressed_waveforms[tmplt_hash] = \
                        pycbc.waveform.compress.CompressedWaveform.from_hdf(f,
                            tmplt_hash, load_now=load_compressed_now)
        else:
            raise ValueError("Unsupported template bank file extension %s" %(
                ext))

        # if approximant is specified, override whatever was in the file
        # (if anything was in the file)
        if approximant is not None:
            # get the approximant for each template
            apprxs = self.parse_approximant(approximant)
            if 'approximant' not in self.table.fieldnames:
                self.table = self.table.add_fields(apprxs, 'approximant')
            else:
                self.table['approximant'] = apprxs
        self.extra_args = kwds
        self.ensure_hash()

    @property
    def parameters(self):
        return self.table.fieldnames

    def ensure_hash(self):
        """Ensure that there is a correctly populated template_hash
        if it doesnt not already exist.
        """
        fields = self.table.fieldnames
        if 'template_hash' in fields:
             return

        # The fields to use in making a template hash
        hash_fields = ['mass1', 'mass2', 'inclination',
                       'spin1x', 'spin1y', 'spin1z',
                       'spin2x', 'spin2y', 'spin2z',]

        fields = [f for f in hash_fields if f in fields]
        template_hash = numpy.array([hash(v) for v in zip(*[self.table[p]
                                    for p in fields])])
        self.table = self.table.add_fields(template_hash, 'template_hash')

    def write_to_hdf(self, filename, start_index=None, stop_index=None,
                     force=False, skip_fields=None):
        """Writes self to the given hdf file.

        Parameters
        ----------
        filename : str
            The name of the file to write to. Must end in '.hdf'.
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

        Returns
        -------
        h5py.File
            The file handler to the output hdf file (left open).
        """
        if not filename.endswith('.hdf'):
            raise ValueError("Unrecoginized file extension")
        if os.path.exists(filename) and not force:
            raise IOError("File %s already exists" %(filename))
        f = h5py.File(filename, 'w')
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
        if self.compressed_waveforms is not None:
            for tmplt_hash in write_tbl.template_hash:
                self.compressed_waveforms[tmplt_hash].write_to_hdf(
                                                        f, tmplt_hash)
        return f

    def end_frequency(self, index):
        """ Return the end frequency of the waveform at the given index value
        """
        from pycbc.waveform.waveform import props

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
        return self.table["approximant"][index]

    def __len__(self):
        return len(self.table)

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
            cached_mem = zeros(wav_len, dtype=numpy.complex64)
        htilde = pycbc.waveform.get_waveform_filter(
            cached_mem, self.table[t_num], approximant=approximant,
            f_lower=low_frequency_cutoff, f_final=max_freq, delta_f=delta_f,
            distance=1./DYN_RANGE_FAC, delta_t=1./(2.*max_freq))
        return htilde

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
        thinning_bank = []
        tau0_temp, _ = pycbc.pnutils.mass1_mass2_to_tau0_tau3(m1, m2, fref)
        indices = []

        for inj in injection_parameters:
            tau0_inj, _ = \
                pycbc.pnutils.mass1_mass2_to_tau0_tau3(inj.mass1, inj.mass2,
                                                       fref)
            inj_indices = np.where(abs(tau0_temp - tau0_inj) <= threshold)[0]
            indices.append(inj_indices)
            indices_combined = np.concatenate(indices)

        indices_unique= np.unique(indices_combined)
        self.table = self.table[indices_unique]   


    def ensure_standard_filter_columns(self, low_frequency_cutoff=None):
        """ Initialize FilterBank common fields
        
        Parameters
        ---------
        low_frequency_cutoff: {float, None}, Optional
            A low frequency cutoff which overrides any given within the
            template bank file.
        """
        
        # Make sure we have a template duration field
        if not hasattr(self.table, 'template_duration'):
            self.table = self.table.add_fields(numpy.zeros(len(self.table),
                                     dtype=numpy.float32), 'template_duration') 

        # Make sure we have a f_lower field
        if low_frequency_cutoff is not None:
            if not hasattr(self.table, 'f_lower'):
                vec = numpy.zeros(len(self.table), dtype=numpy.float32)
                self.table = self.table.add_fields(vec, 'f_lower')
            self.table['f_lower'][:] = low_frequency_cutoff        

class LiveFilterBank(TemplateBank):
    def __init__(self, filename, sample_rate, minimum_buffer,
                       approximant=None, increment=8, parameters=None,
                       load_compressed=True, load_compressed_now=False, 
                       low_frequency_cutoff=None,
                       **kwds):

        self.increment = increment
        self.filename = filename
        self.sample_rate = sample_rate
        self.minimum_buffer = minimum_buffer

        super(LiveFilterBank, self).__init__(filename, approximant=approximant,
                parameters=parameters, load_compressed=load_compressed,
                load_compressed_now=load_compressed_now, **kwds)
        self.ensure_standard_filter_columns(low_frequency_cutoff=low_frequency_cutoff)
        self.table.sort(order='mchirp')
        self.hash_lookup = {}
        for i, p in enumerate(self.table):
            hash_value =  hash((p.mass1, p.mass2, p.spin1z, p.spin2z))
            self.hash_lookup[hash_value] = i

    def round_up(self, num):
        """Determine the length to use for this waveform by rounding.

        Parameters
        ----------
        num : int
            Proposed size of waveform in seconds

        Returns
        -------
        size: int
            The rounded size to use for the waveform buffer in seconds. This
        is calculaed using an internal `increment` attribute, which determines
        the discreteness of the rounding.
        """
        inc = self.increment
        size = numpy.ceil(num / self.sample_rate / inc) * self.sample_rate * inc
        return size

    def getslice(self, sindex):
        instance = copy(self)
        instance.table = self.table[sindex]
        return instance

    def id_from_hash(self, hash_value):
        """Get the index of this template based on its hash value

        Parameters
        ----------
        hash : int
            Value of the template hash

        Returns
        --------
        index : int
            The ordered index that this template has in the template bank.
        """
        return self.hash_lookup[hash_value]

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.getslice(index)

        approximant = self.approximant(index)
        f_end = self.end_frequency(index)
        flow = self.table[index].f_lower        

        # Determine the length of time of the filter, rounded up to
        # nearest power of two
        min_buffer = .5 + self.minimum_buffer

        from pycbc.waveform.waveform import props
        p = props(self.table[index])
        p.pop('approximant')
        buff_size = pycbc.waveform.get_waveform_filter_length_in_time(approximant, **p)
        
        tlen = self.round_up((buff_size + min_buffer) * self.sample_rate)
        flen = tlen / 2 + 1

        delta_f = self.sample_rate / float(tlen)

        if f_end is None or f_end >= (flen * delta_f):
            f_end = (flen-1) * delta_f

        logging.info("Generating %s, %ss, %i" % (approximant, 1.0/delta_f, index))

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        htilde = pycbc.waveform.get_waveform_filter(
            zeros(flen, dtype=numpy.complex64), self.table[index],
            approximant=approximant, f_lower=flow, f_final=f_end,
            delta_f=delta_f, delta_t=1.0/self.sample_rate, distance=distance,
            **self.extra_args)

        # If available, record the total duration (which may
        # include ringdown) and the duration up to merger since they will be
        # erased by the type conversion below.
        ttotal = template_duration = -1
        if hasattr(htilde, 'length_in_time'):
                ttotal = htilde.length_in_time
        if hasattr(htilde, 'chirp_length'):
                template_duration = htilde.chirp_length

        self.table[index].template_duration = template_duration

        htilde = htilde.astype(numpy.complex64)
        htilde.f_lower = flow
        htilde.end_idx = int(htilde.f_end / htilde.delta_f)
        htilde.params = self.table[index]
        htilde.chirp_length = template_duration
        htilde.length_in_time = ttotal
        htilde.approximant = approximant
        htilde.end_frequency = f_end

        # Add sigmasq as a method of this instance
        htilde.sigmasq = types.MethodType(sigma_cached, htilde)
        htilde._sigmasq = {}

        htilde.id = self.id_from_hash(hash((htilde.params.mass1,
                                      htilde.params.mass2,
                                      htilde.params.spin1z,
                                      htilde.params.spin2z)))
        return htilde

class FilterBank(TemplateBank):
    def __init__(self, filename, filter_length, delta_f, dtype,
                 out=None, max_template_length=None,
                 approximant=None, parameters=None,
                 load_compressed=True,
                 load_compressed_now=False,
                 low_frequency_cutoff=None,
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

        super(FilterBank, self).__init__(filename, approximant=approximant,
            parameters=parameters, load_compressed=load_compressed,
            load_compressed_now=load_compressed_now,
            **kwds)
        self.ensure_standard_filter_columns(low_frequency_cutoff=low_frequency_cutoff)

        self.min_f_lower = min(self.table['f_lower'])
        if self.f_lower is None and self.min_f_lower == 0.:
            raise ValueError('Invalid low-frequency cutoff settings')

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
        if self.f_lower is None:
            f_low = self.table[index].f_lower
        elif self.max_template_length is not None:
            f_low = find_variable_start_frequency(approximant,
                                                  self.table[index],
                                                  self.f_lower,
                                                  self.max_template_length)
        else:
            f_low = self.f_lower

        logging.info('%s: generating %s from %s Hz' % (index, approximant, f_low))

        # Clear the storage memory
        poke  = tempout.data
        tempout.clear()

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        htilde = pycbc.waveform.get_waveform_filter(
            tempout[0:self.filter_length], self.table[index],
            approximant=approximant, f_lower=f_low, f_final=f_end,
            delta_f=self.delta_f, delta_t=self.delta_t, distance=distance,
            **self.extra_args)

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
    l = max_length + 1
    f = f_start - delta_f
    while l > max_length:
        f += delta_f
        l = pycbc.waveform.get_waveform_filter_length_in_time(approximant,
                                                      parameters, f_lower=f)
    return f


class FilterBankSkyMax(TemplateBank):
    def __init__(self, filename, filter_length, delta_f, f_lower,
                 dtype, out_plus=None, out_cross=None,
                 max_template_length=None, parameters=None,
                 load_compressed=True, load_compressed_now=False,
                 **kwds):
        self.out_plus = out_plus
        self.out_cross = out_cross
        self.dtype = dtype
        self.f_lower = f_lower
        self.filename = filename
        self.delta_f = delta_f
        self.N = (filter_length - 1 ) * 2
        self.delta_t = 1.0 / (self.N * self.delta_f)
        self.filter_length = filter_length
        self.max_template_length = max_template_length

        super(FilterBankSkyMax, self).__init__(filename, parameters=parameters,
            load_compressed=True, load_compressed_now=False,
            **kwds)

        # add a template duration field if it already doesn't exist
        if not hasattr(self.table, 'template_duration'):
            if dtype == numpy.complex64 or dtype == numpy.float32:
                rdtype = numpy.float32
            else:
                rdtype = numpy.float64
            self.table = self.table.add_fields(numpy.zeros(len(self.table),
                                               dtype=rdtype),
                                               'template_duration')

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
        if self.max_template_length is not None:
            f_low = find_variable_start_frequency(approximant,
                                                  self.table[index],
                                                  self.f_lower,
                                                  self.max_template_length)
        else:
            f_low = self.f_lower

        logging.info('%s: generating %s from %s Hz' % (index, approximant, f_low))

        # What does this do???
        poke1 = tempoutplus.data
        poke2 = tempoutcross.data

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
