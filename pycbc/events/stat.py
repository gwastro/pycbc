# Copyright (C) 2016 Alex Nitz
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
This module contains functions for calculating coincident ranking statistic
values.
"""
import logging
import numpy
import h5py
from . import ranking
from . import coinc_rate
from .eventmgr_cython import logsignalrateinternals_computepsignalbins
from .eventmgr_cython import logsignalrateinternals_compute2detrate


class Stat(object):
    """Base class which should be extended to provide a coincident statistic"""

    def __init__(self, sngl_ranking, files=None, ifos=None, **kwargs):
        """
        Create a statistic class instance

        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.
        files: list of strs, needed for some statistics
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.
        ifos: list of strs, needed for some statistics
            The list of detector names
        """

        self.files = {}
        files = files or []
        for filename in files:
            with h5py.File(filename, 'r') as f:
                stat = f.attrs['stat']
            if hasattr(stat, 'decode'):
                stat = stat.decode()
            if stat in self.files:
                raise RuntimeError("We already have one file with stat attr ="
                                   " %s. Can't provide more than one!" % stat)
            logging.info("Found file %s for stat %s", filename, stat)
            self.files[stat] = filename

        # Provide the dtype of the single detector method's output
        # This is used by background estimation codes that need to maintain
        # a buffer of such values.
        self.single_dtype = numpy.float32
        # True if a larger single detector statistic will produce a larger
        # coincident statistic
        self.single_increasing = True

        self.ifos = ifos or []

        self.sngl_ranking = sngl_ranking
        self.sngl_ranking_kwargs = {}
        for key, value in kwargs.items():
            if key.startswith('sngl_ranking_'):
                self.sngl_ranking_kwargs[key[13:]] = value

    def get_sngl_ranking(self, trigs):
        """
        Returns the ranking for the single detector triggers.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return ranking.get_sngls_ranking_from_trigs(
            trigs,
            self.sngl_ranking,
            **self.sngl_ranking_kwargs
        )

    def single(self, trigs): # pylint:disable=unused-argument
        """
        Calculate the necessary single detector information

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        err_msg = "This function is a stub that should be overridden by the "
        err_msg += "sub-classes. You shouldn't be seeing this error!"
        raise NotImplementedError(err_msg)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for a single detector candidate

        Parameters
        ----------
        single_info: tuple
            Tuple containing two values. The first is the ifo (str) and the
            second is the single detector triggers.

        Returns
        -------
        numpy.ndarray
            The array of single detector statistics
        """
        err_msg = "This function is a stub that should be overridden by the "
        err_msg += "sub-classes. You shouldn't be seeing this error!"
        raise NotImplementedError(err_msg)

    def rank_stat_coinc(self, s, slide, step, to_shift,
                        **kwargs): # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.
        """
        err_msg = "This function is a stub that should be overridden by the "
        err_msg += "sub-classes. You shouldn't be seeing this error!"
        raise NotImplementedError(err_msg)

    def _check_coinc_lim_subclass(self, allowed_names):
        """
        Check that we are not using coinc_lim_for_thresh when not valid.

        coinc_lim_for_thresh is only defined for the statistic it is present
        in. If we subclass, we must check explicitly that it is still valid and
        indicate this in the code. If the code does not have this explicit
        check you will see the failure message here.

        Parameters
        -----------
        allowed_names : list
            list of allowed classes for the specific sub-classed method.
        """
        if type(self).__name__ not in allowed_names:
            err_msg = "This is being called from a subclass which has not "
            err_msg += "been checked for validity with this method. If it is "
            err_msg += "valid for the subclass to come here, include in the "
            err_msg += "list of allowed_names above."
            raise NotImplementedError(err_msg)

    def coinc_lim_for_thresh(self, s, thresh, limifo,
                             **kwargs): # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest

        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.
        """

        err_msg = "This function is a stub that should be overridden by the "
        err_msg += "sub-classes. You shouldn't be seeing this error!"
        raise NotImplementedError(err_msg)


class QuadratureSumStatistic(Stat):
    """Calculate the quadrature sum coincident detection statistic"""

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        Here just the ranking is computed and returned.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return self.get_sngl_ranking(trigs)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for a single detector candidate

        Parameters
        ----------
        single_info: tuple
            Tuple containing two values. The first is the ifo (str) and the
            second is the single detector triggers.

        Returns
        -------
        numpy.ndarray
            The array of single detector statistics
        """
        return self.single(single_info[1])

    def rank_stat_coinc(self, sngls_list, slide, step, to_shift,
                        **kwargs): # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.

        Parameters
        ----------
        sngls_list: list
            List of (ifo, single detector statistic) tuples
        slide: (unused in this statistic)
        step: (unused in this statistic)
        to_shift: list
            List of integers indicating what multiples of the time shift will
            be applied (unused in this statistic)

        Returns
        -------
        numpy.ndarray
            Array of coincident ranking statistic values
        """
        cstat = sum(sngl[1] ** 2. for sngl in sngls_list) ** 0.5
        # For single-detector "cuts" the single ranking is set to -1
        for sngls in sngls_list:
            cstat[sngls == -1] = 0
        return cstat

    def coinc_lim_for_thresh(self, s, thresh, limifo,
                             **kwargs): # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest

        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples for all detectors
            except limifo.
        thresh: float
            The threshold on the coincident statistic.
        limifo: string
            The ifo for which the limit is to be found.

        Returns
        -------
        numpy.ndarray
            Array of limits on the limifo single statistic to
            exceed thresh.
        """
        # Safety against subclassing and not rethinking this
        allowed_names = ['QuadratureSumStatistic']
        self._check_coinc_lim_subclass(allowed_names)

        s0 = thresh ** 2. - sum(sngl[1] ** 2. for sngl in s)
        s0[s0 < 0] = 0
        return s0 ** 0.5


class PhaseTDStatistic(QuadratureSumStatistic):
    """
    Statistic that re-weights combined newsnr using coinc parameters.

    The weighting is based on the PDF of time delays, phase differences and
    amplitude ratios between triggers in different ifos.
    """

    def __init__(self, sngl_ranking, files=None, ifos=None,
                 pregenerate_hist=True, **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.

        files: list of strs, unused here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.

        ifos: list of strs, needed here
            The list of detector names

        pregenerate_hist: bool, optional
            If False, do not pregenerate histogram on class instantiation.
            Default is True.
        """

        QuadratureSumStatistic.__init__(self, sngl_ranking, files=files,
                                        ifos=ifos, **kwargs)

        self.single_dtype = [
            ('snglstat', numpy.float32),
            ('coa_phase', numpy.float32),
            ('end_time', numpy.float64),
            ('sigmasq', numpy.float32),
            ('snr', numpy.float32)
        ]

        # Assign attribute so that it can be replaced with other functions
        self.has_hist = False
        self.hist_ifos = None
        self.ref_snr = 5.0
        self.relsense = {}
        self.swidth = self.pwidth = self.twidth = None
        self.srbmin = self.srbmax = None
        self.max_penalty = None
        self.pdtype = []
        self.weights = {}
        self.param_bin = {}
        self.two_det_flag = (len(ifos) == 2)
        self.two_det_weights = {}
        # Some memory
        self.pdif = numpy.zeros(128, dtype=numpy.float64)
        self.tdif = numpy.zeros(128, dtype=numpy.float64)
        self.sdif = numpy.zeros(128, dtype=numpy.float64)
        self.tbin = numpy.zeros(128, dtype=numpy.int32)
        self.pbin = numpy.zeros(128, dtype=numpy.int32)
        self.sbin = numpy.zeros(128, dtype=numpy.int32)

        if pregenerate_hist and not len(ifos) == 1:
            self.get_hist()

    def get_hist(self, ifos=None):
        """
        Read in a signal density file for the ifo combination

        Parameters
        ----------
        ifos: list
            The list of ifos. Needed if not given when initializing the class
            instance.
        """

        ifos = ifos or self.ifos

        selected = None
        for name in self.files:
            # Pick out the statistic files that provide phase / time/ amp
            # relationships and match to the ifos in use
            if 'phasetd_newsnr' in name:
                ifokey = name.split('_')[2]
                num = len(ifokey) / 2
                if num != len(ifos):
                    continue

                match = [ifo in ifokey for ifo in ifos]
                if False in match:
                    continue
                selected = name
                break

        if selected is None and len(ifos) > 1:
            raise RuntimeError("Couldn't figure out which stat file to use")

        logging.info("Using signal histogram %s for ifos %s", selected, ifos)
        weights = {}
        param = {}

        with h5py.File(self.files[selected], 'r') as histfile:
            self.hist_ifos = histfile.attrs['ifos']

            # Patch for pre-hdf5=3.0 histogram files
            try:
                logging.info("Decoding hist ifos ..")
                self.hist_ifos = [i.decode('UTF-8') for i in self.hist_ifos]
            except (UnicodeDecodeError, AttributeError):
                pass

            # Histogram bin attributes
            self.twidth = histfile.attrs['twidth']
            self.pwidth = histfile.attrs['pwidth']
            self.swidth = histfile.attrs['swidth']
            self.srbmin = histfile.attrs['srbmin']
            self.srbmax = histfile.attrs['srbmax']
            relfac = histfile.attrs['sensitivity_ratios']

            for ifo in self.hist_ifos:
                weights[ifo] = histfile[ifo]['weights'][:]
                param[ifo] = histfile[ifo]['param_bin'][:]

        n_ifos = len(self.hist_ifos)

        bin_volume = (self.twidth * self.pwidth * self.swidth) ** (n_ifos - 1)
        self.hist_max = - 1. * numpy.inf

        # Read histogram for each ifo, to use if that ifo has smallest SNR in
        # the coinc
        for ifo in self.hist_ifos:

            # renormalise to PDF
            self.weights[ifo] = \
                weights[ifo] / (weights[ifo].sum() * bin_volume)

            if param[ifo].dtype == numpy.int8:
                # Older style, incorrectly sorted histogram file
                ncol = param[ifo].shape[1]
                self.pdtype = [('c%s' % i, param[ifo].dtype) for i in range(ncol)]
                self.param_bin[ifo] = numpy.zeros(len(self.weights[ifo]),
                                                  dtype=self.pdtype)
                for i in range(ncol):
                    self.param_bin[ifo]['c%s' % i] = param[ifo][:, i]

                lsort = self.param_bin[ifo].argsort()
                self.param_bin[ifo] = self.param_bin[ifo][lsort]
                self.weights[ifo] = self.weights[ifo][lsort]
            else:
                # New style, efficient histogram file
                # param bin and weights have already been sorted
                self.param_bin[ifo] = param[ifo]
                self.pdtype = self.param_bin[ifo].dtype

            # Max_penalty is a small number to assigned to any bins without
            # histogram entries. All histograms in a given file have the same
            # min entry by design, so use the min of the last one read in.
            self.max_penalty = self.weights[ifo].min()
            self.hist_max = max(self.hist_max, self.weights[ifo].max())

            if self.two_det_flag:
                # The density of signals is computed as a function of 3 binned
                # parameters: time difference (t), phase difference (p) and
                # SNR ratio (s). These are computed for each combination of
                # detectors, so for detectors 6 differences are needed. However
                # many combinations of these parameters are highly unlikely and
                # no instances of these combinations occurred when generating
                # the statistic files. Rather than storing a bunch of 0s, these
                # values are just not stored at all. This reduces the size of
                # the statistic file, but means we have to identify the correct
                # value to read for every trigger. For 2 detectors we can
                # expand the weights lookup table here, basically adding in all
                # the "0" values. This makes looking up a value in the
                # "weights" table a O(N) rather than O(NlogN) operation. It
                # sacrifices RAM to do this, so is a good tradeoff for 2
                # detectors, but not for 3!
                if not hasattr(self, 'c0_size'):
                    self.c0_size = {}
                    self.c1_size = {}
                    self.c2_size = {}

                self.c0_size[ifo] = numpy.int32(
                    2 * (abs(self.param_bin[ifo]['c0']).max() + 1)
                )
                self.c1_size[ifo] = numpy.int32(
                    2 * (abs(self.param_bin[ifo]['c1']).max() + 1)
                )
                self.c2_size[ifo] = numpy.int32(
                    2 * (abs(self.param_bin[ifo]['c2']).max() + 1)
                )

                array_size = [self.c0_size[ifo], self.c1_size[ifo],
                              self.c2_size[ifo]]
                dtypec = self.weights[ifo].dtype
                self.two_det_weights[ifo] = \
                    numpy.zeros(array_size, dtype=dtypec) + self.max_penalty
                id0 = self.param_bin[ifo]['c0'].astype(numpy.int32) \
                    + self.c0_size[ifo] // 2
                id1 = self.param_bin[ifo]['c1'].astype(numpy.int32) \
                    + self.c1_size[ifo] // 2
                id2 = self.param_bin[ifo]['c2'].astype(numpy.int32) \
                    + self.c2_size[ifo] // 2
                self.two_det_weights[ifo][id0, id1, id2] = self.weights[ifo]

        for ifo, sense in zip(self.hist_ifos, relfac):
            self.relsense[ifo] = sense

        self.has_hist = True

    def logsignalrate(self, stats, shift, to_shift):
        """
        Calculate the normalized log rate density of signals via lookup

        Parameters
        ----------
        stats: dict of dicts
            Single-detector quantities for each detector
        shift: numpy array of float
            Time shift vector for each coinc to be ranked
        to_shift: list of ints
            Multiple of the time shift to apply, ordered as self.ifos

        Returns
        -------
        value: log of coinc signal rate density for the given single-ifo
            triggers and time shifts
        """
        # Convert time shift vector to dict, as hist ifos and self.ifos may
        # not be in same order
        to_shift = {ifo: s for ifo, s in zip(self.ifos, to_shift)}

        if not self.has_hist:
            self.get_hist()

        # Figure out which ifo of the contributing ifos has the smallest SNR,
        # to use as reference for choosing the signal histogram.
        snrs = numpy.array([numpy.array(stats[ifo]['snr'], ndmin=1)
                           for ifo in self.ifos])
        smin = snrs.argmin(axis=0)
        # Store a list of the triggers using each ifo as reference
        rtypes = {ifo: numpy.where(smin == j)[0]
                  for j, ifo in enumerate(self.ifos)}

        # Get reference ifo information
        rate = numpy.zeros(len(shift), dtype=numpy.float32)
        ps = {ifo: numpy.array(stats[ifo]['coa_phase'], ndmin=1)
              for ifo in self.ifos}
        ts = {ifo: numpy.array(stats[ifo]['end_time'], ndmin=1)
              for ifo in self.ifos}
        ss = {ifo: numpy.array(stats[ifo]['snr'], ndmin=1)
              for ifo in self.ifos}
        sigs = {ifo: numpy.array(stats[ifo]['sigmasq'], ndmin=1)
                for ifo in self.ifos}
        for ref_ifo in self.ifos:
            rtype = rtypes[ref_ifo]
            pref = ps[ref_ifo]
            tref = ts[ref_ifo]
            sref = ss[ref_ifo]
            sigref = sigs[ref_ifo]
            senseref = self.relsense[self.hist_ifos[0]]

            binned = []
            other_ifos = [ifo for ifo in self.ifos if ifo != ref_ifo]
            for ifo in other_ifos:
                # Assign cached memory
                length = len(rtype)
                while length > len(self.pdif):
                    newlen = len(self.pdif) * 2
                    self.pdif = numpy.zeros(newlen, dtype=numpy.float64)
                    self.tdif = numpy.zeros(newlen, dtype=numpy.float64)
                    self.sdif = numpy.zeros(newlen, dtype=numpy.float64)
                    self.pbin = numpy.zeros(newlen, dtype=numpy.int32)
                    self.tbin = numpy.zeros(newlen, dtype=numpy.int32)
                    self.sbin = numpy.zeros(newlen, dtype=numpy.int32)

                # Calculate differences
                logsignalrateinternals_computepsignalbins(
                    self.pdif,
                    self.tdif,
                    self.sdif,
                    self.pbin,
                    self.tbin,
                    self.sbin,
                    ps[ifo],
                    ts[ifo],
                    ss[ifo],
                    sigs[ifo],
                    pref,
                    tref,
                    sref,
                    sigref,
                    shift,
                    rtype,
                    self.relsense[ifo],
                    senseref,
                    self.twidth,
                    self.pwidth,
                    self.swidth,
                    to_shift[ref_ifo],
                    to_shift[ifo],
                    length
                )

                binned += [
                    self.tbin[:length],
                    self.pbin[:length],
                    self.sbin[:length]
                ]

            # Read signal weight from precalculated histogram
            if self.two_det_flag:
                # High-RAM, low-CPU option for two-det
                logsignalrateinternals_compute2detrate(
                    binned[0],
                    binned[1],
                    binned[2],
                    self.c0_size[ref_ifo],
                    self.c1_size[ref_ifo],
                    self.c2_size[ref_ifo],
                    rate,
                    rtype,
                    sref,
                    self.two_det_weights[ref_ifo],
                    self.max_penalty,
                    self.ref_snr,
                    len(rtype)
                )
            else:
                # Low[er]-RAM, high[er]-CPU option for >two det

                # Convert binned to same dtype as stored in hist
                nbinned = numpy.zeros(len(binned[1]), dtype=self.pdtype)
                for i, b in enumerate(binned):
                    nbinned[f'c{i}'] = b

                loc = numpy.searchsorted(self.param_bin[ref_ifo], nbinned)
                loc[loc == len(self.weights[ref_ifo])] = 0
                rate[rtype] = self.weights[ref_ifo][loc]

                # These weren't in our histogram so give them max penalty
                # instead of random value
                missed = numpy.where(
                    self.param_bin[ref_ifo][loc] != nbinned
                )[0]
                rate[rtype[missed]] = self.max_penalty
                # Scale by signal population SNR
                rate[rtype] *= (sref[rtype] / self.ref_snr) ** -4.0

        return numpy.log(rate)

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        Here the ranking as well as phase, endtime and sigma-squared values.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information. 'snr', 'chisq',
            'chisq_dof', 'coa_phase', 'end_time', and 'sigmasq' are required
            keys.

        Returns
        -------
        numpy.ndarray
            Array of single detector parameter values
        """
        sngl_stat = self.get_sngl_ranking(trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase'][:]
        singles['end_time'] = trigs['end_time'][:]
        singles['sigmasq'] = trigs['sigmasq'][:]
        singles['snr'] = trigs['snr'][:]
        return numpy.array(singles, ndmin=1)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for a single detector candidate

        Parameters
        ----------
        single_info: tuple
            Tuple containing two values. The first is the ifo (str) and the
            second is the single detector triggers.

        Returns
        -------
        numpy.ndarray
            The array of single detector statistics
        """
        return self.single(single_info[1])

    def rank_stat_coinc(self, sngls_list, slide, step, to_shift,
                        **kwargs):  # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic, defined in Eq 2 of
        [Nitz et al, 2017](https://doi.org/10.3847/1538-4357/aa8f50).
        """
        rstat = sum(s[1]['snglstat'] ** 2 for s in sngls_list)
        cstat = rstat + 2. * self.logsignalrate(dict(sngls_list),
                                                slide * step,
                                                to_shift)
        cstat[cstat < 0] = 0
        return cstat ** 0.5

    def coinc_lim_for_thresh(self, sngls_list, thresh, limifo,
                             **kwargs):  # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest.
        Calculate the required single detector statistic to exceed the
        threshold for each of the input triggers.
        """
        # Safety against subclassing and not rethinking this
        allowed_names = ['PhaseTDStatistic']
        self._check_coinc_lim_subclass(allowed_names)

        if not self.has_hist:
            self.get_hist()

        lim_stat = [b['snglstat'] for a, b in sngls_list if a == limifo][0]
        s1 = thresh ** 2. - lim_stat ** 2.
        # Assume best case scenario and use maximum signal rate
        s1 -= 2. * self.hist_max
        s1[s1 < 0] = 0
        return s1 ** 0.5


class ExpFitStatistic(QuadratureSumStatistic):
    """
    Detection statistic using an exponential falloff noise model.

    Statistic approximates the negative log noise coinc rate density per
    template over single-ifo newsnr values.
    """

    def __init__(self, sngl_ranking, files=None, ifos=None, **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.

        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.

        ifos: list of strs, not used here
            The list of detector names
        """

        if not files:
            raise RuntimeError("Statistic files not specified")
        QuadratureSumStatistic.__init__(self, sngl_ranking, files=files,
                                        ifos=ifos, **kwargs)

        # the stat file attributes are hard-coded as '%{ifo}-fit_coeffs'
        parsed_attrs = [f.split('-') for f in self.files.keys()]
        self.bg_ifos = [at[0] for at in parsed_attrs if
                       (len(at) == 2 and at[1] == 'fit_coeffs')]
        if not len(self.bg_ifos):
            raise RuntimeError("None of the statistic files has the required "
                               "attribute called {ifo}-fit_coeffs !")
        self.fits_by_tid = {}
        self.alphamax = {}
        for i in self.bg_ifos:
            self.fits_by_tid[i] = self.assign_fits(i)
            self.get_ref_vals(i)

        self.single_increasing = False

    def assign_fits(self, ifo):
        """
        Extract fits from fit files

        Parameters
        -----------
        ifo: str
            The detector to get fits for.

        Returns
        -------
        rate_dict: dict
            A dictionary containing the fit information in the `alpha`, `rate`
            and `thresh` keys/.
        """
        coeff_file = h5py.File(self.files[f'{ifo}-fit_coeffs'], 'r')
        template_id = coeff_file['template_id'][:]
        # the template_ids and fit coeffs are stored in an arbitrary order
        # create new arrays in template_id order for easier recall
        tid_sort = numpy.argsort(template_id)

        fits_by_tid_dict = {}
        fits_by_tid_dict['smoothed_fit_coeff'] = \
            coeff_file['fit_coeff'][:][tid_sort]
        fits_by_tid_dict['smoothed_rate_above_thresh'] = \
            coeff_file['count_above_thresh'][:][tid_sort].astype(float)
        fits_by_tid_dict['smoothed_rate_in_template'] = \
            coeff_file['count_in_template'][:][tid_sort].astype(float)

        # The by-template fits may have been stored in the smoothed fits file
        if 'fit_by_template' in coeff_file:
            coeff_fbt = coeff_file['fit_by_template']
            fits_by_tid_dict['fit_by_fit_coeff'] = \
                coeff_fbt['fit_coeff'][:][tid_sort]
            fits_by_tid_dict['fit_by_rate_above_thresh'] = \
                coeff_fbt['count_above_thresh'][:][tid_sort].astype(float)
            fits_by_tid_dict['fit_by_rate_in_template'] = \
                coeff_file['count_in_template'][:][tid_sort].astype(float)

        # Keep the fit threshold in fits_by_tid
        fits_by_tid_dict['thresh'] = coeff_file.attrs['stat_threshold']

        coeff_file.close()

        return fits_by_tid_dict

    def get_ref_vals(self, ifo):
        """
        Get the largest `alpha` value over all templates for given ifo.

        This is stored in `self.alphamax[ifo]` in the class instance.

        Parameters
        -----------
        ifo: str
            The detector to get fits for.
        """
        self.alphamax[ifo] = self.fits_by_tid[ifo]['smoothed_fit_coeff'].max()

    def find_fits(self, trigs):
        """
        Get fit coeffs for a specific ifo and template id(s)

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.
            The coincidence executable will always call this using a bunch of
            trigs from a single template, there template_num is stored as an
            attribute and we just return the single value for all templates.
            If multiple templates are in play we must return arrays.

        Returns
        --------
        alphai: float or numpy array
            The alpha fit value(s)
        ratei: float or numpy array
            The rate fit value(s)
        thresh: float or numpy array
            The thresh fit value(s)
        """

        try:
            tnum = trigs.template_num  # exists if accessed via coinc_findtrigs
            ifo = trigs.ifo
        except AttributeError:
            tnum = trigs['template_id']  # exists for SingleDetTriggers
            assert len(self.ifos) == 1
            # Should be exactly one ifo provided
            ifo = self.ifos[0]
        # fits_by_tid is a dictionary of dictionaries of arrays
        # indexed by ifo / coefficient name / template_id
        alphai = self.fits_by_tid[ifo]['smoothed_fit_coeff'][tnum]
        ratei = self.fits_by_tid[ifo]['smoothed_rate_above_thresh'][tnum]
        thresh = self.fits_by_tid[ifo]['thresh']

        return alphai, ratei, thresh

    def lognoiserate(self, trigs):
        """
        Calculate the log noise rate density over single-ifo ranking

        Read in single trigger information, compute the ranking
        and rescale by the fitted coefficients alpha and rate

        Parameters
        -----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        ---------
        lognoisel: numpy.array
            Array of log noise rate density for each input trigger.

        """
        alphai, ratei, thresh = self.find_fits(trigs)
        sngl_stat = self.get_sngl_ranking(trigs)
        # alphai is constant of proportionality between single-ifo newsnr and
        #   negative log noise likelihood in given template
        # ratei is rate of trigs in given template compared to average
        # thresh is stat threshold used in given ifo
        lognoisel = - alphai * (sngl_stat - thresh) + numpy.log(alphai) + \
                      numpy.log(ratei)
        return numpy.array(lognoisel, ndmin=1, dtype=numpy.float32)

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        In this case the ranking rescaled (see the lognoiserate method here).

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """

        return self.lognoiserate(trigs)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for a single detector candidate

        Parameters
        ----------
        single_info: tuple
            Tuple containing two values. The first is the ifo (str) and the
            second is the single detector triggers.

        Returns
        -------
        numpy.ndarray
            The array of single detector statistics
        """
        err_msg = "Sorry! No-one has implemented this method yet! "
        raise NotImplementedError(err_msg)

    def rank_stat_coinc(self, s, slide, step, to_shift,
                        **kwargs): # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.
        """
        err_msg = "Sorry! No-one has implemented this method yet! "
        raise NotImplementedError(err_msg)

    def coinc_lim_for_thresh(self, s, thresh, limifo,
                             **kwargs): # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest
        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.
        """
        err_msg = "Sorry! No-one has implemented this method yet! "
        raise NotImplementedError(err_msg)

    # Keeping this here to help write the new coinc method.
    def coinc_OLD(self, s0, s1, slide, step): # pylint:disable=unused-argument
        """Calculate the final coinc ranking statistic"""

        # Approximate log likelihood ratio by summing single-ifo negative
        # log noise likelihoods
        loglr = - s0 - s1
        # add squares of threshold stat values via idealized Gaussian formula
        threshes = [self.fits_by_tid[i]['thresh'] for i in self.bg_ifos]
        loglr += sum([t**2. / 2. for t in threshes])
        # convert back to a coinc-SNR-like statistic
        # via log likelihood ratio \propto rho_c^2 / 2
        return (2. * loglr) ** 0.5

    # Keeping this here to help write the new coinc_lim method
    def coinc_lim_for_thresh_OLD(self, s0, thresh):
        """Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.

        Parameters
        ----------
        s0: numpy.ndarray
            Single detector ranking statistic for the first detector.
        thresh: float
            The threshold on the coincident statistic.

        Returns
        -------
        numpy.ndarray
            Array of limits on the second detector single statistic to
        exceed thresh.
        """
        s1 = - (thresh ** 2.) / 2. - s0
        threshes = [self.fits_by_tid[i]['thresh'] for i in self.bg_ifos]
        s1 += sum([t**2. / 2. for t in threshes])
        return s1


class ExpFitCombinedSNR(ExpFitStatistic):
    """
    Reworking of ExpFitStatistic designed to resemble network SNR

    Use a monotonic function of the negative log noise rate density which
    approximates combined (new)snr for coincs with similar newsnr in each ifo
    """

    def __init__(self, sngl_ranking, files=None, ifos=None, **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.

        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.

        ifos: list of strs, not used here
            The list of detector names
        """

        ExpFitStatistic.__init__(self, sngl_ranking, files=files, ifos=ifos,
                                 **kwargs)
        # for low-mass templates the exponential slope alpha \approx 6
        self.alpharef = 6.
        self.single_increasing = True

    def use_alphamax(self):
        """
        Compute the reference alpha from the fit files.

        Use the harmonic mean of the maximum individual ifo slopes as the
        reference value of alpha.
        """

        inv_alphas = [1. / self.alphamax[i] for i in self.bg_ifos]
        self.alpharef = 1. / (sum(inv_alphas) / len(inv_alphas))

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """

        logr_n = self.lognoiserate(trigs)
        _, _, thresh = self.find_fits(trigs)
        # shift by log of reference slope alpha
        logr_n += -1. * numpy.log(self.alpharef)
        # add threshold and rescale by reference slope
        stat = thresh - (logr_n / self.alpharef)
        return numpy.array(stat, ndmin=1, dtype=numpy.float32)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for single detector candidates

        Parameters
        ----------
        single_info: tuple
            Tuple containing two values. The first is the ifo (str) and the
            second is the single detector triggers.

        Returns
        -------
        numpy.ndarray
            The array of single detector statistics
        """
        if self.single_increasing:
            sngl_multiifo = single_info[1]['snglstat']
        else:
            sngl_multiifo = -1.0 * single_info[1]['snglstat']
        return sngl_multiifo

    def rank_stat_coinc(self, s, slide, step, to_shift,
                        **kwargs): # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.

        Parameters
        ----------
        sngls_list: list
            List of (ifo, single detector statistic) tuples
        slide: (unused in this statistic)
        step: (unused in this statistic)
        to_shift: list
            List of integers indicating what multiples of the time shift will
            be applied (unused in this statistic)

        Returns
        -------
        numpy.ndarray
            Array of coincident ranking statistic values
        """

        # scale by 1/sqrt(number of ifos) to resemble network SNR
        return sum(sngl[1] for sngl in s) / len(s)**0.5

    def coinc_lim_for_thresh(self, s, thresh, limifo,
                             **kwargs): # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest

        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples for all detectors
            except limifo.
        thresh: float
            The threshold on the coincident statistic.
        limifo: string
            The ifo for which the limit is to be found.

        Returns
        -------
        numpy.ndarray
            Array of limits on the limifo single statistic to
            exceed thresh.
        """
        # Safety against subclassing and not rethinking this
        allowed_names = ['ExpFitCombinedSNR']
        self._check_coinc_lim_subclass(allowed_names)

        return thresh * ((len(s) + 1) ** 0.5) - sum(sngl[1] for sngl in s)


class PhaseTDExpFitStatistic(PhaseTDStatistic, ExpFitCombinedSNR):
    """
    Statistic combining exponential noise model with signal histogram PDF
    """

    def __init__(self, sngl_ranking, files=None, ifos=None, **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.
        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.
        ifos: list of strs, needed here
            The list of detector names
        """

        # read in both foreground PDF and background fit info
        ExpFitCombinedSNR.__init__(self, sngl_ranking, files=files, ifos=ifos,
                                   **kwargs)
        # need the self.single_dtype value from PhaseTDStatistic
        PhaseTDStatistic.__init__(self, sngl_ranking, files=files,
                                  ifos=ifos, **kwargs)

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        In this case the ranking rescaled (see the lognoiserate method here)
        with the phase, end time, sigma and SNR values added in.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """

        # same single-ifo stat as ExpFitCombinedSNR
        sngl_stat = ExpFitCombinedSNR.single(self, trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase'][:]
        singles['end_time'] = trigs['end_time'][:]
        singles['sigmasq'] = trigs['sigmasq'][:]
        singles['snr'] = trigs['snr'][:]
        return numpy.array(singles, ndmin=1)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for a single detector candidate

        Parameters
        ----------
        single_info: tuple
            Tuple containing two values. The first is the ifo (str) and the
            second is the single detector triggers.

        Returns
        -------
        numpy.ndarray
            The array of single detector statistics
        """
        err_msg = "Sorry! No-one has implemented this method yet! "
        raise NotImplementedError(err_msg)

    def rank_stat_coinc(self, s, slide, step, to_shift,
                        **kwargs): # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.
        """
        err_msg = "Sorry! No-one has implemented this method yet! "
        raise NotImplementedError(err_msg)

    def coinc_lim_for_thresh(self, s, thresh, limifo,
                             **kwargs): # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest
        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.
        """
        err_msg = "Sorry! No-one has implemented this method yet! "
        raise NotImplementedError(err_msg)

    # Keeping the old statistic code here for now to help with reimplementing
    def coinc_OLD(self, s0, s1, slide, step):
        # logsignalrate function inherited from PhaseTDStatistic
        logr_s = self.logsignalrate(s0, s1, slide * step)
        # rescale by ExpFitCombinedSNR reference slope as for sngl stat
        cstat = s0['snglstat'] + s1['snglstat'] + logr_s / self.alpharef
        # cut off underflowing and very small values
        cstat[cstat < 8.] = 8.
        # scale to resemble network SNR
        return cstat / (2.**0.5)

    def coinc_lim_for_thresh_OLD(self, s0, thresh):
        # if the threshold is below this value all triggers will
        # pass because of rounding in the coinc method
        if thresh <= (8. / (2.**0.5)):
            return -1. * numpy.ones(len(s0['snglstat'])) * numpy.inf
        if not self.has_hist:
            self.get_hist()
        # Assume best case scenario and use maximum signal rate
        logr_s = self.hist_max
        s1 = (2 ** 0.5) * thresh - s0['snglstat'] - logr_s / self.alpharef
        return s1


class ExpFitBgRateStatistic(ExpFitStatistic):
    """
    Detection statistic using an exponential falloff noise model.

    Statistic calculates the log noise coinc rate for each
    template over single-ifo newsnr values.
    """

    def __init__(self, sngl_ranking, files=None, ifos=None,
                 benchmark_lograte=-14.6, **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.
        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.
        ifos: list of strs, not used here
            The list of detector names
        benchmark_lograte: float, default=-14.6
            benchmark_lograte is log of a representative noise trigger rate.
            The default comes from H1L1 (O2) and is 4.5e-7 Hz.
        """

        super(ExpFitBgRateStatistic, self).__init__(sngl_ranking,
                                                    files=files, ifos=ifos,
                                                    **kwargs)
        self.benchmark_lograte = benchmark_lograte

        # Reassign the rate to be number per time rather than an arbitrarily
        # normalised number
        for ifo in self.bg_ifos:
            self.reassign_rate(ifo)

    def reassign_rate(self, ifo):
        """
        Reassign the rate to be number per time rather

        Reassign the rate to be number per time rather than an arbitrarily
        normalised number.

        Parameters
        -----------
        ifo: str
            The ifo to consider.
        """
        with h5py.File(self.files[f'{ifo}-fit_coeffs'], 'r') as coeff_file:
            analysis_time = float(coeff_file.attrs['analysis_time'])
            fbt = 'fit_by_template' in coeff_file

        self.fits_by_tid[ifo]['smoothed_rate_above_thresh'] /= analysis_time
        self.fits_by_tid[ifo]['smoothed_rate_in_template'] /= analysis_time
        # The by-template fits may have been stored in the smoothed fits file
        if fbt:
            self.fits_by_tid[ifo]['fit_by_rate_above_thresh'] /= analysis_time
            self.fits_by_tid[ifo]['fit_by_rate_in_template'] /= analysis_time

    def rank_stat_coinc(self, s, slide, step, to_shift,
                        **kwargs): # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.

        Parameters
        ----------
        sngls_list: list
            List of (ifo, single detector statistic) tuples
        slide: (unused in this statistic)
        step: (unused in this statistic)
        to_shift: list
            List of integers indicating what multiples of the time shift will
            be applied (unused in this statistic)

        Returns
        -------
        numpy.ndarray
            Array of coincident ranking statistic values
        """

        # ranking statistic is -ln(expected rate density of noise triggers)
        # plus normalization constant
        sngl_dict = {sngl[0]: sngl[1] for sngl in s}
        ln_noise_rate = coinc_rate.combination_noise_lograte(
                                  sngl_dict, kwargs['time_addition'])
        loglr = - ln_noise_rate + self.benchmark_lograte
        return loglr

    def coinc_lim_for_thresh(self, s, thresh, limifo, **kwargs):
        """
        Optimization function to identify coincs too quiet to be of interest

        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples for all detectors
            except limifo.
        thresh: float
            The threshold on the coincident statistic.
        limifo: string
            The ifo for which the limit is to be found.

        Returns
        -------
        numpy.ndarray
            Array of limits on the limifo single statistic to
            exceed thresh.
        """

        # Safety against subclassing and not rethinking this
        allowed_names = ['ExpFitBgRateStatistic']
        self._check_coinc_lim_subclass(allowed_names)

        sngl_dict = {sngl[0]: sngl[1] for sngl in s}
        sngl_dict[limifo] = numpy.zeros(len(s[0][1]))
        ln_noise_rate = coinc_rate.combination_noise_lograte(
                                  sngl_dict, kwargs['time_addition'])
        loglr = - thresh - ln_noise_rate + self.benchmark_lograte
        return loglr


class ExpFitFgBgNormStatistic(PhaseTDStatistic,
                              ExpFitBgRateStatistic):
    """
    Statistic combining PhaseTD, ExpFitBg and additional foreground info.
    """

    def __init__(self, sngl_ranking, files=None, ifos=None,
                 reference_ifos='H1,L1', **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.
        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.
        ifos: list of strs
            The list of detector names
        reference_ifos: string of comma separated ifo prefixes
            Detectors to be used as the reference network for network
            sensitivity comparisons. Each must be in fits_by_tid
        """

        # read in background fit info and store it
        ExpFitBgRateStatistic.__init__(self, sngl_ranking, files=files,
                                       ifos=ifos, **kwargs)
        # if ifos not already set, determine via background fit info
        self.ifos = self.ifos or self.bg_ifos
        # PhaseTD statistic single_dtype plus network sensitivity benchmark
        PhaseTDStatistic.__init__(self, sngl_ranking, files=files,
                                  ifos=self.ifos, **kwargs)
        self.single_dtype.append(('benchmark_logvol', numpy.float32))

        for ifo in self.bg_ifos:
            self.assign_median_sigma(ifo)

        ref_ifos = reference_ifos.split(',')

        # benchmark_logvol is a benchmark sensitivity array over template id
        hl_net_med_sigma = numpy.amin([self.fits_by_tid[ifo]['median_sigma']
                                       for ifo in ref_ifos], axis=0)
        self.benchmark_logvol = 3.0 * numpy.log(hl_net_med_sigma)
        self.single_increasing = False

    def assign_median_sigma(self, ifo):
        """
        Read and sort the median_sigma values from input files.

        Parameters
        ----------
        ifo: str
            The ifo to consider.
        """

        with h5py.File(self.files[f'{ifo}-fit_coeffs'], 'r') as coeff_file:
            template_id = coeff_file['template_id'][:]
            tid_sort = numpy.argsort(template_id)
            self.fits_by_tid[ifo]['median_sigma'] = \
                coeff_file['median_sigma'][:][tid_sort]

    def lognoiserate(self, trigs, alphabelow=6):
        """
        Calculate the log noise rate density over single-ifo ranking

        Read in single trigger information, make the newsnr statistic
        and rescale by the fitted coefficients alpha and rate

        Parameters
        -----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.
        alphabelow: float, default=6
            Use this slope to fit the noise triggers below the point at which
            fits are present in the input files.

        Returns
        ---------
        lognoisel: numpy.array
            Array of log noise rate density for each input trigger.
        """
        alphai, ratei, thresh = self.find_fits(trigs)
        newsnr = self.get_sngl_ranking(trigs)
        # Above the threshold we use the usual fit coefficient (alpha)
        # below threshold use specified alphabelow
        bt = newsnr < thresh
        lognoisel = - alphai * (newsnr - thresh) + numpy.log(alphai) + \
                        numpy.log(ratei)
        lognoiselbt = - alphabelow * (newsnr - thresh) + \
                           numpy.log(alphabelow) + numpy.log(ratei)
        lognoisel[bt] = lognoiselbt[bt]
        return numpy.array(lognoisel, ndmin=1, dtype=numpy.float32)

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        In this case the ranking rescaled (see the lognoiserate method here)
        with the phase, end time, sigma, SNR, template_id and the
        benchmark_logvol values added in.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.
        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """

        # single-ifo stat = log of noise rate
        sngl_stat = self.lognoiserate(trigs)
        # populate other fields to calculate phase/time/amp consistency
        # and sigma comparison
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase'][:]
        singles['end_time'] = trigs['end_time'][:]
        singles['sigmasq'] = trigs['sigmasq'][:]
        singles['snr'] = trigs['snr'][:]
        try:
            tnum = trigs.template_num  # exists if accessed via coinc_findtrigs
        except AttributeError:
            tnum = trigs['template_id']  # exists for SingleDetTriggers
            # Should only be one ifo fit file provided
            assert len(self.ifos) == 1
        # Store benchmark log volume as single-ifo information since the coinc
        # method does not have access to template id
        singles['benchmark_logvol'] = self.benchmark_logvol[tnum]
        return numpy.array(singles, ndmin=1)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for single detector candidates

        Parameters
        ----------
        single_info: tuple
            Tuple containing two values. The first is the ifo (str) and the
            second is the single detector triggers.

        Returns
        -------
        numpy.ndarray
            The array of single detector statistics
        """
        sngls = single_info[1]

        ln_noise_rate = sngls['snglstat']
        ln_noise_rate -= self.benchmark_lograte
        network_sigmasq = sngls['sigmasq']
        network_logvol = 1.5 * numpy.log(network_sigmasq)
        benchmark_logvol = sngls['benchmark_logvol']
        network_logvol -= benchmark_logvol
        ln_s = -4 * numpy.log(sngls['snr'] / self.ref_snr)
        loglr = network_logvol - ln_noise_rate + ln_s
        # cut off underflowing and very small values
        loglr[loglr < -30.] = -30.
        return loglr

    def rank_stat_coinc(self, s, slide, step, to_shift,
                        **kwargs): # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.

        Parameters
        ----------
        sngls_list: list
            List of (ifo, single detector statistic) tuples
        slide: (unused in this statistic)
        step: (unused in this statistic)
        to_shift: list
            List of integers indicating what multiples of the time shift will
            be applied (unused in this statistic)

        Returns
        -------
        numpy.ndarray
            Array of coincident ranking statistic values
        """

        sngl_rates = {sngl[0]: sngl[1]['snglstat'] for sngl in s}
        ln_noise_rate = coinc_rate.combination_noise_lograte(
                                  sngl_rates, kwargs['time_addition'])
        ln_noise_rate -= self.benchmark_lograte

        # Network sensitivity for a given coinc type is approximately
        # determined by the least sensitive ifo
        network_sigmasq = numpy.amin([sngl[1]['sigmasq'] for sngl in s],
                                     axis=0)
        # Volume \propto sigma^3 or sigmasq^1.5
        network_logvol = 1.5 * numpy.log(network_sigmasq)
        # Get benchmark log volume as single-ifo information :
        # benchmark_logvol for a given template is not ifo-dependent, so
        # choose the first ifo for convenience
        benchmark_logvol = s[0][1]['benchmark_logvol']
        network_logvol -= benchmark_logvol

        # Use prior histogram to get Bayes factor for signal vs noise
        # given the time, phase and SNR differences between IFOs

        # First get signal PDF logr_s
        stat = {ifo: st for ifo, st in s}
        logr_s = self.logsignalrate(stat, slide * step, to_shift)

        # Find total volume of phase-time-amplitude space occupied by noise
        # coincs
        # Extent of time-difference space occupied
        noise_twindow = coinc_rate.multiifo_noise_coincident_area(
                            self.hist_ifos, kwargs['time_addition'])
        # Volume is the allowed time difference window, multiplied by 2pi for
        # each phase difference dimension and by allowed range of SNR ratio
        # for each SNR ratio dimension : there are (n_ifos - 1) dimensions
        # for both phase and SNR
        n_ifos = len(self.hist_ifos)
        hist_vol = noise_twindow * \
            (2 * numpy.pi * (self.srbmax - self.srbmin) * self.swidth) ** \
            (n_ifos - 1)
        # Noise PDF is 1/volume, assuming a uniform distribution of noise
        # coincs
        logr_n = - numpy.log(hist_vol)

        # Combine to get final statistic: log of
        # ((rate of signals / rate of noise) * PTA Bayes factor)
        loglr = network_logvol - ln_noise_rate + logr_s - logr_n

        # cut off underflowing and very small values
        loglr[loglr < -30.] = -30.
        return loglr

    def coinc_lim_for_thresh(self, s, thresh, limifo,
                             **kwargs): # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest

        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples for all detectors
            except limifo.
        thresh: float
            The threshold on the coincident statistic.
        limifo: string
            The ifo for which the limit is to be found.

        Returns
        -------
        numpy.ndarray
            Array of limits on the limifo single statistic to
            exceed thresh.
        """

        # Safety against subclassing and not rethinking this
        allowed_names = ['ExpFitFgBgNormStatistic',
                         'ExpFitFgBgNormBBHStatistic',
                         'DQExpFitFgBgNormStatistic',
                         'ExpFitFgBgKDEStatistic']
        self._check_coinc_lim_subclass(allowed_names)

        if not self.has_hist:
            self.get_hist()
        # if the threshold is below this value all triggers will
        # pass because of rounding in the coinc method
        if thresh <= -30:
            return numpy.ones(len(s[0][1]['snglstat'])) * numpy.inf
        sngl_rates = {sngl[0]: sngl[1]['snglstat'] for sngl in s}
        # Add limifo to singles dict so that overlap time is calculated correctly
        sngl_rates[limifo] = numpy.zeros(len(s[0][1]))
        ln_noise_rate = coinc_rate.combination_noise_lograte(
                                  sngl_rates, kwargs['time_addition'])
        ln_noise_rate -= self.benchmark_lograte

        # Assume best case and use the maximum sigma squared from all triggers
        network_sigmasq = numpy.ones(len(s[0][1])) * kwargs['max_sigmasq']
        # Volume \propto sigma^3 or sigmasq^1.5
        network_logvol = 1.5 * numpy.log(network_sigmasq)
        # Get benchmark log volume as single-ifo information :
        # benchmark_logvol for a given template is not ifo-dependent, so
        # choose the first ifo for convenience
        benchmark_logvol = s[0][1]['benchmark_logvol']
        network_logvol -= benchmark_logvol

        # Assume best case scenario and use maximum signal rate
        logr_s = numpy.log(self.hist_max
                           * (kwargs['min_snr'] / self.ref_snr) ** -4.0)

        # Find total volume of phase-time-amplitude space occupied by noise
        # coincs
        # Extent of time-difference space occupied
        noise_twindow = coinc_rate.multiifo_noise_coincident_area(
                            self.hist_ifos, kwargs['time_addition'])
        # Volume is the allowed time difference window, multiplied by 2pi for
        # each phase difference dimension and by allowed range of SNR ratio
        # for each SNR ratio dimension : there are (n_ifos - 1) dimensions
        # for both phase and SNR
        n_ifos = len(self.hist_ifos)
        hist_vol = noise_twindow * \
            (2 * numpy.pi * (self.srbmax - self.srbmin) * self.swidth) ** \
            (n_ifos - 1)
        # Noise PDF is 1/volume, assuming a uniform distribution of noise
        # coincs
        logr_n = - numpy.log(hist_vol)

        loglr = - thresh + network_logvol - ln_noise_rate + logr_s - logr_n
        return loglr


class ExpFitFgBgNormBBHStatistic(ExpFitFgBgNormStatistic):
    """
    The ExpFitFgBgNormStatistic with a mass weighting factor.

    This is the same as the ExpFitFgBgNormStatistic except the likelihood
    is multiplied by a signal rate prior modelled as uniform over chirp mass.
    As templates are distributed roughly according to mchirp^(-11/3) we
    weight by the inverse of this. This ensures that loud events at high mass
    where template density is sparse are not swamped by events at lower masses
    where template density is high.
    """

    def __init__(self, sngl_ranking, files=None, ifos=None,
                 max_chirp_mass=None, **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.
        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.
        ifos: list of strs, not used here
            The list of detector names
        max_chirp_mass: float, default=None
            If given, if a template's chirp mass is above this value it will
            be reweighted as if it had this chirp mass. This is to avoid the
            problem where the distribution fails to be accurate at high mass
            and we can have a case where a single highest-mass template might
            produce *all* the loudest background (and foreground) events.
        """
        ExpFitFgBgNormStatistic.__init__(self, sngl_ranking, files=files,
                                         ifos=ifos, **kwargs)
        self.mcm = max_chirp_mass
        self.curr_mchirp = None

    def logsignalrate(self, stats, shift, to_shift):
        """
        Calculate the normalized log rate density of signals via lookup

        This calls back to the Parent class and then applies the chirp mass
        weighting factor.

        Parameters
        ----------
        stats: list of dicts giving single-ifo quantities, ordered as
            self.ifos
        shift: numpy array of float, size of the time shift vector for each
            coinc to be ranked
        to_shift: list of int, multiple of the time shift to apply ordered
            as self.ifos

        Returns
        -------
        value: log of coinc signal rate density for the given single-ifo
            triggers and time shifts
        """
        # Model signal rate as uniform over chirp mass, background rate is
        # proportional to mchirp^(-11/3) due to density of templates
        logr_s = ExpFitFgBgNormStatistic.logsignalrate(
                    self,
                    stats,
                    shift,
                    to_shift
                    )
        logr_s += numpy.log((self.curr_mchirp / 20.0) ** (11./3.0))
        return logr_s

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        In this case the ranking rescaled (see the lognoiserate method here)
        with the phase, end time, sigma, SNR, template_id and the
        benchmark_logvol values added in. This also stored the current chirp
        mass for use when computing the coinc statistic values.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        from pycbc.conversions import mchirp_from_mass1_mass2
        self.curr_mchirp = mchirp_from_mass1_mass2(trigs.param['mass1'],
                                                   trigs.param['mass2'])
        if self.mcm is not None:
            # Careful - input might be a str, so cast to float
            self.curr_mchirp = min(self.curr_mchirp, float(self.mcm))
        return ExpFitFgBgNormStatistic.single(self, trigs)

    def coinc_lim_for_thresh(self, s, thresh, limifo,
                             **kwargs): # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest

        Calculate the required single detector statistic to exceed
        the threshold for each of the input triggers.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples for all detectors
            except limifo.
        thresh: float
            The threshold on the coincident statistic.
        limifo: string
            The ifo for which the limit is to be found.

        Returns
        -------
        numpy.ndarray
            Array of limits on the limifo single statistic to
            exceed thresh.
        """

        loglr = ExpFitFgBgNormStatistic.coinc_lim_for_thresh(
                    self, s, thresh, limifo, **kwargs)
        loglr += numpy.log((self.curr_mchirp / 20.0) ** (11./3.0))
        return loglr


class ExpFitFgBgKDEStatistic(ExpFitFgBgNormStatistic):
    """
    The ExpFitFgBgNormStatistic with an additional mass and spin weighting
    factor determined by KDE statistic files.

    This is the same as the ExpFitFgBgNormStatistic except the likelihood
    ratio is multiplied by the ratio of signal KDE to template KDE over some
    parameters covering the bank.
    """

    def __init__(self, sngl_ranking, files=None, ifos=None, **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.
        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.
        ifos: list of strs, not used here
            The list of detector names
        """
        ExpFitFgBgNormStatistic.__init__(self, sngl_ranking, files=files,
                                         ifos=ifos, **kwargs)
        # The stat file attributes are hard-coded as 'signal-kde_file'
        # and 'template-kde_file'
        parsed_attrs = [f.split('-') for f in self.files.keys()]
        self.kde_names = [at[0] for at in parsed_attrs if
                       (len(at) == 2 and at[1] == 'kde_file')]
        assert sorted(self.kde_names) == ['signal', 'template'], \
            "Two stat files are required, they should have stat attr " \
            "'signal-kde_file' and 'template-kde_file' respectively"

        self.kde_by_tid = {}
        for kname in self.kde_names:
            self.assign_kdes(kname)
        # This will hold the template ids of the events for the statistic
        # calculation
        self.curr_tnum = None

    def assign_kdes(self, kname):
        """
        Extract values from KDE files

        Parameters
        -----------
        kname: str
            Used to label the kde files.
        """
        with h5py.File(self.files[kname+'-kde_file'], 'r') as kde_file:
            self.kde_by_tid[kname+'_kdevals'] = kde_file['data_kde'][:]

    def single(self, trigs):
        """
        Calculate the necessary single detector information including getting
        template ids from single detector triggers.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        try:
            # template_num exists if accessed via coinc_findtrigs
            self.curr_tnum = trigs.template_num
        except AttributeError:
            # exists for SingleDetTriggers
            self.curr_tnum = trigs['template_id']
        return ExpFitFgBgNormStatistic.single(self, trigs)

    def logsignalrate(self, stats, shift, to_shift):
        """
        Calculate the normalized log rate density of signals via lookup.

        This calls back to the parent class and then applies the ratio_kde
        weighting factor.

        Parameters
        ----------
        stats: list of dicts giving single-ifo quantities, ordered as
            self.ifos
        shift: numpy array of float, size of the time shift vector for each
            coinc to be ranked
        to_shift: list of int, multiple of the time shift to apply ordered
            as self.ifos

        Returns
        -------
        value: log of coinc signal rate density for the given single-ifo
            triggers and time shifts
        """
        logr_s = ExpFitFgBgNormStatistic.logsignalrate(self, stats, shift,
                                                       to_shift)
        signal_kde = self.kde_by_tid["signal_kdevals"][self.curr_tnum]
        template_kde = self.kde_by_tid["template_kdevals"][self.curr_tnum]
        logr_s += numpy.log(signal_kde / template_kde)
        return logr_s

    def coinc_lim_for_thresh(self, s, thresh, limifo, **kwargs):
        """
        Optimization function to identify coincs too quiet to be of interest

        Calculate the required single detector statistic to exceed the
        threshold for each of the input trigers.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples for all detectors
            except limifo.
        thresh: float
            The threshold on the coincident statistic.
        limifo: string
            The ifo for which the limit is to be found.

        Returns
        -------
        numpy.ndarray
            Array of limits on the limifo single statistic to
            exceed thresh.
        """
        loglr = ExpFitFgBgNormStatistic.coinc_lim_for_thresh(
            self, s, thresh, limifo, **kwargs)
        signal_kde = self.kde_by_tid["signal_kdevals"][self.curr_tnum]
        template_kde = self.kde_by_tid["template_kdevals"][self.curr_tnum]
        loglr += numpy.log(signal_kde / template_kde)
        return loglr


class DQExpFitFgBgNormStatistic(ExpFitFgBgNormStatistic):
    """
    The ExpFitFgBgNormStatistic with DQ-based reranking.

    This is the same as the ExpFitFgBgNormStatistic except the likelihood
    ratio is corrected via estimating relative noise trigger rates based on
    the DQ time series.
    """

    def __init__(self, sngl_ranking, files=None, ifos=None,
                 **kwargs):
        """
        Parameters
        ----------
        sngl_ranking: str
            The name of the ranking to use for the single-detector triggers.
        files: list of strs, needed here
            A list containing the filenames of hdf format files used to help
            construct the coincident statistics. The files must have a 'stat'
            attribute which is used to associate them with the appropriate
            statistic class.
        ifos: list of strs, not used here
            The list of detector names
        """
        ExpFitFgBgNormStatistic.__init__(self, sngl_ranking, files=files,
                                         ifos=ifos, **kwargs)
        self.dq_val_by_time = {}
        self.dq_bin_by_id = {}
        for k in self.files.keys():
            parsed_attrs = k.split('-')
            if len(parsed_attrs) < 3:
                continue
            if parsed_attrs[2] == 'dq_ts_reference':
                ifo = parsed_attrs[0]
                dq_type = parsed_attrs[1]
                dq_vals = self.assign_dq_val(k)
                dq_bins = self.assign_bin_id(k)
                if ifo not in self.dq_val_by_time:
                    self.dq_val_by_time[ifo] = {}
                    self.dq_bin_by_id[ifo] = {}
                self.dq_val_by_time[ifo][dq_type] = dq_vals
                self.dq_bin_by_id[ifo][dq_type] = dq_bins

    def assign_bin_id(self, key):
        """
        Assign bin ID values
        Assign each template id to a bin name based on a
        referenced statistic file.

        Parameters
        ----------
        key: str
            statistic file key string

        Returns
        ---------
        bin_dict: dict of strs
            Dictionary containing the bin name for each template id
        """
        ifo = key.split('-')[0]
        with h5py.File(self.files[key], 'r') as dq_file:
            bin_names = dq_file.attrs['names'][:]
            locs = []
            names = []
            for bin_name in bin_names:
                bin_locs = dq_file[ifo + '/locs/' + bin_name][:]
                locs = list(locs)+list(bin_locs.astype(int))
                names = list(names)+list([bin_name]*len(bin_locs))

        bin_dict = dict(zip(locs, names))
        return bin_dict

    def assign_dq_val(self, key):
        """
        Assign dq values to each time for every bin based on a
        referenced statistic file.

        Parameters
        ----------
        key: str
            statistic file key string

        Returns
        ---------
        dq_dict: dict of {time: dq_value} dicts for each bin
            Dictionary containing the mapping between the time
            and the dq value for each individual bin.

        """
        ifo = key.split('-')[0]
        with h5py.File(self.files[key], 'r') as dq_file:
            times = dq_file[ifo+'/times'][:]
            bin_names = dq_file.attrs['names'][:]
            dq_dict = {}
            for bin_name in bin_names:
                dq_vals = dq_file[ifo+'/dq_vals/'+bin_name][:]
                dq_dict[bin_name] = dict(zip(times, dq_vals))

        return dq_dict

    def find_dq_val(self, trigs):
        """Get dq values for a specific ifo and times"""
        time = trigs['end_time'].astype(int)
        try:
            tnum = trigs.template_num
            ifo = trigs.ifo
        except AttributeError:
            tnum = trigs['template_id']
            assert len(self.ifos) == 1
            # Should be exactly one ifo provided
            ifo = self.ifos[0]
        dq_val = numpy.zeros(len(time))
        if ifo in self.dq_val_by_time:
            for (i, t) in enumerate(time):
                for k in self.dq_val_by_time[ifo].keys():
                    if isinstance(tnum, numpy.ndarray):
                        bin_name = self.dq_bin_by_id[ifo][k][tnum[i]]
                    else:
                        bin_name = self.dq_bin_by_id[ifo][k][tnum]
                    val = self.dq_val_by_time[ifo][k][bin_name][int(t)]
                    dq_val[i] = max(dq_val[i], val)
        return dq_val

    def lognoiserate(self, trigs):
        """
        Calculate the log noise rate density over single-ifo ranking

        Read in single trigger information, compute the ranking
        and rescale by the fitted coefficients alpha and rate

        Parameters
        -----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        ---------
        lognoisel: numpy.array
            Array of log noise rate density for each input trigger.
        """
        logr_n = ExpFitFgBgNormStatistic.lognoiserate(
                    self, trigs)
        logr_n += self.find_dq_val(trigs)
        return logr_n


statistic_dict = {
    'quadsum': QuadratureSumStatistic,
    'single_ranking_only': QuadratureSumStatistic,
    'phasetd': PhaseTDStatistic,
    'exp_fit_stat': ExpFitStatistic,
    'exp_fit_csnr': ExpFitCombinedSNR,
    'phasetd_exp_fit_stat': PhaseTDExpFitStatistic,
    'dq_phasetd_exp_fit_fgbg_norm': DQExpFitFgBgNormStatistic,
    'exp_fit_bg_rate': ExpFitBgRateStatistic,
    'phasetd_exp_fit_fgbg_norm': ExpFitFgBgNormStatistic,
    'phasetd_exp_fit_fgbg_bbh_norm': ExpFitFgBgNormBBHStatistic,
    'phasetd_exp_fit_fgbg_kde': ExpFitFgBgKDEStatistic,
}


def get_statistic(stat):
    """
    Error-handling sugar around dict lookup for coincident statistics

    Parameters
    ----------
    stat : string
        Name of the coincident statistic

    Returns
    -------
    class
        Subclass of Stat base class

    Raises
    ------
    RuntimeError
        If the string is not recognized as corresponding to a Stat subclass
    """
    try:
        return statistic_dict[stat]
    except KeyError:
        raise RuntimeError('%s is not an available detection statistic' % stat)


def insert_statistic_option_group(parser, default_ranking_statistic=None):
    """
    Add ranking statistic options to the optparser object.

    Adds the options used to initialize a PyCBC Stat class.

    Parameters
    -----------
    parser : object
        OptionParser instance.
    default_ranking_statisic : str
        Allows setting a default statistic for the '--ranking-statistic'
        option. The option is no longer required if a default is provided.

    Returns
    --------
    strain_opt_group : optparser.argument_group
        The argument group that is added to the parser.
    """

    statistic_opt_group = parser.add_argument_group(
        "Options needed to initialize a PyCBC Stat class for computing the "
        "ranking of events from a PyCBC search."
    )

    statistic_opt_group.add_argument(
        "--ranking-statistic",
        default=default_ranking_statistic,
        choices=statistic_dict.keys(),
        required=True if default_ranking_statistic is None else False,
        help="The coinc ranking statistic to calculate"
    )

    statistic_opt_group.add_argument(
        "--sngl-ranking",
        choices=ranking.sngls_ranking_function_dict.keys(),
        required=True,
        help="The single-detector trigger ranking to use."
    )

    statistic_opt_group.add_argument(
        "--statistic-files",
        nargs='*',
        action='append',
        default=[],
        help="Files containing ranking statistic info"
    )

    statistic_opt_group.add_argument(
        "--statistic-keywords",
        nargs='*',
        default=[],
        help="Provide additional key-word arguments to be sent to "
             "the statistic class when it is initialized. Should "
             "be given in format --statistic-keywords "
             "KWARG1:VALUE1 KWARG2:VALUE2 KWARG3:VALUE3 ..."
    )

    return statistic_opt_group


def parse_statistic_keywords_opt(stat_kwarg_list):
    """
    Parse the list of statistic keywords into an appropriate dictionary.

    Take input from the input argument ["KWARG1:VALUE1", "KWARG2:VALUE2",
    "KWARG3:VALUE3"] and convert into a dictionary.

    Parameters
    ----------
    stat_kwarg_list : list
        Statistic keywords in list format

    Returns
    -------
    stat_kwarg_dict : dict
        Statistic keywords in dict format
    """
    stat_kwarg_dict = {}
    for inputstr in stat_kwarg_list:
        try:
            key, value = inputstr.split(':')
            stat_kwarg_dict[key] = value
        except ValueError:
            err_txt = "--statistic-keywords must take input in the " \
                      "form KWARG1:VALUE1 KWARG2:VALUE2 KWARG3:VALUE3 ... " \
                      "Received {}".format(' '.join(stat_kwarg_list))
            raise ValueError(err_txt)

    return stat_kwarg_dict


def get_statistic_from_opts(opts, ifos):
    """
    Return a Stat class from an optparser object.

    This will assume that the options in the statistic_opt_group are present
    and will use these options to call stat.get_statistic and initialize the
    appropriate Stat subclass with appropriate kwargs.

    Parameters
    ----------
    opts : optparse.OptParser instance
        The command line options
    ifos : list
        The list of detector names

    Returns
    -------
    class
        Subclass of Stat base class
    """
    # Allow None inputs
    if opts.statistic_files is None:
        opts.statistic_files = []
    if opts.statistic_keywords is None:
        opts.statistic_keywords = []

    # flatten the list of lists of filenames to a single list (may be empty)
    opts.statistic_files = sum(opts.statistic_files, [])

    extra_kwargs = parse_statistic_keywords_opt(opts.statistic_keywords)

    stat_class = get_statistic(opts.ranking_statistic)(
        opts.sngl_ranking,
        opts.statistic_files,
        ifos=ifos,
        **extra_kwargs
    )

    return stat_class
