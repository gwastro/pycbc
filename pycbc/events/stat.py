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
from . import ranking
from . import coinc_rate


class Stat(object):
    """Base class which should be extended to provide a coincident statistic"""

    def __init__(self, files=None, ifos=None, **kwargs):
        """Create a statistic class instance

        Parameters
        ----------
        files: list of strs
            A list containing the filenames of hdf format files used to help
        construct the coincident statistics. The files must have a 'stat'
        attribute which is used to associate them with the appropriate
        statistic class.

        ifos: list of detector names, optional
        """
        import h5py

        self.files = {}
        files = files or []
        for filename in files:
            f = h5py.File(filename, 'r')
            stat = (f.attrs['stat']).decode()
            if stat in self.files:
                raise RuntimeError("We already have one file with stat attr ="
                                   " %s. Can't provide more than one!" % stat)
            logging.info("Found file %s for stat %s", filename, stat)
            self.files[stat] = f

        # Provide the dtype of the single detector method's output
        # This is used by background estimation codes that need to maintain
        # a buffer of such values.
        self.single_dtype = numpy.float32
        # True if a larger single detector statistic will produce a larger
        # coincident statistic
        self.single_increasing = True

        self.ifos = ifos or []


class NewSNRStatistic(Stat):
    """Calculate the NewSNR coincident detection statistic"""

    def single(self, trigs):
        """Calculate the single detector statistic, here equal to newsnr

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
            Dictionary-like object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return ranking.get_newsnr(trigs)

    def coinc(self, s0, s1, slide, step): # pylint:disable=unused-argument
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s0: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s1: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: (unused in this statistic)
        step: (unused in this statistic)

        Returns
        -------
        numpy.ndarray
            Array of coincident ranking statistic values
        """
        return (s0 ** 2. + s1 ** 2.) ** 0.5

    def coinc_lim_for_thresh(self, s0, thresh):
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
        s1 = thresh ** 2. - s0 ** 2.
        s1[s1 < 0] = 0
        return s1 ** 0.5

    def coinc_multiifo(self, s, slide, step, to_shift,
                       **kwargs): # pylint:disable=unused-argument
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s: list
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
        return sum(sngl[1] ** 2. for sngl in s) ** 0.5

    def coinc_multiifo_lim_for_thresh(self, s, thresh, limifo,
                                      **kwargs): # pylint:disable=unused-argument
        """Calculate the required single detector statistic to exceed
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
        s0 = thresh ** 2. - sum(sngl[1] ** 2. for sngl in s)
        s0[s0 < 0] = 0
        return s0 ** 0.5


class NewSNRSGStatistic(NewSNRStatistic):
    """Calculate the NewSNRSG coincident detection statistic"""

    def single(self, trigs):
        """Calculate the single detector statistic, here equal to newsnr_sgveto

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
            Dictionary-like object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return ranking.get_newsnr_sgveto(trigs)


class NewSNRSGPSDStatistic(NewSNRSGStatistic):
    """Calculate the NewSNRSGPSD coincident detection statistic"""

    def single(self, trigs):
        """Calculate the single detector statistic, here equal to newsnr
        combined with sgveto and psdvar statistic

        Parameters
        ----------
        trigs: dict of numpy.ndarrays

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return ranking.get_newsnr_sgveto_psdvar(trigs)


class NewSNRSGPSDThresholdStatistic(NewSNRSGStatistic):
    """Calculate the NewSNRSGPSD coincident detection statistic"""

    def single(self, trigs):
        """Calculate the single detector statistic, here equal to newsnr
        combined with sgveto and psdvar statistic

        Parameters
        ----------
        trigs: dict of numpy.ndarrays

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return ranking.get_newsnr_sgveto_psdvar_threshold(trigs)


class NewSNRSGPSDScaledStatistic(NewSNRSGStatistic):
    """Calculate the NewSNRSGPSD coincident detection statistic"""

    def single(self, trigs):
        """Calculate the single detector statistic, here equal to newsnr
        combined with sgveto and psdvar statistic

        Parameters
        ----------
        trigs: dict of numpy.ndarrays

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return ranking.get_newsnr_sgveto_psdvar_scaled(trigs)


class NewSNRSGPSDScaledThresholdStatistic(NewSNRSGStatistic):
    """Calculate the NewSNRSGPSD coincident detection statistic"""

    def single(self, trigs):
        """Calculate the single detector statistic, here equal to newsnr
        combined with sgveto and psdvar statistic

        Parameters
        ----------
        trigs: dict of numpy.ndarrays

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return ranking.get_newsnr_sgveto_psdvar_scaled_threshold(trigs)


class NetworkSNRStatistic(NewSNRStatistic):
    """Same as the NewSNR statistic, but just sum of squares of SNRs"""

    def single(self, trigs):
        return trigs['snr']


class NewSNRCutStatistic(NewSNRStatistic):
    """Same as the NewSNR statistic, but demonstrates a cut of the triggers"""

    def single(self, trigs):
        """Calculate the single detector statistic.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
            Dictionary-like object holding single detector trigger information.

        Returns
        -------
        newsnr: numpy.ndarray
            Array of single detector values
        """
        newsnr = ranking.get_newsnr(trigs)
        rchisq = trigs['chisq'][:] / (2. * trigs['chisq_dof'][:] - 2.)
        newsnr[numpy.logical_and(newsnr < 10, rchisq > 2)] = -1
        return newsnr

    def coinc(self, s0, s1, slide, step): # pylint:disable=unused-argument
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s0: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s1: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: (unused in this statistic)
        step: (unused in this statistic)

        Returns
        -------
        cstat: numpy.ndarray
            Array of coincident ranking statistic values
        """
        cstat = (s0 ** 2. + s1 ** 2.) ** 0.5
        cstat[s0 == -1] = 0
        cstat[s1 == -1] = 0
        return cstat

    def coinc_lim_for_thresh(self, s0, thresh):
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
        s1 = thresh ** 2. - s0 ** 2.
        s1[s0 == -1] = numpy.inf
        s1[s1 < 0] = 0
        return s1 ** 0.5


class PhaseTDNewStatistic(NewSNRStatistic):
    """Statistic that re-weights combined newsnr using coinc parameters.

    The weighting is based on the PDF of time delays, phase differences and
    amplitude ratios between triggers in different ifos.
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        NewSNRStatistic.__init__(self, files=files, ifos=ifos, **kwargs)

        self.single_dtype = [('snglstat', numpy.float32),
                             ('coa_phase', numpy.float32),
                             ('end_time', numpy.float64),
                             ('sigmasq', numpy.float32),
                             ('snr', numpy.float32)
                             ]

        # Assign attribute so that it can be replaced with other functions
        self.get_newsnr = ranking.get_newsnr
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

    def get_hist(self, ifos=None):
        """Read in a signal density file for the ifo combination"""

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

                match = [ifo in name for ifo in ifos]
                if False in match:
                    continue
                else:
                    selected = name
                    break

        if selected is None:
            raise RuntimeError("Couldn't figure out which stat file to use")

        logging.info("Using signal histogram %s for ifos %s", selected, ifos)
        histfile = self.files[selected]
        self.hist_ifos = histfile.attrs['ifos']
        n_ifos = len(self.hist_ifos)

        # Histogram bin attributes
        self.twidth = histfile.attrs['twidth']
        self.pwidth = histfile.attrs['pwidth']
        self.swidth = histfile.attrs['swidth']
        self.srbmin = histfile.attrs['srbmin']
        self.srbmax = histfile.attrs['srbmax']

        bin_volume = (self.twidth * self.pwidth * self.swidth) ** (n_ifos - 1)
        self.hist_max = - 1. * numpy.inf

        # Read histogram for each ifo, to use if that ifo has smallest SNR in
        # the coinc
        for ifo in self.hist_ifos:

            weights = histfile[ifo]['weights'][:]
            # renormalise to PDF
            self.weights[ifo] = weights / (weights.sum() * bin_volume)

            param = histfile[ifo]['param_bin'][:]

            if param.dtype == numpy.int8:
                # Older style, incorrectly sorted histogram file
                ncol = param.shape[1]
                self.pdtype = [('c%s' % i, param.dtype) for i in range(ncol)]
                self.param_bin[ifo] = numpy.zeros(len(self.weights[ifo]),
                                                  dtype=self.pdtype)
                for i in range(ncol):
                    self.param_bin[ifo]['c%s' % i] = param[:, i]

                lsort = self.param_bin[ifo].argsort()
                self.param_bin[ifo] = self.param_bin[ifo][lsort]
                self.weights[ifo] = self.weights[ifo][lsort]
            else:
                # New style, efficient histogram file
                # param bin and weights have already been sorted
                self.param_bin[ifo] = param
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

                self.c0_size[ifo] = 2 * (abs(self.param_bin[ifo]['c0']).max() + 1)
                self.c1_size[ifo] = 2 * (abs(self.param_bin[ifo]['c1']).max() + 1)
                self.c2_size[ifo] = 2 * (abs(self.param_bin[ifo]['c2']).max() + 1)

                array_size = [self.c0_size[ifo], self.c1_size[ifo],
                              self.c2_size[ifo]]
                dtypec = self.weights[ifo].dtype
                self.two_det_weights[ifo] = \
                    numpy.zeros(array_size, dtype=dtypec) + self.max_penalty
                id0 = self.param_bin[ifo]['c0'].astype(numpy.int32) + self.c0_size[ifo] // 2
                id1 = self.param_bin[ifo]['c1'].astype(numpy.int32) + self.c1_size[ifo] // 2
                id2 = self.param_bin[ifo]['c2'].astype(numpy.int32) + self.c2_size[ifo] // 2
                self.two_det_weights[ifo][id0, id1, id2] = self.weights[ifo]

        relfac = histfile.attrs['sensitivity_ratios']
        for ifo, sense in zip(self.hist_ifos, relfac):
            self.relsense[ifo] = sense

        self.has_hist = True

    def single(self, trigs):
        """Calculate the single detector statistic & assemble other parameters

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information. 'snr', 'chisq',
        'chisq_dof', 'coa_phase', 'end_time', and 'sigmasq' are required keys.

        Returns
        -------
        numpy.ndarray
            Array of single detector parameter values
        """
        sngl_stat = self.get_newsnr(trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase'][:]
        singles['end_time'] = trigs['end_time'][:]
        singles['sigmasq'] = trigs['sigmasq'][:]
        singles['snr'] = trigs['snr'][:]
        return numpy.array(singles, ndmin=1)

    def logsignalrate(self, s0, s1, shift):
        to_shift = [-1, 0]
        stats = {self.ifos[0]: s0, self.ifos[1]: s1}
        return self.logsignalrate_multiifo(stats, shift, to_shift)

    def logsignalrate_multiifo(self, stats, shift, to_shift):
        """Calculate the normalized log rate density of signals via lookup

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
        # Convert time shift vector to dict, as hist ifos and self.ifos may
        # not be in same order
        to_shift = {ifo: s for ifo, s in zip(self.ifos, to_shift)}

        if not self.has_hist:
            self.get_hist()

        # Figure out which ifo of the contributing ifos has the smallest SNR,
        # to use as reference for choosing the signal histogram.
        snrs = numpy.array([numpy.array(stats[ifo]['snr'], ndmin=1)
                           for ifo in self.ifos])
        smin = numpy.argmin(snrs, axis=0)
        # Store a list of the triggers using each ifo as reference
        rtypes = {ifo: numpy.where(smin == j)[0]
                  for j, ifo in enumerate(self.ifos)}

        # Get reference ifo information
        rate = numpy.zeros(len(shift), dtype=numpy.float32)
        for ref_ifo in self.ifos:
            rtype = rtypes[ref_ifo]
            ref = stats[ref_ifo]
            pref = numpy.array(ref['coa_phase'], ndmin=1)[rtype]
            tref = numpy.array(ref['end_time'], ndmin=1)[rtype]
            sref = numpy.array(ref['snr'], ndmin=1)[rtype]
            sigref = numpy.array(ref['sigmasq'], ndmin=1) ** 0.5
            sigref = sigref[rtype]
            senseref = self.relsense[self.hist_ifos[0]]

            binned = []
            other_ifos = [ifo for ifo in self.ifos if ifo != ref_ifo]
            for ifo in other_ifos:
                sc = stats[ifo]
                p = numpy.array(sc['coa_phase'], ndmin=1)[rtype]
                t = numpy.array(sc['end_time'], ndmin=1)[rtype]
                s = numpy.array(sc['snr'], ndmin=1)[rtype]

                sense = self.relsense[ifo]
                sig = numpy.array(sc['sigmasq'], ndmin=1) ** 0.5
                sig = sig[rtype]

                # Calculate differences
                pdif = (pref - p) % (numpy.pi * 2.0)
                tdif = shift[rtype] * to_shift[ref_ifo] + \
                    tref - shift[rtype] * to_shift[ifo] - t
                sdif = s / sref * sense / senseref * sigref / sig

                # Put into bins
                tbin = (tdif / self.twidth).astype(numpy.int)
                pbin = (pdif / self.pwidth).astype(numpy.int)
                sbin = (sdif / self.swidth).astype(numpy.int)
                binned += [tbin, pbin, sbin]

            # Convert binned to same dtype as stored in hist
            nbinned = numpy.zeros(len(pbin), dtype=self.pdtype)
            for i, b in enumerate(binned):
                nbinned['c%s' % i] = b

            # Read signal weight from precalculated histogram
            if self.two_det_flag:
                # High-RAM, low-CPU option for two-det
                rate[rtype] = numpy.zeros(len(nbinned)) + self.max_penalty

                id0 = nbinned['c0'].astype(numpy.int32) + self.c0_size[ref_ifo] // 2
                id1 = nbinned['c1'].astype(numpy.int32) + self.c1_size[ref_ifo] // 2
                id2 = nbinned['c2'].astype(numpy.int32) + self.c2_size[ref_ifo] // 2

                # look up keys which are within boundaries
                within = (id0 > 0) & (id0 < self.c0_size[ref_ifo])
                within = within & (id1 > 0) & (id1 < self.c1_size[ref_ifo])
                within = within & (id2 > 0) & (id2 < self.c2_size[ref_ifo])
                within = numpy.where(within)[0]
                rate[rtype[within]] = self.two_det_weights[ref_ifo][id0[within], id1[within], id2[within]]
            else:
                # Low[er]-RAM, high[er]-CPU option for >two det
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
            rate[rtype] *= (sref / self.ref_snr) ** -4.0

        return numpy.log(rate)


class PhaseTDStatistic(NewSNRStatistic):
    """Statistic that re-weights combined newsnr using coinc parameters.

    The weighting is based on the PDF of time delays, phase differences and
    amplitude ratios between triggers in different ifos.
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        NewSNRStatistic.__init__(self, files=files, ifos=ifos, **kwargs)

        self.single_dtype = [('snglstat', numpy.float32),
                             ('coa_phase', numpy.float32),
                             ('end_time', numpy.float64),
                             ('sigmasq', numpy.float32),
                             ('snr', numpy.float32)
                             ]

        # Assign attribute so that it can be replaced with other functions
        self.get_newsnr = ranking.get_newsnr

        self.hist = None
        self.bins = {}
        self.hist_ifos = []

    def get_hist(self, ifos=None, norm='max'):
        """Read in a signal density file for the ifo combination"""

        # default name for old 2-ifo workflow
        if 'phasetd_newsnr' in self.files:
            histfile = self.files['phasetd_newsnr']
        else:
            ifos = ifos or self.ifos  # if None, use the instance attribute
            if len(ifos) != 2:
                raise RuntimeError("Need exactly 2 ifos for the p/t/a "
                                   "statistic! Ifos given were " + ifos)
            matching = [k for k in self.files.keys() if \
                        'phasetd' in k and (ifos[0] in k and ifos[1] in k)]
            if len(matching) == 1:
                histfile = self.files[matching[0]]
            else:
                raise RuntimeError(
                  "%i statistic files had an attribute matching phasetd*%s%s !"
                  "Should be exactly 1" % (len(matching), ifos[0], ifos[1]))
            logging.info("Using signal histogram %s for ifos %s", matching,
                         ifos)

        self.hist = histfile['map'][:]
        self.hist_ifos = ifos

        if norm == 'max':
            # Normalize so that peak of hist is equal to unity
            self.hist = self.hist / float(self.hist.max())
            self.hist = numpy.log(self.hist)
        else:
            raise NotImplementedError("Sorry, we have no other normalizations")

        # Bin boundaries are stored in the hdf file
        self.bins['dt'] = histfile['tbins'][:]
        self.bins['dphi'] = histfile['pbins'][:]
        self.bins['snr'] = histfile['sbins'][:]
        self.bins['sigma_ratio'] = histfile['rbins'][:]
        self.hist_max = self.hist.max()

    def single(self, trigs):
        """Calculate the single detector statistic & assemble other parameters

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information. 'snr', 'chisq',
        'chisq_dof', 'coa_phase', 'end_time', and 'sigmasq' are required keys.

        Returns
        -------
        numpy.ndarray
            Array of single detector parameter values
        """
        sngl_stat = self.get_newsnr(trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase'][:]
        singles['end_time'] = trigs['end_time'][:]
        singles['sigmasq'] = trigs['sigmasq'][:]
        singles['snr'] = trigs['snr'][:]
        return numpy.array(singles, ndmin=1)

    def signal_hist(self, td, pd, sn0, sn1, rd):
        assert self.hist is not None

        # enforce that sigma ratio is < 1 by swapping values
        snr0 = sn0 * 1
        snr1 = sn1 * 1

        snr0[rd > 1] = sn1[rd > 1]
        snr1[rd > 1] = sn0[rd > 1]
        rd[rd > 1] = 1. / rd[rd > 1]

        # Find which bin each coinc falls into
        tv = numpy.searchsorted(self.bins['dt'], td) - 1
        pv = numpy.searchsorted(self.bins['dphi'], pd) - 1
        s0v = numpy.searchsorted(self.bins['snr'], snr0) - 1
        s1v = numpy.searchsorted(self.bins['snr'], snr1) - 1
        rv = numpy.searchsorted(self.bins['sigma_ratio'], rd) - 1

        # Enforce that points fit into the bin boundaries: if a point lies
        # outside the boundaries it is pushed back to the nearest bin.
        for binnum, axis in zip([tv, pv, rv, s0v, s1v],
                                ['dt', 'dphi', 'sigma_ratio', 'snr', 'snr']):
            binend = len(self.bins[axis])
            binnum[binnum < 0] = 0
            binnum[binnum >= binend - 1] = binend - 2

        return self.hist[tv, pv, s0v, s1v, rv]

    def slide_dt(self, singles, shift, slide_vec):
        # Apply time shifts in the multiples specified by slide_vec
        # and return resulting time difference
        assert len(singles) == 2
        assert len(slide_vec) == 2
        dt = singles[0]['end_time'] + shift * slide_vec[0] -\
            (singles[1]['end_time'] + shift * slide_vec[1])
        return dt

    def logsignalrate(self, s0, s1, shift):
        """Calculate the normalized log rate density of signals via lookup"""

        # does not require ifos to be specified, only 1 p/t/a file
        if self.hist is None:
            self.get_hist()

        # for 2-ifo pipeline, add time shift to 2nd ifo ('s1')
        slidevec = [0, 1]
        td = numpy.array(self.slide_dt([s0, s1], shift, slidevec),
                         ndmin=1)
        if numpy.any(td > 1.):
            raise RuntimeError(
              "Time difference bigger than 1 second after applying any time "
              "shifts! This should not happen")
        pd = numpy.array((s0['coa_phase'] - s1['coa_phase']) % \
                         (2. * numpy.pi), ndmin=1)
        sn0 = numpy.array(s0['snr'], ndmin=1)
        sn1 = numpy.array(s1['snr'], ndmin=1)
        rd = numpy.array((s0['sigmasq'] / s1['sigmasq']) ** 0.5, ndmin=1)

        return self.signal_hist(td, pd, sn0, sn1, rd)

    def logsignalrate_multiifo(self, s, shift, to_shift):
        """
        Parameters
        ----------
        s: list, length 2
            List of sets of single-ifo trigger parameter values
        shift: numpy.ndarray
            Array of floats giving the time shifts to be applied with
        multiples given by to_shift
        to_shift: list, length 2
            List of time shift multiples
        """
        assert len(s) == 2
        assert len(to_shift) == 2

        # At present for triples use the H/L signal histogram
        hist_ifos = self.ifos if len(self.ifos) == 2 else ['H1', 'L1']
        if self.hist is None:
            self.get_hist(hist_ifos)
        else:
            assert self.hist_ifos == hist_ifos
            logging.info("Using pre-set signal histogram for %s",
                         self.hist_ifos)

        td = self.slide_dt(s, shift, to_shift)
        if numpy.any(td > 1.):
            raise RuntimeError(
              "Time difference bigger than 1 second after applying any time "
              "shifts! This should not happen")
        pd = numpy.array((s[0]['coa_phase'] - s[1]['coa_phase']) % \
                         (2. * numpy.pi), ndmin=1)
        sn0 = numpy.array(s[0]['snr'], ndmin=1)
        sn1 = numpy.array(s[1]['snr'], ndmin=1)
        rd = numpy.array((s[0]['sigmasq'] / s[1]['sigmasq']) ** 0.5, ndmin=1)

        return self.signal_hist(td, pd, sn0, sn1, rd)

    def coinc(self, s0, s1, slide, step):
        """Calculate the coincident detection statistic.
        Parameters
        ----------
        s0: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s1: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: numpy.ndarray
            Array of ints. These represent the multiple of the timeslide
        interval to bring a pair of single detector triggers into coincidence.
        step: float
            The timeslide interval in seconds.
        Returns
        -------
        coinc_stat: numpy.ndarray
            An array of the coincident ranking statistic values
        """
        rstat = s0['snglstat'] ** 2. + s1['snglstat'] ** 2.
        cstat = rstat + 2. * self.logsignalrate(s0, s1, slide * step)
        cstat[cstat < 0] = 0
        return cstat ** 0.5

    def coinc_lim_for_thresh(self, s0, thresh):
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
        if self.hist is None:
            self.get_hist()
        s1 = thresh ** 2. - s0['snglstat'] ** 2.
        # Assume best case scenario and use maximum signal rate
        s1 -= 2. * self.hist_max
        s1[s1 < 0] = 0
        return s1 ** 0.5


class PhaseTDSGStatistic(PhaseTDStatistic):
    """PhaseTDStatistic but with sine-Gaussian veto added to the

    single-detector ranking
    """
    def __init__(self, files=None, ifos=None, **kwargs):
        PhaseTDStatistic.__init__(self, files=files, ifos=ifos, **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto


class ExpFitStatistic(NewSNRStatistic):
    """Detection statistic using an exponential falloff noise model.

    Statistic approximates the negative log noise coinc rate density per
    template over single-ifo newsnr values.
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        if not len(files):
            raise RuntimeError("Can't find any statistic files !")
        NewSNRStatistic.__init__(self, files=files, ifos=ifos, **kwargs)

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

        self.get_newsnr = ranking.get_newsnr
        self.single_increasing = False

    def assign_fits(self, ifo):
        coeff_file = self.files[ifo+'-fit_coeffs']
        template_id = coeff_file['template_id'][:]
        alphas = coeff_file['fit_coeff'][:]
        rates = coeff_file['count_above_thresh'][:]
        # the template_ids and fit coeffs are stored in an arbitrary order
        # create new arrays in template_id order for easier recall
        tid_sort = numpy.argsort(template_id)
        return {'alpha': alphas[tid_sort],
                'rate': rates[tid_sort],
                'thresh': coeff_file.attrs['stat_threshold']
                }

    def get_ref_vals(self, ifo):
        self.alphamax[ifo] = self.fits_by_tid[ifo]['alpha'].max()

    def find_fits(self, trigs):
        """Get fit coeffs for a specific ifo and template id(s)"""

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
        alphai = self.fits_by_tid[ifo]['alpha'][tnum]
        ratei = self.fits_by_tid[ifo]['rate'][tnum]
        thresh = self.fits_by_tid[ifo]['thresh']
        return alphai, ratei, thresh

    def lognoiserate(self, trigs):
        """Calculate the log noise rate density over single-ifo newsnr

        Read in single trigger information, make the newsnr statistic
        and rescale by the fitted coefficients alpha and rate
        """
        alphai, ratei, thresh = self.find_fits(trigs)
        newsnr = self.get_newsnr(trigs)
        # alphai is constant of proportionality between single-ifo newsnr and
        #   negative log noise likelihood in given template
        # ratei is rate of trigs in given template compared to average
        # thresh is stat threshold used in given ifo
        lognoisel = - alphai * (newsnr - thresh) + numpy.log(alphai) + \
                      numpy.log(ratei)
        return numpy.array(lognoisel, ndmin=1, dtype=numpy.float32)

    def single(self, trigs):
        """Single-detector statistic, here just equal to the log noise rate"""

        return self.lognoiserate(trigs)

    def coinc(self, s0, s1, slide, step): # pylint:disable=unused-argument
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

    def coinc_lim_for_thresh(self, s0, thresh):
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
    """Reworking of ExpFitStatistic designed to resemble network SNR

    Use a monotonic function of the negative log noise rate density which
    approximates combined (new)snr for coincs with similar newsnr in each ifo
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        ExpFitStatistic.__init__(self, files=files, ifos=ifos, **kwargs)
        # for low-mass templates the exponential slope alpha \approx 6
        self.alpharef = 6.
        self.single_increasing = True

    def use_alphamax(self):
        # take reference slope as the harmonic mean of individual ifo slopes
        inv_alphas = [1. / self.alphamax[i] for i in self.bg_ifos]
        self.alpharef = 1. / (sum(inv_alphas) / len(inv_alphas))

    def single(self, trigs):
        logr_n = self.lognoiserate(trigs)
        _, _, thresh = self.find_fits(trigs)
        # shift by log of reference slope alpha
        logr_n += -1. * numpy.log(self.alpharef)
        # add threshold and rescale by reference slope
        stat = thresh - (logr_n / self.alpharef)
        return numpy.array(stat, ndmin=1, dtype=numpy.float32)

    def single_multiifo(self, s):
        if self.single_increasing:
            sngl_multiifo = s[1]['snglstat']
        else:
            sngl_multiifo = -1.0 * s[1]['snglstat']
        return sngl_multiifo

    def coinc(self, s0, s1, slide, step): # pylint:disable=unused-argument
        # scale by 1/sqrt(2) to resemble network SNR
        return (s0 + s1) / 2.**0.5

    def coinc_lim_for_thresh(self, s0, thresh):
        return thresh * (2. ** 0.5) - s0

    def coinc_multiifo(self, s, slide, step, to_shift,
                       **kwargs): # pylint:disable=unused-argument
        # scale by 1/sqrt(number of ifos) to resemble network SNR
        return sum(sngl[1] for sngl in s) / len(s)**0.5

    def coinc_multiifo_lim_for_thresh(self, s, thresh,
                                      limifo, **kwargs): # pylint:disable=unused-argument
        return thresh * ((len(s) + 1) ** 0.5) - sum(sngl[1] for sngl in s)


class ExpFitSGCombinedSNR(ExpFitCombinedSNR):
    """ExpFitCombinedSNR but with sine-Gaussian veto added to the single

    detector ranking
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        ExpFitCombinedSNR.__init__(self, files=files, ifos=ifos, **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto


class ExpFitSGPSDCombinedSNR(ExpFitCombinedSNR):
    """ExpFitCombinedSNR but with sine-Gaussian veto and PSD variation added to

    the single detector ranking
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        ExpFitCombinedSNR.__init__(self, files=files, ifos=ifos, **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar


class PhaseTDExpFitStatistic(PhaseTDStatistic, ExpFitCombinedSNR):
    """Statistic combining exponential noise model with signal histogram PDF"""

    # default is 2-ifo operation with exactly 1 'phasetd' file
    def __init__(self, files=None, ifos=None, **kwargs):
        # read in both foreground PDF and background fit info
        ExpFitCombinedSNR.__init__(self, files=files, ifos=ifos, **kwargs)
        # need the self.single_dtype value from PhaseTDStatistic
        PhaseTDStatistic.__init__(self, files=files, ifos=ifos, **kwargs)

    def single(self, trigs):
        # same single-ifo stat as ExpFitCombinedSNR
        sngl_stat = ExpFitCombinedSNR.single(self, trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase'][:]
        singles['end_time'] = trigs['end_time'][:]
        singles['sigmasq'] = trigs['sigmasq'][:]
        singles['snr'] = trigs['snr'][:]
        return numpy.array(singles, ndmin=1)

    def coinc(self, s0, s1, slide, step):
        # logsignalrate function inherited from PhaseTDStatistic
        logr_s = self.logsignalrate(s0, s1, slide * step)
        # rescale by ExpFitCombinedSNR reference slope as for sngl stat
        cstat = s0['snglstat'] + s1['snglstat'] + logr_s / self.alpharef
        # cut off underflowing and very small values
        cstat[cstat < 8.] = 8.
        # scale to resemble network SNR
        return cstat / (2.**0.5)

    def coinc_lim_for_thresh(self, s0, thresh):
        # if the threshold is below this value all triggers will
        # pass because of rounding in the coinc method
        if thresh <= (8. / (2.**0.5)):
            return -1. * numpy.ones(len(s0['snglstat'])) * numpy.inf
        if self.hist is None:
            self.get_hist()
        # Assume best case scenario and use maximum signal rate
        logr_s = self.hist_max
        s1 = (2. ** 0.5) * thresh - s0['snglstat'] - logr_s / self.alpharef
        return s1


class PhaseTDNewExpFitStatistic(PhaseTDNewStatistic, ExpFitCombinedSNR):
    """Statistic combining exponential noise model with signal histogram PDF"""

    # default is 2-ifo operation with exactly 1 'phasetd' file
    def __init__(self, files=None, ifos=None, **kwargs):
        # read in both foreground PDF and background fit info
        ExpFitCombinedSNR.__init__(self, files=files, ifos=ifos, **kwargs)
        # need the self.single_dtype value from PhaseTDStatistic
        PhaseTDNewStatistic.__init__(self, files=files, ifos=ifos, **kwargs)

    def single(self, trigs):
        # same single-ifo stat as ExpFitCombinedSNR
        sngl_stat = ExpFitCombinedSNR.single(self, trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase'][:]
        singles['end_time'] = trigs['end_time'][:]
        singles['sigmasq'] = trigs['sigmasq'][:]
        singles['snr'] = trigs['snr'][:]
        return numpy.array(singles, ndmin=1)

    def coinc(self, s0, s1, slide, step):
        # logsignalrate function inherited from PhaseTDStatistic
        logr_s = self.logsignalrate(s0, s1, slide * step)
        # rescale by ExpFitCombinedSNR reference slope as for sngl stat
        cstat = s0['snglstat'] + s1['snglstat'] + logr_s / self.alpharef
        # cut off underflowing and very small values
        cstat[cstat < 8.] = 8.
        # scale to resemble network SNR
        return cstat / (2.**0.5)

    def coinc_lim_for_thresh(self, s0, thresh):
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


class PhaseTDExpFitSGStatistic(PhaseTDExpFitStatistic):
    """Statistic combining exponential noise model with signal histogram PDF

    adding the sine-Gaussian veto to the single detector ranking
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        PhaseTDExpFitStatistic.__init__(self, files=files, ifos=ifos, **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto


class PhaseTDNewExpFitSGStatistic(PhaseTDNewExpFitStatistic):
    """Statistic combining exponential noise model with signal histogram PDF

    adding the sine-Gaussian veto to the single detector ranking
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        PhaseTDNewExpFitStatistic.__init__(self, files=files, ifos=ifos,
                                           **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto


class PhaseTDExpFitSGPSDStatistic(PhaseTDExpFitSGStatistic):
    """Statistic combining exponential noise model with signal histogram PDF

    adding the sine-Gaussian veto and PSD variation statistic to the
    single detector ranking
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        PhaseTDExpFitSGStatistic.__init__(self, files=files, ifos=ifos,
                                          **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar


class PhaseTDExpFitSGPSDScaledStatistic(PhaseTDExpFitSGStatistic):
    """Statistic combining exponential noise model with signal histogram PDF

    adding the sine-Gaussian veto and PSD variation statistic to the
    single detector ranking
    """

    def __init__(self, files=None, ifos=None, **kwargs):
        PhaseTDExpFitSGStatistic.__init__(self, files=files, ifos=ifos,
                                          **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar_scaled


class MaxContTradNewSNRStatistic(NewSNRStatistic):
    """Combination of NewSNR with the power chisq and auto chisq"""

    def single(self, trigs):
        """Calculate the single detector statistic.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
            Dictionary-like object holding single detector trigger information.
            'snr', 'cont_chisq', 'cont_chisq_dof', 'chisq_dof' and 'chisq'
            are required keys for this statistic.

        Returns
        -------
        stat: numpy.ndarray
            The array of single detector values
        """
        chisq_newsnr = ranking.get_newsnr(trigs)
        rautochisq = trigs['cont_chisq'][:] / trigs['cont_chisq_dof'][:]
        autochisq_newsnr = ranking.newsnr(trigs['snr'][:], rautochisq)
        return numpy.array(numpy.minimum(chisq_newsnr, autochisq_newsnr,
                           dtype=numpy.float32), ndmin=1, copy=False)


class ExpFitSGBgRateStatistic(ExpFitStatistic):
    """Detection statistic using an exponential falloff noise model.

    Statistic calculates the log noise coinc rate for each
    template over single-ifo newsnr values.
    """

    def __init__(self, files=None, ifos=None, benchmark_lograte=-14.6,
                 **kwargs):
        # benchmark_lograte is log of a representative noise trigger rate
        # This comes from H1L1 (O2) and is 4.5e-7 Hz
        super(ExpFitSGBgRateStatistic, self).__init__(files=files, ifos=ifos,
                                                      **kwargs)
        self.benchmark_lograte = benchmark_lograte
        self.get_newsnr = ranking.get_newsnr_sgveto

        # Reassign the rate to be number per time rather than an arbitrarily
        # normalised number
        for ifo in self.bg_ifos:
            self.reassign_rate(ifo)

    def reassign_rate(self, ifo):
        coeff_file = self.files[ifo+'-fit_coeffs']
        template_id = coeff_file['template_id'][:]
        # create arrays in template_id order for easier recall
        tid_sort = numpy.argsort(template_id)
        self.fits_by_tid[ifo]['rate'] = \
            coeff_file['count_above_thresh'][:][tid_sort] / \
            float(coeff_file.attrs['analysis_time'])

    def coinc_multiifo(self, s, slide, step, to_shift,
                       **kwargs): # pylint:disable=unused-argument
        # ranking statistic is -ln(expected rate density of noise triggers)
        # plus normalization constant
        sngl_dict = {sngl[0]: sngl[1] for sngl in s}
        ln_noise_rate = coinc_rate.combination_noise_lograte(
                                  sngl_dict, kwargs['time_addition'])
        loglr = - ln_noise_rate + self.benchmark_lograte
        return loglr

    def coinc_multiifo_lim_for_thresh(self, s, thresh, limifo, **kwargs):
        sngl_dict = {sngl[0]: sngl[1] for sngl in s}
        sngl_dict[limifo] = numpy.zeros(len(s[0][1]))
        ln_noise_rate = coinc_rate.combination_noise_lograte(
                                  sngl_dict, kwargs['time_addition'])
        loglr = - thresh - ln_noise_rate + self.benchmark_lograte
        return loglr


class ExpFitSGFgBgRateStatistic(PhaseTDStatistic, ExpFitSGBgRateStatistic):

    def __init__(self, files=None, ifos=None, **kwargs):
        # read in background fit info and store it
        ExpFitSGBgRateStatistic.__init__(self, files=files, ifos=ifos,
                                         **kwargs)
        # if ifos not already set, determine via background fit info
        self.ifos = self.ifos or self.bg_ifos
        # PhaseTD statistic single_dtype plus network sensitivity benchmark
        PhaseTDStatistic.__init__(self, files=files, ifos=self.ifos, **kwargs)
        self.single_dtype.append(('benchmark_logvol', numpy.float32))

        self.get_newsnr = ranking.get_newsnr_sgveto

        for ifo in self.bg_ifos:
            self.assign_median_sigma(ifo)
        # benchmark_logvol is a benchmark sensitivity array over template id
        hl_net_med_sigma = numpy.amin([self.fits_by_tid[ifo]['median_sigma']
                                       for ifo in ['H1', 'L1']], axis=0)
        self.benchmark_logvol = 3.0 * numpy.log(hl_net_med_sigma)
        self.single_increasing = False

    def assign_median_sigma(self, ifo):
        coeff_file = self.files[ifo + '-fit_coeffs']
        template_id = coeff_file['template_id'][:]
        tid_sort = numpy.argsort(template_id)
        self.fits_by_tid[ifo]['median_sigma'] = \
            coeff_file['median_sigma'][:][tid_sort]

    def single(self, trigs):
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
        # store benchmark log volume as single-ifo information since the coinc
        # method does not have access to template id
        singles['benchmark_logvol'] = self.benchmark_logvol[tnum]
        return numpy.array(singles, ndmin=1)

    def coinc_multiifo(self, s, slide, step, to_shift,
                       **kwargs): # pylint:disable=unused-argument
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
        # Get benchmark log volume as single-ifo information
        # NB benchmark logvol for a given template is not ifo-dependent
        # - choose the first ifo for convenience
        benchmark_logvol = s[0][1]['benchmark_logvol']
        network_logvol -= benchmark_logvol

        coincifos = [sngl[0] for sngl in s]
        # logsignalrate function from PhaseTDStatistic
        if ('H1' in coincifos and 'L1' in coincifos):
            # apply HL hist for HL & HLV coincs, keep only H/L info
            s_hl = [sngl[1] for sngl in s if sngl[0] in ['H1', 'L1']]
            shift_hl = [sh for sngl, sh in zip(s, to_shift) if \
                        sngl[0] in ['H1', 'L1']]
            logr_s = self.logsignalrate_multiifo(s_hl, slide * step, shift_hl)
        else:
            logr_s = self.logsignalrate_multiifo([sngl[1] for sngl in s],
                                                 slide * step, to_shift)

        loglr = logr_s + network_logvol - ln_noise_rate
        # cut off underflowing and very small values
        loglr[loglr < -30.] = -30.
        return loglr

    def coinc_multiifo_lim_for_thresh(self, s, thresh, limifo,
                                      **kwargs): # pylint:disable=unused-argument
        if self.hist is None:
            self.get_hist()
        # if the threshold is below this value all triggers will
        # pass because of rounding in the coinc method
        if thresh <= -30.:
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
        # Get benchmark log volume as single-ifo information
        # NB benchmark logvol for a given template is not ifo-dependent
        # - choose the first ifo for convenience
        benchmark_logvol = s[0][1]['benchmark_logvol']
        network_logvol -= benchmark_logvol

        loglr = - thresh + self.hist_max + network_logvol - ln_noise_rate
        return loglr


class ExpFitSGFgBgNormNewStatistic(PhaseTDNewStatistic,
                                   ExpFitSGBgRateStatistic):

    def __init__(self, files=None, ifos=None, **kwargs):
        # read in background fit info and store it
        ExpFitSGBgRateStatistic.__init__(self, files=files, ifos=ifos,
                                         **kwargs)
        # if ifos not already set, determine via background fit info
        self.ifos = self.ifos or self.bg_ifos
        # PhaseTD statistic single_dtype plus network sensitivity benchmark
        PhaseTDNewStatistic.__init__(self, files=files, ifos=self.ifos,
                                     **kwargs)
        self.single_dtype.append(('benchmark_logvol', numpy.float32))

        self.get_newsnr = ranking.get_newsnr_sgveto

        for ifo in self.bg_ifos:
            self.assign_median_sigma(ifo)
        # benchmark_logvol is a benchmark sensitivity array over template id
        hl_net_med_sigma = numpy.amin([self.fits_by_tid[ifo]['median_sigma']
                                       for ifo in ['H1', 'L1']], axis=0)
        self.benchmark_logvol = 3.0 * numpy.log(hl_net_med_sigma)
        self.single_increasing = False

    def assign_median_sigma(self, ifo):
        coeff_file = self.files[ifo + '-fit_coeffs']
        template_id = coeff_file['template_id'][:]
        tid_sort = numpy.argsort(template_id)
        self.fits_by_tid[ifo]['median_sigma'] = \
            coeff_file['median_sigma'][:][tid_sort]

    def lognoiserate(self, trigs, alphabelow=6):
        """Calculate the log noise rate density over single-ifo newsnr

        Read in single trigger information, make the newsnr statistic
        and rescale by the fitted coefficients alpha and rate
        """
        alphai, ratei, thresh = self.find_fits(trigs)
        newsnr = self.get_newsnr(trigs)
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

    def single_multiifo(self, s):
        ln_noise_rate = s[1]['snglstat']
        ln_noise_rate -= self.benchmark_lograte
        network_sigmasq = s[1]['sigmasq']
        network_logvol = 1.5 * numpy.log(network_sigmasq)
        benchmark_logvol = s[1]['benchmark_logvol']
        network_logvol -= benchmark_logvol
        ln_s = -4 * numpy.log(s[1]['snr'] / self.ref_snr)
        loglr = network_logvol - ln_noise_rate + ln_s
        # cut off underflowing and very small values
        loglr[loglr < -30.] = -30.
        return loglr

    def coinc_multiifo(self, s, slide, step, to_shift,
                       **kwargs): # pylint:disable=unused-argument
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
        logr_s = self.logsignalrate_multiifo(stat,
                                             slide * step, to_shift)

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

    def coinc_multiifo_lim_for_thresh(self, s, thresh, limifo,
                                      **kwargs): # pylint:disable=unused-argument
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


class ExpFitSGPSDFgBgNormStatistic(ExpFitSGFgBgNormNewStatistic):
    def __init__(self, files=None, ifos=None, **kwargs):
        ExpFitSGFgBgNormNewStatistic.__init__(self, files=files, ifos=ifos,
                                              **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar


class ExpFitSGPSDScaledFgBgNormStatistic(ExpFitSGFgBgNormNewStatistic):
    def __init__(self, files=None, ifos=None, **kwargs):
        ExpFitSGFgBgNormNewStatistic.__init__(self, files=files, ifos=ifos,
                                              **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar_scaled


class ExpFitSGPSDFgBgNormThreshStatistic(ExpFitSGFgBgNormNewStatistic):
    def __init__(self, files=None, ifos=None, **kwargs):
        ExpFitSGFgBgNormNewStatistic.__init__(self, files=files, ifos=ifos,
                                              **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar_threshold


class ExpFitSGPSDFgBgNormBBHStatistic(ExpFitSGFgBgNormNewStatistic):
    def __init__(self, files=None, ifos=None, max_chirp_mass=None, **kwargs):
        ExpFitSGFgBgNormNewStatistic.__init__(self, files=files, ifos=ifos,
                                              **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar
        self.mcm = max_chirp_mass
        self.curr_mchirp = None

    def single(self, trigs):
        from pycbc.conversions import mchirp_from_mass1_mass2
        self.curr_mchirp = mchirp_from_mass1_mass2(trigs.param['mass1'],
                                                   trigs.param['mass2'])
        if self.mcm is not None:
            # Careful - input might be a str, so cast to float
            self.curr_mchirp = min(self.curr_mchirp, float(self.mcm))
        return ExpFitSGFgBgNormNewStatistic.single(self, trigs)

    def logsignalrate_multiifo(self, stats, shift, to_shift):
        # model signal rate as uniform over chirp mass, background rate is
        # proportional to mchirp^(-11/3) due to density of templates
        logr_s = ExpFitSGFgBgNormNewStatistic.logsignalrate_multiifo(
                                                  self, stats, shift, to_shift)
        logr_s += numpy.log((self.curr_mchirp / 20.0) ** (11./3.0))
        return logr_s

    def coinc_multiifo_lim_for_thresh(self, s, thresh, limifo,
                                      **kwargs): # pylint:disable=unused-argument
        loglr = ExpFitSGFgBgNormNewStatistic.coinc_multiifo_lim_for_thresh(
                    self, s, thresh, limifo, **kwargs)
        loglr += numpy.log((self.curr_mchirp / 20.0) ** (11./3.0))
        return loglr


class ExpFitSGPSDFgBgNormBBHThreshStatistic(ExpFitSGPSDFgBgNormBBHStatistic):
    def __init__(self, files=None, ifos=None, max_chirp_mass=None, **kwargs):
        ExpFitSGPSDFgBgNormBBHStatistic.__init__(self, files=files, ifos=ifos,
                                                 max_chirp_mass=None, **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar_threshold


class ExpFitSGPSDSTFgBgNormBBHStatistic(ExpFitSGPSDFgBgNormBBHStatistic):
    def __init__(self, files=None, ifos=None, max_chirp_mass=None, **kwargs):
        ExpFitSGPSDFgBgNormBBHStatistic.__init__(self, files=files, ifos=ifos,
                                                 max_chirp_mass=None, **kwargs)
        self.get_newsnr = ranking.get_newsnr_sgveto_psdvar_scaled_threshold


statistic_dict = {
    'newsnr': NewSNRStatistic,
    'network_snr': NetworkSNRStatistic,
    'newsnr_cut': NewSNRCutStatistic,
    'phasetd_newsnr': PhaseTDStatistic,
    'phasetd_newsnr_sgveto': PhaseTDSGStatistic,
    'exp_fit_stat': ExpFitStatistic,
    'exp_fit_csnr': ExpFitCombinedSNR,
    'exp_fit_sg_csnr': ExpFitSGCombinedSNR,
    'exp_fit_sg_csnr_psdvar': ExpFitSGPSDCombinedSNR,
    'phasetd_exp_fit_stat': PhaseTDExpFitStatistic,
    'max_cont_trad_newsnr': MaxContTradNewSNRStatistic,
    'phasetd_exp_fit_stat_sgveto': PhaseTDExpFitSGStatistic,
    'phasetd_new_exp_fit_stat_sgveto': PhaseTDNewExpFitSGStatistic,
    'newsnr_sgveto': NewSNRSGStatistic,
    'newsnr_sgveto_psdvar': NewSNRSGPSDStatistic,
    'phasetd_exp_fit_stat_sgveto_psdvar': PhaseTDExpFitSGPSDStatistic,
    'phasetd_exp_fit_stat_sgveto_psdvar_scaled':
        PhaseTDExpFitSGPSDScaledStatistic,
    'exp_fit_sg_bg_rate': ExpFitSGBgRateStatistic,
    'exp_fit_sg_fgbg_rate': ExpFitSGFgBgRateStatistic,
    'exp_fit_sg_fgbg_norm_new': ExpFitSGFgBgNormNewStatistic,
    '2ogc': ExpFitSGPSDScaledFgBgNormStatistic, # backwards compatible
    '2ogcbbh': ExpFitSGPSDSTFgBgNormBBHStatistic, # backwards compatible
    'exp_fit_sg_fgbg_norm_psdvar': ExpFitSGPSDFgBgNormStatistic,
    'exp_fit_sg_fgbg_norm_psdvar_thresh': ExpFitSGPSDFgBgNormThreshStatistic,
    'exp_fit_sg_fgbg_norm_psdvar_bbh': ExpFitSGPSDFgBgNormBBHStatistic,
    'exp_fit_sg_fgbg_norm_psdvar_bbh_thresh':
        ExpFitSGPSDFgBgNormBBHThreshStatistic
}

sngl_statistic_dict = {
    'newsnr': NewSNRStatistic,
    'new_snr': NewSNRStatistic, # For backwards compatibility
    'snr': NetworkSNRStatistic,
    'newsnr_cut': NewSNRCutStatistic,
    'exp_fit_csnr': ExpFitCombinedSNR,
    'exp_fit_sg_csnr': ExpFitSGCombinedSNR,
    'max_cont_trad_newsnr': MaxContTradNewSNRStatistic,
    'newsnr_sgveto': NewSNRSGStatistic,
    'newsnr_sgveto_psdvar': NewSNRSGPSDStatistic,
    'newsnr_sgveto_psdvar_threshold': NewSNRSGPSDThresholdStatistic,
    'newsnr_sgveto_psdvar_scaled': NewSNRSGPSDScaledStatistic,
    'newsnr_sgveto_psdvar_scaled_threshold':
        NewSNRSGPSDScaledThresholdStatistic,
    'exp_fit_sg_csnr_psdvar': ExpFitSGPSDCombinedSNR
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


def get_sngl_statistic(stat):
    """
    Error-handling sugar around dict lookup for single-detector statistics

    Parameters
    ----------
    stat : string
        Name of the single-detector statistic

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
        return sngl_statistic_dict[stat]
    except KeyError:
        raise RuntimeError('%s is not an available detection statistic' % stat)
