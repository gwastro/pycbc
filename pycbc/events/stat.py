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
from hashlib import sha1
from datetime import datetime as dt
import numpy
import h5py

from . import ranking
from . import coinc_rate
from .eventmgr_cython import logsignalrateinternals_computepsignalbins
from .eventmgr_cython import logsignalrateinternals_compute2detrate

logger = logging.getLogger("pycbc.events.stat")

_allowed_statistic_features = [
    "phasetd",
    "kde",
    "dq",
    "chirp_mass",
    "sensitive_volume",
    "normalize_fit_rate",
]


class Stat(object):
    """Base class which should be extended to provide a statistic"""

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
            with h5py.File(filename, "r") as f:
                stat = f.attrs["stat"]
            if hasattr(stat, "decode"):
                stat = stat.decode()
            if stat in self.files:
                raise RuntimeError(
                    "We already have one file with stat attr = "
                    "%s. Can't provide more than one!" % stat
                )
            logger.info("Found file %s for stat %s", filename, stat)
            self.files[stat] = filename
        # Keep track of when stat files hashes so it can be
        # reloaded if it has changed
        self.file_hashes = self.get_file_hashes()

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
        self.kwargs = {}
        for key, value in kwargs.items():
            if key.startswith("sngl_ranking_"):
                # Note that all the sngl_ranking keywords are floats,
                # so this conversion is safe - if that changes in the future,
                # we may need something more sophisticated
                self.sngl_ranking_kwargs[key[13:]] = float(value)
            else:
                self.kwargs[key] = value

    def get_file_hashes(self):
        """
        Get sha1 hashes for all the files
        """
        logger.debug(
            "Getting file hashes"
        )
        start = dt.now()
        file_hashes = {}
        for stat, filename in self.files.items():
            with open(filename, 'rb') as file_binary:
                file_hashes[stat] = sha1(file_binary.read()).hexdigest()
        logger.debug(
            "Got file hashes for %d files, took %.3es",
            len(self.files),
            (dt.now() - start).total_seconds()
        )
        return file_hashes

    def files_changed(self):
        """
        Compare hashes of files now with the ones we have cached
        """
        changed_file_hashes = self.get_file_hashes()
        for stat, old_hash in self.file_hashes.items():
            if changed_file_hashes[stat] != old_hash:
                logger.info(
                    "%s statistic file %s has changed",
                    ''.join(self.ifos),
                    stat,
                )
            else:
                # Remove the dataset from the dictionary of hashes
                del changed_file_hashes[stat]

        if changed_file_hashes == {}:
            logger.debug(
                "No %s statistic files have changed",
                ''.join(self.ifos)
            )

        return list(changed_file_hashes.keys())

    def check_update_files(self):
        """
        Check whether files associated with the statistic need updated,
        then do so for each file which needs changing
        """
        files_changed = self.files_changed()
        for file_key in files_changed:
            self.update_file(file_key)
        self.file_hashes = self.get_file_hashes()

    def update_file(self, key):
        """
        Update file used in this statistic referenced by key.
        """
        err_msg = "This function is a stub that should be overridden by the "
        err_msg += "sub-classes. You shouldn't be seeing this error!"
        raise NotImplementedError(err_msg)

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

    def single(self, trigs):  # pylint:disable=unused-argument
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

    def rank_stat_single(self, single_info):  # pylint:disable=unused-argument
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

    def rank_stat_coinc(
        self, s, slide, step, to_shift, **kwargs
    ):  # pylint:disable=unused-argument
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

    def coinc_lim_for_thresh(
        self, s, thresh, limifo, **kwargs
    ):  # pylint:disable=unused-argument
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

        For this statistic this is just passing through the
        single value, which will be the second entry in the tuple.

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
        return single_info[1]

    def rank_stat_coinc(
        self, sngls_list, slide, step, to_shift, **kwargs
    ):  # pylint:disable=unused-argument
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

    def coinc_lim_for_thresh(
        self, s, thresh, limifo, **kwargs
    ):  # pylint:disable=unused-argument
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
        allowed_names = ["QuadratureSumStatistic"]
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

    def __init__(
        self,
        sngl_ranking,
        files=None,
        ifos=None,
        pregenerate_hist=True,
        **kwargs,
    ):
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
        QuadratureSumStatistic.__init__(
            self, sngl_ranking, files=files, ifos=ifos, **kwargs
        )

        self.single_dtype = [
            ("snglstat", numpy.float32),
            ("coa_phase", numpy.float32),
            ("end_time", numpy.float64),
            ("sigmasq", numpy.float32),
            ("snr", numpy.float32),
        ]

        # Assign attribute so that it can be replaced with other functions
        self.has_hist = False
        self.hist_ifos = None
        self.ref_snr = 5.
        self.relsense = {}
        self.swidth = self.pwidth = self.twidth = None
        self.srbmin = self.srbmax = None
        self.max_penalty = None
        self.pdtype = []
        self.weights = {}
        self.param_bin = {}
        self.two_det_flag = len(ifos) == 2
        self.two_det_weights = {}
        # Some memory
        self.pdif = numpy.zeros(128, dtype=numpy.float64)
        self.tdif = numpy.zeros(128, dtype=numpy.float64)
        self.sdif = numpy.zeros(128, dtype=numpy.float64)
        self.tbin = numpy.zeros(128, dtype=numpy.int32)
        self.pbin = numpy.zeros(128, dtype=numpy.int32)
        self.sbin = numpy.zeros(128, dtype=numpy.int32)

        # Is the histogram needed to be pre-generated?
        hist_needed = pregenerate_hist
        hist_needed &= not len(ifos) == 1
        hist_needed &= (type(self).__name__ == "PhaseTD" or self.kwargs["phasetd"])

        if hist_needed:
            self.get_hist()
        elif len(ifos) == 1:
            # remove all phasetd files from self.files and self.file_hashes,
            # as they are not needed
            for k in list(self.files.keys()):
                if 'phasetd_newsnr' in k:
                    del self.files[k]
                    del self.file_hashes[k]

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
            if "phasetd_newsnr" in name:
                ifokey = name.split("_")[2]
                num = len(ifokey) / 2
                if num != len(ifos):
                    continue

                match = [ifo in ifokey for ifo in ifos]
                if False in match:
                    continue
                selected = name
                break

        # If there are other phasetd_newsnr files, they aren't needed.
        # So tidy them out of the self.files dictionary
        rejected = [key for key in self.files.keys()
                    if 'phasetd_newsnr' in key and not key == selected]
        for k in rejected:
            del self.files[k]
            del self.file_hashes[k]

        if selected is None and len(ifos) > 1:
            raise RuntimeError("Couldn't figure out which stat file to use")
        if len(ifos) == 1:
            # We dont need the histogram file, but we are trying to get one
            # just skip it in this case
            return

        logger.info("Using signal histogram %s for ifos %s", selected, ifos)
        weights = {}
        param = {}

        with h5py.File(self.files[selected], "r") as histfile:
            self.hist_ifos = histfile.attrs["ifos"]

            # Patch for pre-hdf5=3.0 histogram files
            try:
                logger.info("Decoding hist ifos ..")
                self.hist_ifos = [i.decode("UTF-8") for i in self.hist_ifos]
            except (UnicodeDecodeError, AttributeError):
                pass

            # Histogram bin attributes
            self.twidth = histfile.attrs["twidth"]
            self.pwidth = histfile.attrs["pwidth"]
            self.swidth = histfile.attrs["swidth"]
            self.srbmin = histfile.attrs["srbmin"]
            self.srbmax = histfile.attrs["srbmax"]
            relfac = histfile.attrs["sensitivity_ratios"]

            for ifo in self.hist_ifos:
                weights[ifo] = histfile[ifo]["weights"][:]
                param[ifo] = histfile[ifo]["param_bin"][:]

        n_ifos = len(self.hist_ifos)

        bin_volume = (self.twidth * self.pwidth * self.swidth) ** (n_ifos - 1)
        self.hist_max = -1. * numpy.inf

        # Read histogram for each ifo, to use if that ifo has smallest SNR in
        # the coinc
        for ifo in self.hist_ifos:

            # renormalise to PDF
            self.weights[ifo] = \
                (weights[ifo] / (weights[ifo].sum() * bin_volume))
            self.weights[ifo] = self.weights[ifo].astype(numpy.float32)

            if param[ifo].dtype == numpy.int8:
                # Older style, incorrectly sorted histogram file
                ncol = param[ifo].shape[1]
                self.pdtype = [
                    ("c%s" % i, param[ifo].dtype) for i in range(ncol)
                ]
                self.param_bin[ifo] = numpy.zeros(
                    len(self.weights[ifo]), dtype=self.pdtype
                )
                for i in range(ncol):
                    self.param_bin[ifo]["c%s" % i] = param[ifo][:, i]

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
                if not hasattr(self, "c0_size"):
                    self.c0_size = {}
                    self.c1_size = {}
                    self.c2_size = {}

                self.c0_size[ifo] = numpy.int32(
                    2 * (abs(self.param_bin[ifo]["c0"]).max() + 1)
                )
                self.c1_size[ifo] = numpy.int32(
                    2 * (abs(self.param_bin[ifo]["c1"]).max() + 1)
                )
                self.c2_size[ifo] = numpy.int32(
                    2 * (abs(self.param_bin[ifo]["c2"]).max() + 1)
                )

                array_size = [
                    self.c0_size[ifo],
                    self.c1_size[ifo],
                    self.c2_size[ifo],
                ]
                dtypec = self.weights[ifo].dtype
                self.two_det_weights[ifo] = (
                    numpy.zeros(array_size, dtype=dtypec) + self.max_penalty
                )
                id0 = (
                    self.param_bin[ifo]["c0"].astype(numpy.int32)
                    + self.c0_size[ifo] // 2
                )
                id1 = (
                    self.param_bin[ifo]["c1"].astype(numpy.int32)
                    + self.c1_size[ifo] // 2
                )
                id2 = (
                    self.param_bin[ifo]["c2"].astype(numpy.int32)
                    + self.c2_size[ifo] // 2
                )
                self.two_det_weights[ifo][id0, id1, id2] = self.weights[ifo]

        for ifo, sense in zip(self.hist_ifos, relfac):
            self.relsense[ifo] = sense

        self.has_hist = True

    def update_file(self, key):
        """
        Update file used in this statistic.
        If others are used (i.e. this statistic is inherited), they will
        need updated separately
        """
        if 'phasetd_newsnr' in key and not len(self.ifos) == 1:
            if ''.join(sorted(self.ifos)) not in key:
                logger.debug(
                    "%s file is not used for %s statistic",
                    key,
                    ''.join(self.ifos)
                )
                return False
            logger.info(
                "Updating %s statistic %s file",
                ''.join(self.ifos),
                key
            )
            # This is a PhaseTDStatistic file which needs updating
            self.get_hist()
            return True
        return False

    def logsignalrate(self, stats, shift, to_shift):
        """
        Calculate the normalized log rate density of coinc signals via lookup

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
        snrs = numpy.array(
            [numpy.array(stats[ifo]["snr"], ndmin=1) for ifo in self.ifos]
        )
        smin = snrs.argmin(axis=0)
        # Store a list of the triggers using each ifo as reference
        rtypes = {
            ifo: numpy.where(smin == j)[0] for j, ifo in enumerate(self.ifos)
        }

        # Get reference ifo information
        rate = numpy.zeros(len(shift), dtype=numpy.float32)
        ps = {ifo: numpy.array(stats[ifo]['coa_phase'],
                               dtype=numpy.float32, ndmin=1)
              for ifo in self.ifos}
        ts = {ifo: numpy.array(stats[ifo]['end_time'],
                               dtype=numpy.float64, ndmin=1)
              for ifo in self.ifos}
        ss = {ifo: numpy.array(stats[ifo]['snr'],
                               dtype=numpy.float32, ndmin=1)
              for ifo in self.ifos}
        sigs = {ifo: numpy.array(stats[ifo]['sigmasq'],
                                 dtype=numpy.float32, ndmin=1)
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
                    length,
                )

                binned += [
                    self.tbin[:length],
                    self.pbin[:length],
                    self.sbin[:length],
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
                    len(rtype),
                )
            else:
                # Low[er]-RAM, high[er]-CPU option for >two det

                # Convert binned to same dtype as stored in hist
                nbinned = numpy.zeros(len(binned[1]), dtype=self.pdtype)
                for i, b in enumerate(binned):
                    nbinned[f"c{i}"] = b

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
                rate[rtype] *= (sref[rtype] / self.ref_snr) ** -4.

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
        singles["snglstat"] = sngl_stat
        singles["coa_phase"] = trigs["coa_phase"][:]
        singles["end_time"] = trigs["end_time"][:]
        singles["sigmasq"] = trigs["sigmasq"][:]
        singles["snr"] = trigs["snr"][:]
        return numpy.array(singles, ndmin=1)

    def rank_stat_single(self, single_info):
        """
        Calculate the statistic for a single detector candidate

        For this statistic this is just passing through the
        single value, which will be the second entry in the tuple.

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
        return single_info[1]["snglstat"]

    def rank_stat_coinc(
        self, sngls_list, slide, step, to_shift, **kwargs
    ):  # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic, defined in Eq 2 of
        [Nitz et al, 2017](https://doi.org/10.3847/1538-4357/aa8f50).
        """
        rstat = sum(s[1]["snglstat"] ** 2 for s in sngls_list)
        cstat = rstat + 2. * self.logsignalrate(
            dict(sngls_list), slide * step, to_shift
        )
        cstat[cstat < 0] = 0
        return cstat**0.5

    def coinc_lim_for_thresh(
        self, sngls_list, thresh, limifo, **kwargs
    ):  # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest.
        Calculate the required single detector statistic to exceed the
        threshold for each of the input triggers.
        """
        # Safety against subclassing and not rethinking this
        allowed_names = ["PhaseTDStatistic"]
        self._check_coinc_lim_subclass(allowed_names)

        if not self.has_hist:
            self.get_hist()

        fixed_stat_sq = sum(
            [b["snglstat"] ** 2 for a, b in sngls_list if a != limifo]
        )
        s1 = thresh ** 2. - fixed_stat_sq
        # Assume best case scenario and use maximum signal rate
        s1 -= 2. * self.hist_max
        s1[s1 < 0] = 0
        return s1**0.5


class ExpFitStatistic(PhaseTDStatistic):
    """
    Detection statistic using an exponential falloff noise model.

    Statistic approximates the negative log noise coinc rate density per
    template over single-ifo newsnr values.

    Extra features can be added by supplying keyword arguments when
    initialising.
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
        kwargs: values and features needed for the statistic
        """
        if not files:
            raise RuntimeError("Files not specified")

        PhaseTDStatistic.__init__(
            self, sngl_ranking, files=files, ifos=ifos, **kwargs
        )

        # Get the single-detector rates fit files
        # the stat file attributes are hard-coded as '%{ifo}-fit_coeffs'
        parsed_attrs = [f.split("-") for f in self.files.keys()]
        self.bg_ifos = [
            at[0]
            for at in parsed_attrs
            if (len(at) == 2 and at[1] == "fit_coeffs")
        ]
        if not len(self.bg_ifos):
            raise RuntimeError(
                "None of the statistic files has the required "
                "attribute called {ifo}-fit_coeffs !"
            )

        # Get the single-detector rates fit values
        self.fits_by_tid = {}
        for i in self.bg_ifos:
            self.fits_by_tid[i] = self.assign_fits(i)
            if self.kwargs["normalize_fit_rate"]:
                self.reassign_rate(i)

        # These are important for the coinc_lim_for_thresh method
        # Is the single threshold a < or > limit?
        self.single_increasing = False

        # Things for calculating the best-case scenario for signal rate
        self.max_sigmasq = -numpy.inf
        self.min_snr = numpy.inf

        # Some modifiers for the statistic to get it into a nice range
        self.benchmark_lograte = float(
            self.kwargs.get("benchmark_lograte", -14.6)
        )
        self.min_stat = float(
            self.kwargs.get("minimum_statistic_cutoff", -30.)
        )

        # Modifier to get a sensible value of the fit slope below threshold
        self.alphabelow = float(
            self.kwargs.get("alpha_below_thresh", numpy.inf)
        )

        # This will be used to keep track of the template number being used
        self.curr_tnum = None

        # Applies a constant offset to all statistic values in a given instance.
        # This can be used to e.g. change relative rankings between different
        # event types. Default is zero offset.
        self.stat_correction = float(
            self.kwargs.get("statistic_correction", 0)
        )

        # Go through the keywords and add class information as needed:
        if self.kwargs["sensitive_volume"]:
            # Add network sensitivity benchmark
            self.single_dtype.append(("benchmark_logvol", numpy.float32))
            # benchmark_logvol is a benchmark sensitivity array
            # over template id
            ref_ifos = self.kwargs.get("reference_ifos", "H1,L1").split(",")
            hl_net_med_sigma = numpy.nanmin(
                [self.fits_by_tid[ifo]["median_sigma"] for ifo in ref_ifos],
                axis=0,
            )
            self.benchmark_logvol = 3. * numpy.log(hl_net_med_sigma)

        if self.kwargs["dq"]:
            # Reweight the noise rate by the dq reweighting factor
            self.dq_rates_by_state = {}
            self.dq_bin_by_tid = {}
            self.dq_state_segments = None
            self.low_latency = False
            self.single_dtype.append(('dq_state', int))

            for ifo in self.ifos:
                key = f"{ifo}-dq_stat_info"
                if key in self.files.keys():
                    self.dq_rates_by_state[ifo] = self.assign_dq_rates(key)
                    self.dq_bin_by_tid[ifo] = self.assign_template_bins(key)
                    self.check_low_latency(key)
                    if not self.low_latency:
                        if self.dq_state_segments is None:
                            self.dq_state_segments = {}
                        self.dq_state_segments[ifo] = self.setup_segments(key)

        if self.kwargs["chirp_mass"]:
            # Reweight the signal rate by the chirp mass of the template
            # This may be stored as a float, so cast just in case
            self.mcm = float(self.kwargs.get("max_chirp_mass", numpy.inf))
            self.curr_mchirp = None

        if self.kwargs["kde"]:
            # Reweight the signal rate by a weighting factor from the KDE of
            # a signal population normalised by expected relative rate of noise
            # triggers for a template
            self.kde_names = []
            self.find_kdes()
            self.kde_by_tid = {}
            for kname in self.kde_names:
                self.assign_kdes(kname)

    def assign_template_bins(self, key):
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
        ifo = key.split("-")[0]
        with h5py.File(self.files[key], "r") as dq_file:
            tids = []
            bin_nums = []
            bin_grp = dq_file[f"{ifo}/bins"]
            for bin_name in bin_grp.keys():
                bin_tids = bin_grp[f"{bin_name}/tids"][:]
                tids = list(tids) + list(bin_tids.astype(int))
                bin_nums = list(bin_nums) + list([bin_name] * len(bin_tids))

        bin_dict = dict(zip(tids, bin_nums))
        return bin_dict

    def assign_dq_rates(self, key):
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
        ifo = key.split("-")[0]
        with h5py.File(self.files[key], "r") as dq_file:
            bin_grp = dq_file[f"{ifo}/bins"]
            dq_dict = {}
            for bin_name in bin_grp.keys():
                dq_dict[bin_name] = bin_grp[f"{bin_name}/dq_rates"][:]

        return dq_dict

    def setup_segments(self, key):
        """
        Store segments from stat file
        """
        ifo = key.split("-")[0]
        with h5py.File(self.files[key], "r") as dq_file:
            ifo_grp = dq_file[ifo]
            dq_state_segs_dict = {}
            for k in ifo_grp["dq_segments"].keys():
                seg_dict = {}
                seg_dict["start"] = \
                    ifo_grp[f"dq_segments/{k}/segment_starts"][:]
                seg_dict["end"] = ifo_grp[f"dq_segments/{k}/segment_ends"][:]
                dq_state_segs_dict[k] = seg_dict

        return dq_state_segs_dict

    def find_dq_noise_rate(self, trigs, dq_state):
        """Get dq values for a specific ifo and dq states"""

        dq_val = numpy.ones(len(dq_state))

        if self.curr_ifo in self.dq_rates_by_state:
            for i, st in enumerate(dq_state):
                if isinstance(self.curr_tnum, numpy.ndarray):
                    bin_name = self.dq_bin_by_tid[self.curr_ifo][self.curr_tnum[i]]
                else:
                    bin_name = self.dq_bin_by_tid[self.curr_ifo][self.curr_tnum]
                dq_val[i] = self.dq_rates_by_state[self.curr_ifo][bin_name][st]
        return dq_val

    def find_dq_state_by_time(self, ifo, times):
        """Get the dq state for an ifo at times"""
        dq_state = numpy.zeros(len(times), dtype=numpy.uint8)
        if ifo in self.dq_state_segments:
            from pycbc.events.veto import indices_within_times

            for k in self.dq_state_segments[ifo]:
                starts = self.dq_state_segments[ifo][k]["start"]
                ends = self.dq_state_segments[ifo][k]["end"]
                inds = indices_within_times(times, starts, ends)
                # states are named in file as 'dq_state_N', need to extract N
                dq_state[inds] = int(k[9:])
        return dq_state

    def check_low_latency(self, key):
        """
        Check if the statistic file indicates low latency mode.
        Parameters
        ----------
        key: str
            Statistic file key string.
        Returns
        -------
        None
        """
        ifo = key.split('-')[0]
        with h5py.File(self.files[key], 'r') as dq_file:
            ifo_grp = dq_file[ifo]
            if 'dq_segments' not in ifo_grp.keys():
                # if segs are not in file, we must be in LL
                if self.dq_state_segments is not None:
                    raise ValueError(
                        'Either all dq stat files must have segments or none'
                    )
                self.low_latency = True
            elif self.low_latency:
                raise ValueError(
                    'Either all dq stat files must have segments or none'
                )

    def reassign_rate(self, ifo):
        """
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

    def assign_fits(self, ifo):
        """
        Extract fits from single-detector rate fit files

        Parameters
        -----------
        ifo: str
            The detector to get fits for.

        Returns
        -------
        rate_dict: dict
            A dictionary containing the fit information in the `alpha`, `rate`
            and `thresh` keys.
        """
        coeff_file = h5py.File(self.files[f"{ifo}-fit_coeffs"], "r")
        template_id = coeff_file["template_id"][:]
        # the template_ids and fit coeffs are stored in an arbitrary order
        # create new arrays in template_id order for easier recall
        tid_sort = numpy.argsort(template_id)

        fits_by_tid_dict = {}
        fits_by_tid_dict["smoothed_fit_coeff"] = coeff_file["fit_coeff"][:][
            tid_sort
        ]
        fits_by_tid_dict["smoothed_rate_above_thresh"] = coeff_file[
            "count_above_thresh"
        ][:][tid_sort].astype(float)
        fits_by_tid_dict["smoothed_rate_in_template"] = coeff_file[
            "count_in_template"
        ][:][tid_sort].astype(float)
        if self.kwargs["sensitive_volume"]:
            fits_by_tid_dict["median_sigma"] = coeff_file["median_sigma"][:][
                tid_sort
            ].astype(float)

        # The by-template fits may have been stored in the smoothed fits file
        if "fit_by_template" in coeff_file:
            coeff_fbt = coeff_file["fit_by_template"]
            fits_by_tid_dict["fit_by_fit_coeff"] = coeff_fbt["fit_coeff"][:][
                tid_sort
            ]
            fits_by_tid_dict["fit_by_rate_above_thresh"] = coeff_fbt[
                "count_above_thresh"
            ][:][tid_sort].astype(float)
            fits_by_tid_dict["fit_by_rate_in_template"] = coeff_file[
                "count_in_template"
            ][:][tid_sort].astype(float)


        # Keep the fit threshold in fits_by_tid
        fits_by_tid_dict["thresh"] = coeff_file.attrs["stat_threshold"]

        coeff_file.close()

        return fits_by_tid_dict

    def update_file(self, key):
        """
        Update file used in this statistic.
        If others are used (i.e. this statistic is inherited), they will
        need updated separately
        """
        # First - check if the phasetd file has been updated, this is
        # inherited from the PhaseTDStatistic
        if PhaseTDStatistic.update_file(self, key):
            return True

        if key.endswith('-fit_coeffs'):
            # This is a ExpFitStatistic file which needs updating
            # Which ifo is it?
            ifo = key[:2]
            self.fits_by_tid[ifo] = self.assign_fits(ifo)
            if self.kwargs["normalize_fit_rate"]:
                self.reassign_rate(ifo)
            self.get_ref_vals(ifo)
            logger.info(
                "Updating %s statistic %s file",
                ''.join(self.ifos),
                key
            )
            return True

        # Is the key a KDE statistic file that we update here?
        if key.endswith('kde_file'):
            logger.info(
                "Updating %s statistic %s file",
                ''.join(self.ifos),
                key
            )
            kde_style = key.split('-')[0]
            self.assign_kdes(kde_style)
            return True

        # We also need to check if the DQ files have updated
        if key.endswith('dq_stat_info'):
            ifo = key.split('-')[0]
            logger.info(
                "Updating %s statistic %s file",
                ifo,
                key
            )
            self.dq_rates_by_state[ifo] = self.assign_dq_rates(key)
            self.dq_bin_by_tid[ifo] = self.assign_template_bins(key)
            return True

        return False

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
            ifo = trigs.ifo
        except AttributeError:
            ifo = trigs.get('ifo', None)
            if ifo is None:
                ifo = self.ifos[0]
            assert ifo in self.ifos

        # fits_by_tid is a dictionary of dictionaries of arrays
        # indexed by ifo / coefficient name / template_id
        alphai = self.fits_by_tid[ifo]["smoothed_fit_coeff"][self.curr_tnum]
        ratei = self.fits_by_tid[ifo]["smoothed_rate_above_thresh"][self.curr_tnum]
        thresh = self.fits_by_tid[ifo]["thresh"]

        return alphai, ratei, thresh

    def find_kdes(self):
        """
        Find which associated files are for the KDE reweighting
        """
        # The stat file attributes are hard-coded as 'signal-kde_file'
        # and 'template-kde_file'
        parsed_attrs = [f.split("-") for f in self.files.keys()]
        self.kde_names = [
            at[0]
            for at in parsed_attrs
            if (len(at) == 2 and at[1] == "kde_file")
        ]
        assert sorted(self.kde_names) == ["signal", "template"], (
            "Two KDE stat files are required, they should have stat attr "
            "'signal-kde_file' and 'template-kde_file' respectively"
        )

    def assign_kdes(self, kname):
        """
        Extract values from KDE files

        Parameters
        -----------
        kname: str
            Used to label the kde files.
        """
        with h5py.File(self.files[kname + "-kde_file"], "r") as kde_file:
            self.kde_by_tid[kname + "_kdevals"] = kde_file["data_kde"][:]

    def kde_ratio(self):
        """
        Calculate the weighting factor according to the ratio of the
        signal and template KDE lookup tables
        """
        signal_kde = self.kde_by_tid["signal_kdevals"][self.curr_tnum]
        template_kde = self.kde_by_tid["template_kdevals"][self.curr_tnum]

        return numpy.log(signal_kde / template_kde)

    def lognoiserate(self, trigs):
        """
        Calculate the log noise rate density over single-ifo ranking

        Read in single trigger information, make the sngl_stat
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
        # What is the template number currently being used?
        try:
            # exists if accessed via coinc_findtrigs, this is a number
            self.curr_tnum = trigs.template_num
        except AttributeError:
            # exists for SingleDetTriggers & pycbc_live get_coinc,
            # this is a numpy array
            self.curr_tnum = trigs["template_id"]

        alphai, ratei, thresh = self.find_fits(trigs)
        sngl_stat = self.get_sngl_ranking(trigs)
        lognoisel = (
            -alphai * (sngl_stat - thresh)
            + numpy.log(alphai)
            + numpy.log(ratei)
        )

        if not numpy.isinf(self.alphabelow):
            # Above the threshold we use the usual fit coefficient (alphai)
            # below threshold use specified alphabelow
            bt = sngl_stat < thresh
            lognoiselbt = (
                -self.alphabelow * (sngl_stat - thresh)
                + numpy.log(self.alphabelow)
                + numpy.log(ratei)
            )
            lognoisel[bt] = lognoiselbt[bt]

        return numpy.array(lognoisel, ndmin=1, dtype=numpy.float32)

    def single(self, trigs):
        """
        Calculate the necessary single detector information

        In this case the ranking rescaled (see the lognoiserate method here).
        with the phase, end time, SNR values added in.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays, h5py group or similar dict-like object
            Object holding single detector trigger information.

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """

        try:
            self.curr_ifo = trigs.ifo
        except AttributeError:
            self.curr_ifo = trigs.get('ifo', None)
            if self.curr_ifo is None:
                self.curr_ifo = self.ifos[0]
            assert self.curr_ifo in self.ifos
            # Should be exactly one ifo provided
            assert len(numpy.unique(self.curr_ifo)) == 1

        # single-ifo stat = log of noise rate
        sngl_stat = self.lognoiserate(trigs)
        # populate other fields to calculate phase/time/amp consistency
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles["coa_phase"] = trigs["coa_phase"][:]
        singles["end_time"] = trigs["end_time"][:]
        singles["snr"] = trigs["snr"][:]
        # Save info about best-case scenario for use later
        if singles["snr"].size:
            self.min_snr = min(singles["snr"].min(), self.min_snr)

        if self.kwargs["sensitive_volume"]:
            # populate fields to allow sensitive volume factor calculation
            singles["sigmasq"] = trigs["sigmasq"][:]
            # Store benchmark log volume as single-ifo information since
            # the ranking methods do not have access to template id
            singles["benchmark_logvol"] = self.benchmark_logvol[self.curr_tnum]

            # Save info about the best-case scenario for use later:
            if singles["sigmasq"].size:
                max_sigsq = numpy.max(singles["sigmasq"])
                self.max_sigmasq = max(max_sigsq, self.max_sigmasq)

        if self.kwargs["chirp_mass"]:
            from pycbc.conversions import mchirp_from_mass1_mass2

            try:
                mass1 = trigs.param['mass1']
                mass2 = trigs.param['mass2']
            except AttributeError:
                mass1 = trigs['mass1']
                mass2 = trigs['mass2']
            self.curr_mchirp = mchirp_from_mass1_mass2(mass1, mass2)

        if self.kwargs["dq"]:
            if self.low_latency:
                # trigs should already have a dq state assigned
                singles['dq_state'] = trigs['dq_state'][:]
            else:
                singles['dq_state'] = self.find_dq_state_by_time(
                    self.curr_ifo, trigs['end_time'][:]
                )
            dq_rate = self.find_dq_noise_rate(trigs, singles['dq_state'])
            dq_rate = numpy.maximum(dq_rate, 1)
            sngl_stat += numpy.log(dq_rate)

        singles["snglstat"] = sngl_stat
        return numpy.array(singles, ndmin=1)

    def sensitive_volume_factor(self, sngls):
        """
        Determine the factor for different network sensitivities

        Assuming a homogeneous distribution of sources, the signal rate
        should be given by the volume of the sphere to which we are
        sensitive.

        Parameters
        ----------
        sngls: tuple
            Single-detector information, tuples of the (ifo, sngl_object),
            where sngl_object is a ReadByTemplate object or similar.

        Returns
        -------
        network_logvol: numpy.array
            Array of values for the network log-volume factor. This is the
            log of the cube of the sensitive distance (sigma), divided by
            a benchmark volume.
        """
        # Get benchmark log volume as single-ifo information :
        # benchmark_logvol for a given template is not ifo-dependent, so
        # choose the first ifo for convenience
        benchmark_logvol = sngls[0][1]["benchmark_logvol"]

        # Benchmark log volume will be the same for all triggers, so if
        # any are nan, they are all nan
        if any(numpy.isnan(benchmark_logvol)):
            # This can be the case in pycbc live if there are no triggers
            # from this template in the trigger fits file. If so, assume 
            # that sigma for the triggers being ranked is
            # representative of the benchmark network.
            return 0

        # Network sensitivity for a given coinc type is approximately
        # determined by the least sensitive ifo
        network_sigmasq = numpy.amin(
            [sngl[1]["sigmasq"] for sngl in sngls], axis=0
        )
        # Volume \propto sigma^3 or sigmasq^1.5
        network_logvol = 1.5 * numpy.log(network_sigmasq) - benchmark_logvol

        return network_logvol

    def logsignalrate_shared(self, sngls_info):
        """
        Calculate the parts of the log signal rate which are shared by
        both the single and coinc ranking statistics

        Parameters
        ----------
        sngls_info: tuple
            Single-detector information, tuples of the (ifo, sngl_object),
            where sngl_object is a ReadByTemplate object or similar.

        Returns
        -------
        sr_factor: numpy.ndarray
            Array of values to be added to the ranking statistic when
            taking various signal rate factors into account
        """
        # Other features affecting the signal rate
        sr_factor = 0
        if self.kwargs["sensitive_volume"]:
            # Sensitive volume - this is proportional to signal rate
            # assuming a homogeneous universe
            sr_factor += self.sensitive_volume_factor(sngls_info)

        if self.kwargs["chirp_mass"]:
            # chirp mass reweighting
            if isinstance(self.curr_mchirp, numpy.ndarray):
                mchirp = numpy.minimum(self.curr_mchirp, self.mcm)
            else:
                # curr_mchirp will be a number
                mchirp = min(self.curr_mchirp, self.mcm)
            sr_factor += numpy.log((mchirp / 20.) ** (11. / 3.))

        if self.kwargs["kde"]:
            # KDE reweighting
            sr_factor += self.kde_ratio()

        return sr_factor

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

        # Basic noise rate: the exp fit rate from the single statistic
        ln_noise_rate = sngls["snglstat"]
        ln_noise_rate -= self.benchmark_lograte

        # Basic signal rate: snr ** -4
        ln_s = -4 * numpy.log(sngls["snr"] / self.ref_snr)
        # Add in the feature-dependent terms to the signal rate
        ln_s += self.logsignalrate_shared((single_info,))

        # Combine the signal and noise rates
        loglr = ln_s - ln_noise_rate

        # Apply statistic correction
        loglr += self.stat_correction

        # cut off underflowing and very small values
        loglr[loglr < self.min_stat] = self.min_stat
        return loglr

    def rank_stat_coinc(
        self, s, slide, step, to_shift, **kwargs
    ):  # pylint:disable=unused-argument
        """
        Calculate the coincident detection statistic.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples
        slide: numpy array
            The number of steps by which to shift each trigger
        step: float
            The size of each step, to be multiplied by to_shift
        to_shift: list
            List of integers indicating what multiples of the time shift will
            be applied
        kwargs: various
            Key-word arguments to define different features and tunable
            values in the statistic. Only time_addition is used here

        Returns
        -------
        loglr: numpy.ndarray
            The ranking statistic for each trigger based on the supplied
            triggers, requested features and keyword arguments

        """
        # ranking statistic is -ln(expected rate density of noise triggers)
        # plus normalization constant
        sngl_dict = {sngl[0]: sngl[1]["snglstat"] for sngl in s}

        # Basic noise rate: sum of noise rates multiplied by the
        # window they can form coincidences in
        ln_noise_rate = coinc_rate.combination_noise_lograte(
            sngl_dict,
            kwargs["time_addition"],
            dets=kwargs.get('dets', None),
        )

        ln_noise_rate -= self.benchmark_lograte

        # Basic option is not to have any signal-based assumptions,
        # so this is zero to begin with
        ln_s = 0

        if self.kwargs["phasetd"]:
            if not self.has_hist:
                self.get_hist()
            # Find total volume of phase-time-amplitude space occupied by
            # noise coincs, so that the logsignalrate function is properly
            # normalized
            # Extent of time-difference space occupied
            noise_twindow = coinc_rate.multiifo_noise_coincident_area(
                self.hist_ifos,
                kwargs["time_addition"],
                dets=kwargs.get('dets', None),
            )
            # Volume is the allowed time difference window, multiplied by 2pi
            # for each phase difference dimension and by allowed range of SNR
            # ratio for each SNR ratio dimension : there are (n_ifos - 1)
            # dimensions for both phase and SNR
            n_ifos = len(self.hist_ifos)
            snr_range = (self.srbmax - self.srbmin) * self.swidth
            hist_vol = noise_twindow * (2. * numpy.pi * snr_range) ** (
                n_ifos - 1
            )
            # Noise PDF is 1/volume, assuming a uniform distribution of noise
            # coincs
            ln_noise_rate -= numpy.log(hist_vol)

            # What is the signal pdf?
            stat = {ifo: st for ifo, st in s}
            ln_s += self.logsignalrate(stat, slide * step, to_shift)

        # Add in the feature-dependent terms to the signal rate
        ln_s += self.logsignalrate_shared(s)

        # Combine the signal and noise rates
        loglr = ln_s - ln_noise_rate

        # Apply statistic correction
        loglr += self.stat_correction

        # cut off underflowing and very small values
        loglr[loglr < self.min_stat] = self.min_stat

        return loglr

    def coinc_lim_for_thresh(
        self, s, thresh, limifo, **kwargs
    ):  # pylint:disable=unused-argument
        """
        Optimization function to identify coincs too quiet to be of interest

        We are trying to get rid of as many events here at the point where
        we can be confident that they will not meet ranking statistic
        thresholds.

        The calculation here is "What is the minimum required snglstat in
        the pivot IFO which could possibly pass the threshold?"

        There are a couple of points to be wary of here, e.g. in the signal
        rate calculation, we take the best-case scenario. By using the
        best-case for signal rate in this calculation, some events are kept
        at this point which are hopeless.

        Parameters
        ----------
        s: list
            List of (ifo, single detector statistic) tuples
        thresh: float
            The threshold value to be checked against
        limifo: string
            The pivot detector, i.e. the one in which the sngl stat must
            reach the value output by ths function
        kwargs: various
            Key-word arguments to define different features and tunable
            values in the statistic. Only time_addition is used here

        Returns
        -------
        numpy.array
            The minimum snglstat values required in the 'pivot' detector
            in order to reach the threshold specified
        """
        # Safety against subclassing and not rethinking this
        allowed_names = ["ExpFitStatistic"]
        self._check_coinc_lim_subclass(allowed_names)

        if thresh <= self.min_stat:
            return numpy.ones(len(s[0][1]["snglstat"])) * numpy.inf

        # Modify the sngls so that the pivot ifo snglstats are zero
        sngl_dict = {sngl[0]: sngl[1]["snglstat"] for sngl in s}
        sngl_dict[limifo] = numpy.zeros(len(s[0][1]))

        # Noise rate calculated as normal. Because of the modification
        # above, this is the rank_stat_coinc noise rate calculation
        # minus the pivot ifo's snglstat
        ln_noise_rate = coinc_rate.combination_noise_lograte(
            sngl_dict, kwargs["time_addition"]
        )
        ln_noise_rate -= self.benchmark_lograte

        # Basic option is not to have any signal-based assumptions,
        # so this is zero to begin with
        ln_s = 0

        if self.kwargs["phasetd"]:
            if not self.has_hist:
                self.get_hist()
            # Assume best-case scenario and use maximum signal rate
            ln_s = numpy.log(
                self.hist_max * (self.min_snr / self.ref_snr) ** -4.
            )

        # Shared info is the same as in the coinc calculation
        ln_s += self.logsignalrate_shared(s)

        # Combine the signal and noise rates
        loglr = ln_s - ln_noise_rate

        # Apply statistic correction
        loglr += self.stat_correction

        # From this combined rate, what is the minimum snglstat value
        # in the pivot IFO needed to reach the threshold?
        return loglr - thresh


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
        ExpFitStatistic.__init__(
            self, sngl_ranking, files=files, ifos=ifos, **kwargs
        )
        # for low-mass templates the exponential slope alpha \approx 6
        self.alpharef = 6.
        self.single_increasing = True
        self.single_dtype = numpy.float32

        # Modifier to get a sensible value of the fit slope below threshold
        self.alphabelow = float(
            self.kwargs.get("alpha_below_thresh", numpy.inf)
        )

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
            sngl_multiifo = single_info[1]
        else:
            sngl_multiifo = -1. * single_info[1]
        return sngl_multiifo

    def rank_stat_coinc(
        self, s, slide, step, to_shift, **kwargs
    ):  # pylint:disable=unused-argument
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
        return sum(sngl[1] for sngl in s) / (len(s) ** 0.5)

    def coinc_lim_for_thresh(
        self, s, thresh, limifo, **kwargs
    ):  # pylint:disable=unused-argument
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
        allowed_names = ["ExpFitCombinedSNR"]
        self._check_coinc_lim_subclass(allowed_names)

        return thresh * ((len(s) + 1) ** 0.5) - sum(sngl[1] for sngl in s)


statistic_dict = {
    "quadsum": QuadratureSumStatistic,
    "single_ranking_only": QuadratureSumStatistic,
    "phasetd": PhaseTDStatistic,
    "exp_fit_csnr": ExpFitCombinedSNR,
    "exp_fit": ExpFitStatistic,
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
        raise RuntimeError("%s is not an available detection statistic" % stat)


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
        help="The coinc ranking statistic to calculate",
    )

    statistic_opt_group.add_argument(
        "--sngl-ranking",
        choices=ranking.sngls_ranking_function_dict.keys(),
        required=True,
        help="The single-detector trigger ranking to use.",
    )

    statistic_opt_group.add_argument(
        "--statistic-files",
        nargs="*",
        action="append",
        default=[],
        help="Files containing ranking statistic info",
    )

    statistic_opt_group.add_argument(
        "--statistic-keywords",
        nargs="*",
        default=[],
        help="Provide additional key-word arguments to be sent to "
        "the statistic class when it is initialized. Should "
        "be given in format --statistic-keywords "
        "KWARG1:VALUE1 KWARG2:VALUE2 KWARG3:VALUE3 ...",
    )

    statistic_opt_group.add_argument(
        "--statistic-features",
        nargs="*",
        default=[],
        choices=_allowed_statistic_features,
        help="Provide additional arguments to include particular "
        "features in the ranking statistic.",
    )

    return statistic_opt_group


def parse_statistic_feature_options(stat_features, stat_kwarg_list):
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

    # Check that the statistic keywords are allowed
    for feature in stat_features:
        if feature not in _allowed_statistic_features:
            err_msg = f"--statistic-feature {feature} not recognised"
            raise NotImplementedError(err_msg)

    # Set values for each feature key to a boolean of whether we want them
    for feature in _allowed_statistic_features:
        stat_kwarg_dict[feature] = feature in stat_features

    for inputstr in stat_kwarg_list:
        try:
            key, value = inputstr.split(":")
            stat_kwarg_dict[key] = value
        except ValueError:
            err_txt = (
                "--statistic-keywords must take input in the "
                "form KWARG1:VALUE1 KWARG2:VALUE2 KWARG3:VALUE3 ... "
                "Received {}".format(" ".join(stat_kwarg_list))
            )
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
    # if needed (e.g. not calling get_statistic_from_opts in a loop)
    if len(opts.statistic_files) > 0 and \
            isinstance(opts.statistic_files[0], list):
        opts.statistic_files = sum(opts.statistic_files, [])

    extra_kwargs = parse_statistic_feature_options(
        opts.statistic_features,
        opts.statistic_keywords,
    )

    stat_class = get_statistic(opts.ranking_statistic)(
        opts.sngl_ranking, opts.statistic_files, ifos=ifos, **extra_kwargs
    )

    return stat_class
