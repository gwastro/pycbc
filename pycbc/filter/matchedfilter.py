# Copyright (C) 2012  Alex Nitz
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
This modules provides functions for matched filtering along with associated
utilities.
"""

import logging
from math import sqrt
from pycbc.types import TimeSeries, FrequencySeries, zeros, Array
from pycbc.types import complex_same_precision_as, real_same_precision_as
from pycbc.fft import fft, ifft, IFFT
import pycbc.scheme
from pycbc import events
from pycbc.events import ranking
import pycbc
import numpy

BACKEND_PREFIX="pycbc.filter.matchedfilter_"

@pycbc.scheme.schemed(BACKEND_PREFIX)
def correlate(x, y, z):
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)


class BatchCorrelator(object):
    """ Create a batch correlation engine
    """
    def __init__(self, xs, zs, size):
        """ Correlate x and y, store in z. Arrays need not be equal length, but
        must be at least size long and of the same dtype. No error checking
        will be performed, so be careful. All dtypes must be complex64.
        Note, must be created within the processing context that it will be used in.
        """
        self.size = int(size)
        self.dtype = xs[0].dtype
        self.num_vectors = len(xs)

        # keep reference to arrays
        self.xs = xs
        self.zs = zs

        # Store each pointer as in integer array
        self.x = Array([v.ptr for v in xs], dtype=numpy.int)
        self.z = Array([v.ptr for v in zs], dtype=numpy.int)

    @pycbc.scheme.schemed(BACKEND_PREFIX)
    def batch_correlate_execute(self, y):
        pass

    execute = batch_correlate_execute


@pycbc.scheme.schemed(BACKEND_PREFIX)
def _correlate_factory(x, y, z):
    err_msg = "This class is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)


class Correlator(object):
    """ Create a correlator engine

    Parameters
    ---------
    x : complex64
      Input pycbc.types.Array (or subclass); it will be conjugated
    y : complex64
      Input pycbc.types.Array (or subclass); it will not be conjugated
    z : complex64
      Output pycbc.types.Array (or subclass).
      It will contain conj(x) * y, element by element

    The addresses in memory of the data of all three parameter vectors
    must be the same modulo pycbc.PYCBC_ALIGNMENT
    """
    def __new__(cls, *args, **kwargs):
        real_cls = _correlate_factory(*args, **kwargs)
        return real_cls(*args, **kwargs) # pylint:disable=not-callable


# The class below should serve as the parent for all schemed classes.
# The intention is that this class serves simply as the location for
# all documentation of the class and its methods, though that is not
# yet implemented.  Perhaps something along the lines of:
#
#    http://stackoverflow.com/questions/2025562/inherit-docstrings-in-python-class-inheritance
#
# will work? Is there a better way?
class _BaseCorrelator(object):
    def correlate(self):
        """
        Compute the correlation of the vectors specified at object
        instantiation, writing into the output vector given when the
        object was instantiated. The intention is that this method
        should be called many times, with the contents of those vectors
        changing between invocations, but not their locations in memory
        or length.
        """
        pass


class MatchedFilterControl(object):
    def __init__(self, low_frequency_cutoff, high_frequency_cutoff, snr_threshold, tlen,
                 delta_f, dtype, segment_list, template_output, use_cluster,
                 downsample_factor=1, upsample_threshold=1, upsample_method='pruned_fft',
                 gpu_callback_method='none', cluster_function='symmetric'):
        """ Create a matched filter engine.

        Parameters
        ----------
        low_frequency_cutoff : {None, float}, optional
            The frequency to begin the filter calculation. If None, begin at the
            first frequency after DC.
        high_frequency_cutoff : {None, float}, optional
            The frequency to stop the filter calculation. If None, continue to the
            the nyquist frequency.
        snr_threshold : float
            The minimum snr to return when filtering
        segment_list : list
            List of FrequencySeries that are the Fourier-transformed data segments
        template_output : complex64
            Array of memory given as the 'out' parameter to waveform.FilterBank
        use_cluster : boolean
            If true, cluster triggers above threshold using a window; otherwise,
            only apply a threshold.
        downsample_factor : {1, int}, optional
            The factor by which to reduce the sample rate when doing a heirarchical
            matched filter
        upsample_threshold : {1, float}, optional
            The fraction of the snr_threshold to trigger on the subsampled filter.
        upsample_method : {pruned_fft, str}
            The method to upsample or interpolate the reduced rate filter.
        cluster_function : {symmetric, str}, optional
            Which method is used to cluster triggers over time. If 'findchirp', a
            sliding forward window; if 'symmetric', each window's peak is compared
            to the windows before and after it, and only kept as a trigger if larger
            than both.
        """
        # Assuming analysis time is constant across templates and segments, also
        # delta_f is constant across segments.
        self.tlen = tlen
        self.flen = self.tlen / 2 + 1
        self.delta_f = delta_f
        self.delta_t = 1.0/(self.delta_f * self.tlen)
        self.dtype = dtype
        self.snr_threshold = snr_threshold
        self.flow = low_frequency_cutoff
        self.fhigh = high_frequency_cutoff
        self.gpu_callback_method = gpu_callback_method
        if cluster_function not in ['symmetric', 'findchirp']:
            raise ValueError("MatchedFilter: 'cluster_function' must be either 'symmetric' or 'findchirp'")
        self.cluster_function = cluster_function
        self.segments = segment_list
        self.htilde = template_output

        if downsample_factor == 1:
            self.snr_mem = zeros(self.tlen, dtype=self.dtype)
            self.corr_mem = zeros(self.tlen, dtype=self.dtype)

            if use_cluster and (cluster_function == 'symmetric'):
                self.matched_filter_and_cluster = self.full_matched_filter_and_cluster_symm
                # setup the threasholding/clustering operations for each segment
                self.threshold_and_clusterers = []
                for seg in self.segments:
                    thresh = events.ThresholdCluster(self.snr_mem[seg.analyze])
                    self.threshold_and_clusterers.append(thresh)
            elif use_cluster and (cluster_function == 'findchirp'):
                self.matched_filter_and_cluster = self.full_matched_filter_and_cluster_fc
            else:
                self.matched_filter_and_cluster = self.full_matched_filter_thresh_only

            # Assuming analysis time is constant across templates and segments, also
            # delta_f is constant across segments.
            self.kmin, self.kmax = get_cutoff_indices(self.flow, self.fhigh,
                                                      self.delta_f, self.tlen)

            # Set up the correlation operations for each analysis segment
            corr_slice = slice(self.kmin, self.kmax)
            self.correlators = []
            for seg in self.segments:
                corr = Correlator(self.htilde[corr_slice],
                                  seg[corr_slice],
                                  self.corr_mem[corr_slice])
                self.correlators.append(corr)

            # setup up the ifft we will do
            self.ifft = IFFT(self.corr_mem, self.snr_mem)

        elif downsample_factor >= 1:
            self.matched_filter_and_cluster = self.heirarchical_matched_filter_and_cluster
            self.downsample_factor = downsample_factor
            self.upsample_method = upsample_method
            self.upsample_threshold = upsample_threshold

            N_full = self.tlen
            N_red = N_full / downsample_factor
            self.kmin_full, self.kmax_full = get_cutoff_indices(self.flow,
                                              self.fhigh, self.delta_f, N_full)

            self.kmin_red, _ = get_cutoff_indices(self.flow,
                                                  self.fhigh, self.delta_f, N_red)

            if self.kmax_full < N_red:
                self.kmax_red = self.kmax_full
            else:
                self.kmax_red = N_red - 1

            self.snr_mem = zeros(N_red, dtype=self.dtype)
            self.corr_mem_full = FrequencySeries(zeros(N_full, dtype=self.dtype), delta_f=self.delta_f)
            self.corr_mem = Array(self.corr_mem_full[0:N_red], copy=False)
            self.inter_vec = zeros(N_full, dtype=self.dtype)

        else:
            raise ValueError("Invalid downsample factor")

    def full_matched_filter_and_cluster_symm(self, segnum, template_norm, window, epoch=None):
        """ Returns the complex snr timeseries, normalization of the complex snr,
        the correlation vector frequency series, the list of indices of the
        triggers, and the snr values at the trigger locations. Returns empty
        lists for these for points that are not above the threshold.

        Calculated the matched filter, threshold, and cluster.

        Parameters
        ----------
        segnum : int
            Index into the list of segments at MatchedFilterControl construction
            against which to filter.
        template_norm : float
            The htilde, template normalization factor.
        window : int
            Size of the window over which to cluster triggers, in samples

        Returns
        -------
        snr : TimeSeries
            A time series containing the complex snr.
        norm : float
            The normalization of the complex snr.
        corrrelation: FrequencySeries
            A frequency series containing the correlation vector.
        idx : Array
            List of indices of the triggers.
        snrv : Array
            The snr values at the trigger locations.
        """
        norm = (4.0 * self.delta_f) / sqrt(template_norm)
        self.correlators[segnum].correlate()
        self.ifft.execute()
        snrv, idx = self.threshold_and_clusterers[segnum].threshold_and_cluster(self.snr_threshold / norm, window)

        if len(idx) == 0:
            return [], [], [], [], []

        logging.info("%s points above threshold" % str(len(idx)))

        snr = TimeSeries(self.snr_mem, epoch=epoch, delta_t=self.delta_t, copy=False)
        corr = FrequencySeries(self.corr_mem, delta_f=self.delta_f, copy=False)
        return snr, norm, corr, idx, snrv

    def full_matched_filter_and_cluster_fc(self, segnum, template_norm, window, epoch=None):
        """ Returns the complex snr timeseries, normalization of the complex snr,
        the correlation vector frequency series, the list of indices of the
        triggers, and the snr values at the trigger locations. Returns empty
        lists for these for points that are not above the threshold.

        Calculated the matched filter, threshold, and cluster.

        Parameters
        ----------
        segnum : int
            Index into the list of segments at MatchedFilterControl construction
            against which to filter.
        template_norm : float
            The htilde, template normalization factor.
        window : int
            Size of the window over which to cluster triggers, in samples

        Returns
        -------
        snr : TimeSeries
            A time series containing the complex snr.
        norm : float
            The normalization of the complex snr.
        corrrelation: FrequencySeries
            A frequency series containing the correlation vector.
        idx : Array
            List of indices of the triggers.
        snrv : Array
            The snr values at the trigger locations.
        """
        norm = (4.0 * self.delta_f) / sqrt(template_norm)
        self.correlators[segnum].correlate()
        self.ifft.execute()
        idx, snrv = events.threshold(self.snr_mem[self.segments[segnum].analyze],
                                     self.snr_threshold / norm)
        idx, snrv = events.cluster_reduce(idx, snrv, window)

        if len(idx) == 0:
            return [], [], [], [], []

        logging.info("%s points above threshold" % str(len(idx)))

        snr = TimeSeries(self.snr_mem, epoch=epoch, delta_t=self.delta_t, copy=False)
        corr = FrequencySeries(self.corr_mem, delta_f=self.delta_f, copy=False)
        return snr, norm, corr, idx, snrv

    def full_matched_filter_thresh_only(self, segnum, template_norm, window=None, epoch=None):
        """ Returns the complex snr timeseries, normalization of the complex snr,
        the correlation vector frequency series, the list of indices of the
        triggers, and the snr values at the trigger locations. Returns empty
        lists for these for points that are not above the threshold.

        Calculated the matched filter, threshold, and cluster.

        Parameters
        ----------
        segnum : int
            Index into the list of segments at MatchedFilterControl construction
            against which to filter.
        template_norm : float
            The htilde, template normalization factor.
        window : int
            Size of the window over which to cluster triggers, in samples.
            This is IGNORED by this function, and provided only for API compatibility.

        Returns
        -------
        snr : TimeSeries
            A time series containing the complex snr.
        norm : float
            The normalization of the complex snr.
        corrrelation: FrequencySeries
            A frequency series containing the correlation vector.
        idx : Array
            List of indices of the triggers.
        snrv : Array
            The snr values at the trigger locations.
        """
        norm = (4.0 * self.delta_f) / sqrt(template_norm)
        self.correlators[segnum].correlate()
        self.ifft.execute()
        idx, snrv = events.threshold_only(self.snr_mem[self.segments[segnum].analyze],
                                          self.snr_threshold / norm)
        logging.info("%s points above threshold" % str(len(idx)))

        snr = TimeSeries(self.snr_mem, epoch=epoch, delta_t=self.delta_t, copy=False)
        corr = FrequencySeries(self.corr_mem, delta_f=self.delta_f, copy=False)
        return snr, norm, corr, idx, snrv

    def heirarchical_matched_filter_and_cluster(self, segnum, template_norm, window):
        """ Returns the complex snr timeseries, normalization of the complex snr,
        the correlation vector frequency series, the list of indices of the
        triggers, and the snr values at the trigger locations. Returns empty
        lists for these for points that are not above the threshold.

        Calculated the matched filter, threshold, and cluster.

        Parameters
        ----------
        segnum : int
            Index into the list of segments at MatchedFilterControl construction
        template_norm : float
            The htilde, template normalization factor.
        window : int
            Size of the window over which to cluster triggers, in samples

        Returns
        -------
        snr : TimeSeries
            A time series containing the complex snr at the reduced sample rate.
        norm : float
            The normalization of the complex snr.
        corrrelation: FrequencySeries
            A frequency series containing the correlation vector.
        idx : Array
            List of indices of the triggers.
        snrv : Array
            The snr values at the trigger locations.
        """
        from pycbc.fft.fftw_pruned import pruned_c2cifft, fft_transpose
        htilde = self.htilde
        stilde = self.segments[segnum]

        norm = (4.0 * stilde.delta_f) / sqrt(template_norm)

        correlate(htilde[self.kmin_red:self.kmax_red],
                  stilde[self.kmin_red:self.kmax_red],
                  self.corr_mem[self.kmin_red:self.kmax_red])

        ifft(self.corr_mem, self.snr_mem)

        if not hasattr(stilde, 'red_analyze'):
            stilde.red_analyze = \
                             slice(stilde.analyze.start/self.downsample_factor,
                                   stilde.analyze.stop/self.downsample_factor)


        idx_red, snrv_red = events.threshold(self.snr_mem[stilde.red_analyze],
                                self.snr_threshold / norm * self.upsample_threshold)
        if len(idx_red) == 0:
            return [], None, [], [], []

        idx_red, _ = events.cluster_reduce(idx_red, snrv_red, window / self.downsample_factor)
        logging.info("%s points above threshold at reduced resolution"\
                      %(str(len(idx_red)),))

        # The fancy upsampling is here
        if self.upsample_method=='pruned_fft':
            idx = (idx_red + stilde.analyze.start/self.downsample_factor)\
                   * self.downsample_factor

            idx = smear(idx, self.downsample_factor)

            # cache transposed  versions of htilde and stilde
            if not hasattr(self.corr_mem_full, 'transposed'):
                self.corr_mem_full.transposed = zeros(len(self.corr_mem_full), dtype=self.dtype)

            if not hasattr(htilde, 'transposed'):
                htilde.transposed = zeros(len(self.corr_mem_full), dtype=self.dtype)
                htilde.transposed[self.kmin_full:self.kmax_full] = htilde[self.kmin_full:self.kmax_full]
                htilde.transposed = fft_transpose(htilde.transposed)

            if not hasattr(stilde, 'transposed'):
                stilde.transposed = zeros(len(self.corr_mem_full), dtype=self.dtype)
                stilde.transposed[self.kmin_full:self.kmax_full] = stilde[self.kmin_full:self.kmax_full]
                stilde.transposed = fft_transpose(stilde.transposed)

            correlate(htilde.transposed, stilde.transposed, self.corr_mem_full.transposed)
            snrv = pruned_c2cifft(self.corr_mem_full.transposed, self.inter_vec, idx, pretransposed=True)
            idx = idx - stilde.analyze.start
            idx2, snrv = events.threshold(Array(snrv, copy=False), self.snr_threshold / norm)

            if len(idx2) > 0:
                correlate(htilde[self.kmax_red:self.kmax_full],
                          stilde[self.kmax_red:self.kmax_full],
                          self.corr_mem_full[self.kmax_red:self.kmax_full])
                idx, snrv = events.cluster_reduce(idx[idx2], snrv, window)
            else:
                idx, snrv = [], []

            logging.info("%s points at full rate and clustering" % len(idx))
            return self.snr_mem, norm, self.corr_mem_full, idx, snrv
        else:
            raise ValueError("Invalid upsample method")


def compute_max_snr_over_sky_loc_stat(hplus, hcross, hphccorr,
                                                      hpnorm=None, hcnorm=None,
                                                      out=None, thresh=0,
                                                      analyse_slice=None):
    """
    Compute the maximized over sky location statistic.

    Parameters
    -----------
    hplus : TimeSeries
        This is the IFFTed complex SNR time series of (h+, data). If not
        normalized, supply the normalization factor so this can be done!
        It is recommended to normalize this before sending through this
        function
    hcross : TimeSeries
        This is the IFFTed complex SNR time series of (hx, data). If not
        normalized, supply the normalization factor so this can be done!
    hphccorr : float
        The real component of the overlap between the two polarizations
        Re[(h+, hx)]. Note that the imaginary component does not enter the
        detection statistic. This must be normalized and is sign-sensitive.
    thresh : float
        Used for optimization. If we do not care about the value of SNR
        values below thresh we can calculate a quick statistic that will
        always overestimate SNR and then only calculate the proper, more
        expensive, statistic at points where the quick SNR is above thresh.
    hpsigmasq : float
        The normalization factor (h+, h+). Default = None (=1, already
        normalized)
    hcsigmasq : float
        The normalization factor (hx, hx). Default = None (=1, already
        normalized)
    out : TimeSeries (optional, default=None)
        If given, use this array to store the output.

    Returns
    --------
    det_stat : TimeSeries
        The SNR maximized over sky location
    """
    # NOTE: Not much optimization has been done here! This may need to be
    # Cythonized.

    if out is None:
        out = zeros(len(hplus))
        out.non_zero_locs = numpy.array([], dtype=out.dtype)
    else:
        if not hasattr(out, 'non_zero_locs'):
            # Doing this every time is not a zero-cost operation
            out.data[:] = 0
            out.non_zero_locs = numpy.array([], dtype=out.dtype)
        else:
            # Only set non zero locations to zero
            out.data[out.non_zero_locs] = 0


    # If threshold is given we can limit the points at which to compute the
    # full statistic
    if thresh:
        # This is the statistic that always overestimates the SNR...
        # It allows some unphysical freedom that the full statistic does not
        idx_p, _ = events.threshold_only(hplus[analyse_slice],
                                                    thresh / (2**0.5 * hpnorm))
        idx_c, _ = events.threshold_only(hcross[analyse_slice],
                                                    thresh / (2**0.5 * hcnorm))
        idx_p = idx_p + analyse_slice.start
        idx_c = idx_c + analyse_slice.start
        hp_red = hplus[idx_p] * hpnorm
        hc_red = hcross[idx_p] * hcnorm
        stat_p = hp_red.real**2 + hp_red.imag**2 + \
                     hc_red.real**2 + hc_red.imag**2
        locs_p = idx_p[stat_p > (thresh*thresh)]
        hp_red = hplus[idx_c] * hpnorm
        hc_red = hcross[idx_c] * hcnorm
        stat_c = hp_red.real**2 + hp_red.imag**2 + \
                     hc_red.real**2 + hc_red.imag**2
        locs_c = idx_c[stat_c > (thresh*thresh)]
        locs = numpy.unique(numpy.concatenate((locs_p, locs_c)))

        hplus = hplus[locs]
        hcross = hcross[locs]

    hplus = hplus * hpnorm
    hcross = hcross * hcnorm


    # Calculate and sanity check the denominator
    denom = 1 - hphccorr*hphccorr
    if denom < 0:
        if hphccorr > 1:
            err_msg = "Overlap between hp and hc is given as %f. " %(hphccorr)
            err_msg += "How can an overlap be bigger than 1?"
            raise ValueError(err_msg)
        else:
            err_msg = "There really is no way to raise this error!?! "
            err_msg += "If you're seeing this, it is bad."
            raise ValueError(err_msg)
    if denom == 0:
        # This case, of hphccorr==1, makes the statistic degenerate
        # This case should not physically be possible luckily.
        err_msg = "You have supplied a real overlap between hp and hc of 1. "
        err_msg += "Ian is reasonably certain this is physically impossible "
        err_msg += "so why are you seeing this?"
        raise ValueError(err_msg)

    assert(len(hplus) == len(hcross))

    # Now the stuff where comp. cost may be a problem
    hplus_magsq = numpy.real(hplus) * numpy.real(hplus) + \
                       numpy.imag(hplus) * numpy.imag(hplus)
    hcross_magsq = numpy.real(hcross) * numpy.real(hcross) + \
                       numpy.imag(hcross) * numpy.imag(hcross)
    rho_pluscross = numpy.real(hplus) * numpy.real(hcross) + numpy.imag(hplus)*numpy.imag(hcross)

    sqroot = (hplus_magsq - hcross_magsq)**2
    sqroot += 4 * (hphccorr * hplus_magsq - rho_pluscross) * \
                  (hphccorr * hcross_magsq - rho_pluscross)
    # Sometimes this can be less than 0 due to numeric imprecision, catch this.
    if (sqroot < 0).any():
        indices = numpy.arange(len(sqroot))[sqroot < 0]
        # This should not be *much* smaller than 0 due to numeric imprecision
        if (sqroot[indices] < -0.0001).any():
            err_msg = "Square root has become negative. Something wrong here!"
            raise ValueError(err_msg)
        sqroot[indices] = 0
    sqroot = numpy.sqrt(sqroot)
    det_stat_sq = 0.5 * (hplus_magsq + hcross_magsq - \
                         2 * rho_pluscross*hphccorr + sqroot) / denom

    det_stat = numpy.sqrt(det_stat_sq)

    if thresh:
        out.data[locs] = det_stat
        out.non_zero_locs = locs
        return out
    else:
        return Array(det_stat, copy=False)

def compute_u_val_for_sky_loc_stat(hplus, hcross, hphccorr,
                                 hpnorm=None, hcnorm=None, indices=None):
    """The max-over-sky location detection statistic maximizes over a phase,
    an amplitude and the ratio of F+ and Fx, encoded in a variable called u.
    Here we return the value of u for the given indices.
    """
    if indices is not None:
        hplus = hplus[indices]
        hcross = hcross[indices]

    if hpnorm is not None:
        hplus = hplus * hpnorm
    if hcnorm is not None:
        hcross = hcross * hcnorm

    # Sanity checking in func. above should already have identified any points
    # which are bad, and should be used to construct indices for input here
    hplus_magsq = numpy.real(hplus) * numpy.real(hplus) + \
                       numpy.imag(hplus) * numpy.imag(hplus)
    hcross_magsq = numpy.real(hcross) * numpy.real(hcross) + \
                       numpy.imag(hcross) * numpy.imag(hcross)
    rho_pluscross = numpy.real(hplus) * numpy.real(hcross) + \
                       numpy.imag(hplus)*numpy.imag(hcross)

    a = hphccorr * hplus_magsq - rho_pluscross
    b = hplus_magsq - hcross_magsq
    c = rho_pluscross - hphccorr * hcross_magsq

    sq_root = b*b - 4*a*c
    sq_root = sq_root**0.5
    sq_root = -sq_root
    # Catch the a->0 case
    bad_lgc = (a == 0)
    dbl_bad_lgc = numpy.logical_and(c == 0, b == 0)
    dbl_bad_lgc = numpy.logical_and(bad_lgc, dbl_bad_lgc)
    # Initialize u
    u = sq_root * 0.
    # In this case u is completely degenerate, so set it to 1
    u[dbl_bad_lgc] = 1.
    # If a->0 avoid overflow by just setting to a large value
    u[bad_lgc & ~dbl_bad_lgc] = 1E17
    # Otherwise normal statistic
    u[~bad_lgc] = (-b[~bad_lgc] + sq_root[~bad_lgc]) / (2*a[~bad_lgc])

    snr_cplx = hplus * u + hcross
    coa_phase = numpy.angle(snr_cplx)

    return u, coa_phase

def compute_max_snr_over_sky_loc_stat_no_phase(hplus, hcross, hphccorr,
                                               hpnorm=None, hcnorm=None,
                                               out=None, thresh=0,
                                               analyse_slice=None):
    """
    Compute the match maximized over polarization phase.

    In contrast to compute_max_snr_over_sky_loc_stat_no_phase this function
    performs no maximization over orbital phase, treating that as an intrinsic
    parameter. In the case of aligned-spin 2,2-mode only waveforms, this
    collapses to the normal statistic (at twice the computational cost!)

    Parameters
    -----------
    hplus : TimeSeries
        This is the IFFTed complex SNR time series of (h+, data). If not
        normalized, supply the normalization factor so this can be done!
        It is recommended to normalize this before sending through this
        function
    hcross : TimeSeries
        This is the IFFTed complex SNR time series of (hx, data). If not
        normalized, supply the normalization factor so this can be done!
    hphccorr : float
        The real component of the overlap between the two polarizations
        Re[(h+, hx)]. Note that the imaginary component does not enter the
        detection statistic. This must be normalized and is sign-sensitive.
    thresh : float
        Used for optimization. If we do not care about the value of SNR
        values below thresh we can calculate a quick statistic that will
        always overestimate SNR and then only calculate the proper, more
        expensive, statistic at points where the quick SNR is above thresh.
    hpsigmasq : float
        The normalization factor (h+, h+). Default = None (=1, already
        normalized)
    hcsigmasq : float
        The normalization factor (hx, hx). Default = None (=1, already
        normalized)
    out : TimeSeries (optional, default=None)
        If given, use this array to store the output.

    Returns
    --------
    det_stat : TimeSeries
        The SNR maximized over sky location
    """
    # NOTE: Not much optimization has been done here! This may need to be
    # Cythonized.

    if out is None:
        out = zeros(len(hplus))
        out.non_zero_locs = numpy.array([], dtype=out.dtype)
    else:
        if not hasattr(out, 'non_zero_locs'):
            # Doing this every time is not a zero-cost operation
            out.data[:] = 0
            out.non_zero_locs = numpy.array([], dtype=out.dtype)
        else:
            # Only set non zero locations to zero
            out.data[out.non_zero_locs] = 0

    # If threshold is given we can limit the points at which to compute the
    # full statistic
    if thresh:
        # This is the statistic that always overestimates the SNR...
        # It allows some unphysical freedom that the full statistic does not
        #
        # For now this is copied from the max-over-phase statistic. One could
        # probably make this faster by removing the imaginary components of
        # the matched filter, as these are not used here.
        idx_p, _ = events.threshold_only(hplus[analyse_slice],
                                                    thresh / (2**0.5 * hpnorm))
        idx_c, _ = events.threshold_only(hcross[analyse_slice],
                                                    thresh / (2**0.5 * hcnorm))
        idx_p = idx_p + analyse_slice.start
        idx_c = idx_c + analyse_slice.start
        hp_red = hplus[idx_p] * hpnorm
        hc_red = hcross[idx_p] * hcnorm
        stat_p = hp_red.real**2 + hp_red.imag**2 + \
                     hc_red.real**2 + hc_red.imag**2
        locs_p = idx_p[stat_p > (thresh*thresh)]
        hp_red = hplus[idx_c] * hpnorm
        hc_red = hcross[idx_c] * hcnorm
        stat_c = hp_red.real**2 + hp_red.imag**2 + \
                     hc_red.real**2 + hc_red.imag**2
        locs_c = idx_c[stat_c > (thresh*thresh)]
        locs = numpy.unique(numpy.concatenate((locs_p, locs_c)))

        hplus = hplus[locs]
        hcross = hcross[locs]

    hplus = hplus * hpnorm
    hcross = hcross * hcnorm


    # Calculate and sanity check the denominator
    denom = 1 - hphccorr*hphccorr
    if denom < 0:
        if hphccorr > 1:
            err_msg = "Overlap between hp and hc is given as %f. " %(hphccorr)
            err_msg += "How can an overlap be bigger than 1?"
            raise ValueError(err_msg)
        else:
            err_msg = "There really is no way to raise this error!?! "
            err_msg += "If you're seeing this, it is bad."
            raise ValueError(err_msg)
    if denom == 0:
        # This case, of hphccorr==1, makes the statistic degenerate
        # This case should not physically be possible luckily.
        err_msg = "You have supplied a real overlap between hp and hc of 1. "
        err_msg += "Ian is reasonably certain this is physically impossible "
        err_msg += "so why are you seeing this?"
        raise ValueError(err_msg)

    assert(len(hplus) == len(hcross))

    # Now the stuff where comp. cost may be a problem
    hplus_magsq = numpy.real(hplus) * numpy.real(hplus)
    hcross_magsq = numpy.real(hcross) * numpy.real(hcross)
    rho_pluscross = numpy.real(hplus) * numpy.real(hcross)

    det_stat_sq = (hplus_magsq + hcross_magsq - 2 * rho_pluscross*hphccorr)

    det_stat = numpy.sqrt(det_stat_sq / denom)

    if thresh:
        out.data[locs] = det_stat
        out.non_zero_locs = locs
        return out
    else:
        return Array(det_stat, copy=False)

def compute_u_val_for_sky_loc_stat_no_phase(hplus, hcross, hphccorr,
                                 hpnorm=None , hcnorm=None, indices=None):
    """The max-over-sky location (no phase) detection statistic maximizes over
    an amplitude and the ratio of F+ and Fx, encoded in a variable called u.
    Here we return the value of u for the given indices.


    """
    if indices is not None:
        hplus = hplus[indices]
        hcross = hcross[indices]

    if hpnorm is not None:
        hplus = hplus * hpnorm
    if hcnorm is not None:
        hcross = hcross * hcnorm

    rhoplusre=numpy.real(hplus)
    rhocrossre=numpy.real(hcross)
    overlap=numpy.real(hphccorr)

    denom = (-rhocrossre+overlap*rhoplusre)
    # Initialize tan_kappa array
    u_val = denom * 0.
    # Catch the denominator -> 0 case
    bad_lgc = (denom == 0)
    u_val[bad_lgc] = 1E17
    # Otherwise do normal statistic
    u_val[~bad_lgc] = (-rhoplusre+overlap*rhocrossre) / \
        (-rhocrossre+overlap*rhoplusre)
    coa_phase = numpy.zeros(len(indices), dtype=numpy.float32)

    return u_val, coa_phase


class MatchedFilterSkyMaxControl(object):
    # FIXME: This seems much more simplistic than the aligned-spin class.
    #        E.g. no correlators. Is this worth updating?
    def __init__(self, low_frequency_cutoff, high_frequency_cutoff,
                snr_threshold, tlen, delta_f, dtype):
        """
        Create a matched filter engine.

        Parameters
        ----------
        low_frequency_cutoff : {None, float}, optional
            The frequency to begin the filter calculation. If None, begin
            at the first frequency after DC.
        high_frequency_cutoff : {None, float}, optional
            The frequency to stop the filter calculation. If None, continue
            to the nyquist frequency.
        snr_threshold : float
            The minimum snr to return when filtering
        """
        self.tlen = tlen
        self.delta_f = delta_f
        self.dtype = dtype
        self.snr_threshold = snr_threshold
        self.flow = low_frequency_cutoff
        self.fhigh = high_frequency_cutoff

        self.matched_filter_and_cluster = \
                                    self.full_matched_filter_and_cluster
        self.snr_plus_mem = zeros(self.tlen, dtype=self.dtype)
        self.corr_plus_mem = zeros(self.tlen, dtype=self.dtype)
        self.snr_cross_mem = zeros(self.tlen, dtype=self.dtype)
        self.corr_cross_mem = zeros(self.tlen, dtype=self.dtype)
        self.snr_mem = zeros(self.tlen, dtype=self.dtype)
        self.cached_hplus_hcross_correlation = None
        self.cached_hplus_hcross_hplus = None
        self.cached_hplus_hcross_hcross = None
        self.cached_hplus_hcross_psd = None


    def full_matched_filter_and_cluster(self, hplus, hcross, hplus_norm,
                                        hcross_norm, psd, stilde, window):
        """
        Return the complex snr and normalization.

        Calculated the matched filter, threshold, and cluster.

        Parameters
        ----------
        h_quantities : Various
            FILL ME IN
        stilde : FrequencySeries
            The strain data to be filtered.
        window : int
            The size of the cluster window in samples.

        Returns
        -------
        snr : TimeSeries
            A time series containing the complex snr.
        norm : float
            The normalization of the complex snr.
        correlation: FrequencySeries
            A frequency series containing the correlation vector.
        idx : Array
            List of indices of the triggers.
        snrv : Array
            The snr values at the trigger locations.
        """

        I_plus, Iplus_corr, Iplus_norm = matched_filter_core(hplus, stilde,
                                          h_norm=hplus_norm,
                                          low_frequency_cutoff=self.flow,
                                          high_frequency_cutoff=self.fhigh,
                                          out=self.snr_plus_mem,
                                          corr_out=self.corr_plus_mem)


        I_cross, Icross_corr, Icross_norm = matched_filter_core(hcross,
                                          stilde, h_norm=hcross_norm,
                                          low_frequency_cutoff=self.flow,
                                          high_frequency_cutoff=self.fhigh,
                                          out=self.snr_cross_mem,
                                          corr_out=self.corr_cross_mem)

        # The information on the complex side of this overlap is important
        # we may want to use this in the future.
        if not id(hplus) == self.cached_hplus_hcross_hplus:
            self.cached_hplus_hcross_correlation = None
        if not id(hcross) == self.cached_hplus_hcross_hcross:
            self.cached_hplus_hcross_correlation = None
        if not id(psd) == self.cached_hplus_hcross_psd:
            self.cached_hplus_hcross_correlation = None
        if self.cached_hplus_hcross_correlation is None:
            hplus_cross_corr = overlap_cplx(hplus, hcross, psd=psd,
                                           low_frequency_cutoff=self.flow,
                                           high_frequency_cutoff=self.fhigh,
                                           normalized=False)
            hplus_cross_corr = numpy.real(hplus_cross_corr)
            hplus_cross_corr = hplus_cross_corr / (hcross_norm*hplus_norm)**0.5
            self.cached_hplus_hcross_correlation = hplus_cross_corr
            self.cached_hplus_hcross_hplus = id(hplus)
            self.cached_hplus_hcross_hcross = id(hcross)
            self.cached_hplus_hcross_psd = id(psd)
        else:
            hplus_cross_corr = self.cached_hplus_hcross_correlation

        snr = self._maximized_snr(I_plus,I_cross,
                                  hplus_cross_corr,
                                  hpnorm=Iplus_norm,
                                  hcnorm=Icross_norm,
                                  out=self.snr_mem,
                                  thresh=self.snr_threshold,
                                  analyse_slice=stilde.analyze)
        # FIXME: This should live further down
        # Convert output to pycbc TimeSeries
        delta_t = 1.0 / (self.tlen * stilde.delta_f)

        snr = TimeSeries(snr, epoch=stilde.start_time, delta_t=delta_t,
                         copy=False)

        idx, snrv = events.threshold_real_numpy(snr[stilde.analyze],
                                                self.snr_threshold)

        if len(idx) == 0:
            return [], 0, 0, [], [], [], [], 0, 0, 0
        logging.info("%s points above threshold", str(len(idx)))


        idx, snrv = events.cluster_reduce(idx, snrv, window)
        logging.info("%s clustered points", str(len(idx)))
        # erased self.
        u_vals, coa_phase = self._maximized_extrinsic_params\
            (I_plus.data, I_cross.data, hplus_cross_corr,
             indices=idx+stilde.analyze.start, hpnorm=Iplus_norm,
             hcnorm=Icross_norm)



        return snr, Iplus_corr, Icross_corr, idx, snrv, u_vals, coa_phase,\
                                      hplus_cross_corr, Iplus_norm, Icross_norm

    def _maximized_snr(self, hplus, hcross, hphccorr, **kwargs):
        return compute_max_snr_over_sky_loc_stat(hplus, hcross, hphccorr,
                                                 **kwargs)

    def _maximized_extrinsic_params(self, hplus, hcross, hphccorr, **kwargs):
        return compute_u_val_for_sky_loc_stat(hplus, hcross, hphccorr,
                                              **kwargs)


class MatchedFilterSkyMaxControlNoPhase(MatchedFilterSkyMaxControl):
    # Basically the same as normal SkyMaxControl, except we use a slight
    # variation in the internal SNR functions.
    def _maximized_snr(self, hplus, hcross, hphccorr, **kwargs):
        return compute_max_snr_over_sky_loc_stat_no_phase(hplus, hcross,
                                                          hphccorr, **kwargs)

    def _maximized_extrinsic_params(self, hplus, hcross, hphccorr, **kwargs):
        return compute_u_val_for_sky_loc_stat_no_phase(hplus, hcross, hphccorr,
                                                       **kwargs)

def make_frequency_series(vec):
    """Return a frequency series of the input vector.

    If the input is a frequency series it is returned, else if the input
    vector is a real time series it is fourier transformed and returned as a
    frequency series.

    Parameters
    ----------
    vector : TimeSeries or FrequencySeries

    Returns
    -------
    Frequency Series: FrequencySeries
        A frequency domain version of the input vector.
    """
    if isinstance(vec, FrequencySeries):
        return vec
    if isinstance(vec, TimeSeries):
        N = len(vec)
        n = N // 2 + 1
        delta_f = 1.0 / N / vec.delta_t
        vectilde =  FrequencySeries(zeros(n, dtype=complex_same_precision_as(vec)),
                                    delta_f=delta_f, copy=False)
        fft(vec, vectilde)
        return vectilde
    else:
        raise TypeError("Can only convert a TimeSeries to a FrequencySeries")

def sigmasq_series(htilde, psd=None, low_frequency_cutoff=None,
            high_frequency_cutoff=None):
    """Return a cumulative sigmasq frequency series.

    Return a frequency series containing the accumulated power in the input
    up to that frequency.

    Parameters
    ----------
    htilde : TimeSeries or FrequencySeries
        The input vector
    psd : {None, FrequencySeries}, optional
        The psd used to weight the accumulated power.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin accumulating power. If None, start at the beginning
        of the vector.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop considering accumulated power. If None, continue
        until the end of the input vector.

    Returns
    -------
    Frequency Series: FrequencySeries
        A frequency series containing the cumulative sigmasq.
    """
    htilde = make_frequency_series(htilde)
    N = (len(htilde)-1) * 2
    norm = 4.0 * htilde.delta_f
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, htilde.delta_f, N)

    sigma_vec = FrequencySeries(zeros(len(htilde), dtype=real_same_precision_as(htilde)),
                                delta_f = htilde.delta_f, copy=False)

    mag = htilde.squared_norm()

    if psd is not None:
        mag /= psd

    sigma_vec[kmin:kmax] = mag[kmin:kmax].cumsum()

    return sigma_vec*norm


def sigmasq(htilde, psd = None, low_frequency_cutoff=None,
            high_frequency_cutoff=None):
    """Return the loudness of the waveform. This is defined (see Duncan
    Brown's thesis) as the unnormalized matched-filter of the input waveform,
    htilde, with itself. This quantity is usually referred to as (sigma)^2
    and is then used to normalize matched-filters with the data.

    Parameters
    ----------
    htilde : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    psd : {None, FrequencySeries}, optional
        The psd used to weight the accumulated power.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin considering waveform power.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop considering waveform power.

    Returns
    -------
    sigmasq: float
    """
    htilde = make_frequency_series(htilde)
    N = (len(htilde)-1) * 2
    norm = 4.0 * htilde.delta_f
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, htilde.delta_f, N)
    ht = htilde[kmin:kmax]

    if psd:
        try:
            numpy.testing.assert_almost_equal(ht.delta_f, psd.delta_f)
        except AssertionError:
            raise ValueError('Waveform does not have same delta_f as psd')

    if psd is None:
        sq = ht.inner(ht)
    else:
        sq = ht.weighted_inner(ht, psd[kmin:kmax])

    return sq.real * norm

def sigma(htilde, psd = None, low_frequency_cutoff=None,
        high_frequency_cutoff=None):
    """ Return the sigma of the waveform. See sigmasq for more details.

    Parameters
    ----------
    htilde : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    psd : {None, FrequencySeries}, optional
        The psd used to weight the accumulated power.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin considering waveform power.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop considering waveform power.

    Returns
    -------
    sigmasq: float
    """
    return sqrt(sigmasq(htilde, psd, low_frequency_cutoff, high_frequency_cutoff))

def get_cutoff_indices(flow, fhigh, df, N):
    """
    Gets the indices of a frequency series at which to stop an overlap
    calculation.

    Parameters
    ----------
    flow: float
        The frequency (in Hz) of the lower index.
    fhigh: float
        The frequency (in Hz) of the upper index.
    df: float
        The frequency step (in Hz) of the frequency series.
    N: int
        The number of points in the **time** series. Can be odd
        or even.

    Returns
    -------
    kmin: int
    kmax: int
    """
    if flow:
        kmin = int(flow / df)
        if kmin < 0:
            err_msg = "Start frequency cannot be negative. "
            err_msg += "Supplied value and kmin {} and {}".format(flow, kmin)
            raise ValueError(err_msg)
    else:
        kmin = 1
    if fhigh:
        kmax = int(fhigh / df )
        if kmax > int((N + 1)/2.):
            kmax = int((N + 1)/2.)
    else:
        # int() truncates towards 0, so this is
        # equivalent to the floor of the float
        kmax = int((N + 1)/2.)

    if kmax <= kmin:
        err_msg = "Kmax cannot be less than or equal to kmin. "
        err_msg += "Provided values of freqencies (min,max) were "
        err_msg += "{} and {} ".format(flow, fhigh)
        err_msg += "corresponding to (kmin, kmax) of "
        err_msg += "{} and {}.".format(kmin, kmax)
        raise ValueError(err_msg)

    return kmin,kmax

def matched_filter_core(template, data, psd=None, low_frequency_cutoff=None,
                  high_frequency_cutoff=None, h_norm=None, out=None, corr_out=None):
    """ Return the complex snr and normalization.

    Return the complex snr, along with its associated normalization of the template,
    matched filtered against the data.

    Parameters
    ----------
    template : TimeSeries or FrequencySeries
        The template waveform
    data : TimeSeries or FrequencySeries
        The strain data to be filtered.
    psd : {FrequencySeries}, optional
        The noise weighting of the filter.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin the filter calculation. If None, begin at the
        first frequency after DC.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop the filter calculation. If None, continue to the
        the nyquist frequency.
    h_norm : {None, float}, optional
        The template normalization. If none, this value is calculated internally.
    out : {None, Array}, optional
        An array to use as memory for snr storage. If None, memory is allocated
        internally.
    corr_out : {None, Array}, optional
        An array to use as memory for correlation storage. If None, memory is allocated
        internally. If provided, management of the vector is handled externally by the
        caller. No zero'ing is done internally.

    Returns
    -------
    snr : TimeSeries
        A time series containing the complex snr.
    corrrelation: FrequencySeries
        A frequency series containing the correlation vector.
    norm : float
        The normalization of the complex snr.
    """
    htilde = make_frequency_series(template)
    stilde = make_frequency_series(data)

    if len(htilde) != len(stilde):
        raise ValueError("Length of template and data must match")

    N = (len(stilde)-1) * 2
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, stilde.delta_f, N)

    if corr_out is not None:
        qtilde = corr_out
    else:
        qtilde = zeros(N, dtype=complex_same_precision_as(data))

    if out is None:
        _q = zeros(N, dtype=complex_same_precision_as(data))
    elif (len(out) == N) and type(out) is Array and out.kind =='complex':
        _q = out
    else:
        raise TypeError('Invalid Output Vector: wrong length or dtype')

    correlate(htilde[kmin:kmax], stilde[kmin:kmax], qtilde[kmin:kmax])

    if psd is not None:
        if isinstance(psd, FrequencySeries):
            try:
                numpy.testing.assert_almost_equal(stilde.delta_f, psd.delta_f)
            except AssertionError:
                raise ValueError("PSD delta_f does not match data")
            qtilde[kmin:kmax] /= psd[kmin:kmax]
        else:
            raise TypeError("PSD must be a FrequencySeries")

    ifft(qtilde, _q)

    if h_norm is None:
        h_norm = sigmasq(htilde, psd, low_frequency_cutoff, high_frequency_cutoff)

    norm = (4.0 * stilde.delta_f) / sqrt( h_norm)

    return (TimeSeries(_q, epoch=stilde._epoch, delta_t=stilde.delta_t, copy=False),
           FrequencySeries(qtilde, epoch=stilde._epoch, delta_f=stilde.delta_f, copy=False),
           norm)

def smear(idx, factor):
    """
    This function will take as input an array of indexes and return every
    unique index within the specified factor of the inputs.

    E.g.: smear([5,7,100],2) = [3,4,5,6,7,8,9,98,99,100,101,102]

    Parameters
    -----------
    idx : numpy.array of ints
        The indexes to be smeared.
    factor : idx
        The factor by which to smear out the input array.

    Returns
    --------
    new_idx : numpy.array of ints
        The smeared array of indexes.
    """


    s = [idx]
    for i in range(factor+1):
        a = i - factor/2
        s += [idx + a]
    return numpy.unique(numpy.concatenate(s))

def matched_filter(template, data, psd=None, low_frequency_cutoff=None,
                  high_frequency_cutoff=None, sigmasq=None):
    """ Return the complex snr.

    Return the complex snr, along with its associated normalization of the
    template, matched filtered against the data.

    Parameters
    ----------
    template : TimeSeries or FrequencySeries
        The template waveform
    data : TimeSeries or FrequencySeries
        The strain data to be filtered.
    psd : FrequencySeries
        The noise weighting of the filter.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin the filter calculation. If None, begin at the
        first frequency after DC.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop the filter calculation. If None, continue to the
        the nyquist frequency.
    sigmasq : {None, float}, optional
        The template normalization. If none, this value is calculated
        internally.

    Returns
    -------
    snr : TimeSeries
        A time series containing the complex snr.
    """
    snr, _, norm = matched_filter_core(template, data, psd=psd,
            low_frequency_cutoff=low_frequency_cutoff,
            high_frequency_cutoff=high_frequency_cutoff, h_norm=sigmasq)
    return snr * norm

_snr = None
def match(vec1, vec2, psd=None, low_frequency_cutoff=None,
          high_frequency_cutoff=None, v1_norm=None, v2_norm=None):
    """ Return the match between the two TimeSeries or FrequencySeries.

    Return the match between two waveforms. This is equivalent to the overlap
    maximized over time and phase.

    Parameters
    ----------
    vec1 : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    vec2 : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    psd : Frequency Series
        A power spectral density to weight the overlap.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin the match.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop the match.
    v1_norm : {None, float}, optional
        The normalization of the first waveform. This is equivalent to its
        sigmasq value. If None, it is internally calculated.
    v2_norm : {None, float}, optional
        The normalization of the second waveform. This is equivalent to its
        sigmasq value. If None, it is internally calculated.

    Returns
    -------
    match: float
    index: int
        The number of samples to shift to get the match.
    """

    htilde = make_frequency_series(vec1)
    stilde = make_frequency_series(vec2)

    N = (len(htilde)-1) * 2

    global _snr
    if _snr is None or _snr.dtype != htilde.dtype or len(_snr) != N:
        _snr = zeros(N,dtype=complex_same_precision_as(vec1))
    snr, _, snr_norm = matched_filter_core(htilde,stilde,psd,low_frequency_cutoff,
                             high_frequency_cutoff, v1_norm, out=_snr)
    maxsnr, max_id = snr.abs_max_loc()
    if v2_norm is None:
        v2_norm = sigmasq(stilde, psd, low_frequency_cutoff, high_frequency_cutoff)
    return maxsnr * snr_norm / sqrt(v2_norm), max_id

def overlap(vec1, vec2, psd=None, low_frequency_cutoff=None,
          high_frequency_cutoff=None, normalized=True):
    """ Return the overlap between the two TimeSeries or FrequencySeries.

    Parameters
    ----------
    vec1 : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    vec2 : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    psd : Frequency Series
        A power spectral density to weight the overlap.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin the overlap.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop the overlap.
    normalized : {True, boolean}, optional
        Set if the overlap is normalized. If true, it will range from 0 to 1.

    Returns
    -------
    overlap: float
    """

    return overlap_cplx(vec1, vec2, psd=psd, \
            low_frequency_cutoff=low_frequency_cutoff,\
            high_frequency_cutoff=high_frequency_cutoff,\
            normalized=normalized).real

def overlap_cplx(vec1, vec2, psd=None, low_frequency_cutoff=None,
          high_frequency_cutoff=None, normalized=True):
    """Return the complex overlap between the two TimeSeries or FrequencySeries.

    Parameters
    ----------
    vec1 : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    vec2 : TimeSeries or FrequencySeries
        The input vector containing a waveform.
    psd : Frequency Series
        A power spectral density to weight the overlap.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin the overlap.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop the overlap.
    normalized : {True, boolean}, optional
        Set if the overlap is normalized. If true, it will range from 0 to 1.

    Returns
    -------
    overlap: complex
    """
    htilde = make_frequency_series(vec1)
    stilde = make_frequency_series(vec2)

    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
            high_frequency_cutoff, stilde.delta_f, (len(stilde)-1) * 2)

    if psd:
        inner = (htilde[kmin:kmax]).weighted_inner(stilde[kmin:kmax], psd[kmin:kmax])
    else:
        inner = (htilde[kmin:kmax]).inner(stilde[kmin:kmax])

    if normalized:
        sig1 = sigma(vec1, psd=psd, low_frequency_cutoff=low_frequency_cutoff,
                     high_frequency_cutoff=high_frequency_cutoff)
        sig2 = sigma(vec2, psd=psd, low_frequency_cutoff=low_frequency_cutoff,
                     high_frequency_cutoff=high_frequency_cutoff)
        norm = 1 / sig1 / sig2
    else:
        norm = 1

    return 4 * htilde.delta_f * inner * norm

def quadratic_interpolate_peak(left, middle, right):
    """ Interpolate the peak and offset using a quadratic approximation

    Parameters
    ----------
    left : numpy array
        Values at a relative bin value of [-1]
    middle : numpy array
        Values at a relative bin value of [0]
    right : numpy array
        Values at a relative bin value of [1]

    Returns
    -------
    bin_offset : numpy array
        Array of bins offsets, each in the range [-1/2, 1/2]
    peak_values : numpy array
        Array of the estimated peak values at the interpolated offset
    """
    bin_offset = 1.0/2.0 * (left - right) / (left - 2 * middle + right)
    peak_value = middle + 0.25 * (left - right) * bin_offset
    return bin_offset, peak_value


class LiveBatchMatchedFilter(object):

    """Calculate SNR and signal consistency tests in a batched progression"""

    def __init__(self, templates, snr_threshold, chisq_bins, sg_chisq,
                 maxelements=2**27,
                 snr_abort_threshold=None,
                 newsnr_threshold=None,
                 max_triggers_in_batch=None):
        """Create a batched matchedfilter instance

        Parameters
        ----------
        templates: list of `FrequencySeries`
            List of templates from the FilterBank class.
        snr_threshold: float
            Minimum value to record peaks in the SNR time series.
        chisq_bins: str
            Str that determines how the number of chisq bins varies as a
            function of the template bank parameters.
        sg_chisq: pycbc.vetoes.SingleDetSGChisq
            Instance of the sg_chisq class to calculate sg_chisq with.
        maxelements: {int, 2**27}
            Maximum size of a batched fourier transform.
        snr_abort_threshold: {float, None}
            If the SNR is above this threshold, do not record any triggers.
        newsnr_threshold: {float, None}
            Only record triggers that have a re-weighted NewSNR above this
        threshold.
        max_triggers_in_batch: {int, None}
            Record X number of the loudest triggers by newsnr in each mpi
        process group. Signal consistency values will also only be calculated
        for these triggers.
        """
        self.snr_threshold = snr_threshold
        self.snr_abort_threshold = snr_abort_threshold
        self.newsnr_threshold = newsnr_threshold
        self.max_triggers_in_batch = max_triggers_in_batch

        from pycbc import vetoes
        self.power_chisq = vetoes.SingleDetPowerChisq(chisq_bins, None)
        self.sg_chisq = sg_chisq

        durations = numpy.array([1.0 / t.delta_f for t in templates])

        lsort = durations.argsort()
        durations = durations[lsort]
        templates = [templates[li] for li in lsort]

        # Figure out how to chunk together the templates into groups to process
        _, counts = numpy.unique(durations, return_counts=True)
        tsamples = [(len(t) - 1) * 2 for t in templates]
        grabs = maxelements / numpy.unique(tsamples)

        chunks = numpy.array([])
        num = 0
        for count, grab in zip(counts, grabs):
            chunks = numpy.append(chunks, numpy.arange(num, count + num, grab))
            chunks = numpy.append(chunks, [count + num])
            num += count
        chunks = numpy.unique(chunks).astype(numpy.uint32)

        # We now have how many templates to grab at a time.
        self.chunks = chunks[1:] - chunks[0:-1]

        self.out_mem = {}
        self.cout_mem = {}
        self.ifts = {}
        chunk_durations = [durations[i] for i in chunks[:-1]]
        self.chunk_tsamples = [tsamples[int(i)] for i in chunks[:-1]]
        samples = self.chunk_tsamples * self.chunks

        # Create workspace memory for correlate and snr
        mem_ids = [(a, b) for a, b in zip(chunk_durations, self.chunks)]
        mem_types = set(zip(mem_ids, samples))

        self.tgroups, self.mids = [], []
        for i, size in mem_types:
            dur, count = i
            self.out_mem[i] = zeros(size, dtype=numpy.complex64)
            self.cout_mem[i] = zeros(size, dtype=numpy.complex64)
            self.ifts[i] = IFFT(self.cout_mem[i], self.out_mem[i],
                                nbatch=count,
                                size=len(self.cout_mem[i]) // count)

        # Split the templates into their processing groups
        for dur, count in mem_ids:
            tgroup = templates[0:count]
            self.tgroups.append(tgroup)
            self.mids.append((dur, count))
            templates = templates[count:]

        # Associate the snr and corr memory block to each template
        self.corr = []
        for i, tgroup in enumerate(self.tgroups):
            psize = self.chunk_tsamples[i]
            s = 0
            e = psize
            mid = self.mids[i]
            for htilde in tgroup:
                htilde.out = self.out_mem[mid][s:e]
                htilde.cout = self.cout_mem[mid][s:e]
                s += psize
                e += psize
            self.corr.append(BatchCorrelator(tgroup, [t.cout for t in tgroup], len(tgroup[0])))

    def set_data(self, data):
        """Set the data reader object to use"""
        self.data = data
        self.block_id = 0

    def combine_results(self, results):
        """Combine results from different batches of filtering"""
        result = {}
        for key in results[0]:
            result[key] = numpy.concatenate([r[key] for r in results])
        return result

    def process_data(self, data_reader):
        """Process the data for all of the templates"""
        self.set_data(data_reader)
        return self.process_all()

    def process_all(self):
        """Process every batch group and return as single result"""
        results = []
        veto_info = []
        while 1:
            result, veto = self._process_batch()
            if result is False: return False
            if result is None: break
            results.append(result)
            veto_info += veto

        result = self.combine_results(results)

        if self.max_triggers_in_batch:
            sort = result['snr'].argsort()[::-1][:self.max_triggers_in_batch]
            for key in result:
                result[key] = result[key][sort]

            tmp = veto_info
            veto_info = [tmp[i] for i in sort]

        result = self._process_vetoes(result, veto_info)
        return result

    def _process_vetoes(self, results, veto_info):
        """Calculate signal based vetoes"""
        chisq = numpy.array(numpy.zeros(len(veto_info)), numpy.float32, ndmin=1)
        dof = numpy.array(numpy.zeros(len(veto_info)), numpy.uint32, ndmin=1)
        sg_chisq = numpy.array(numpy.zeros(len(veto_info)), numpy.float32,
                               ndmin=1)
        results['chisq'] = chisq
        results['chisq_dof'] = dof
        results['sg_chisq'] = sg_chisq

        keep = []
        for i, (snrv, norm, l, htilde, stilde) in enumerate(veto_info):
            correlate(htilde, stilde, htilde.cout)
            c, d = self.power_chisq.values(htilde.cout, snrv,
                                           norm, stilde.psd, [l], htilde)
            chisq[i] = c[0] / d[0]
            dof[i] = d[0]

            sgv = self.sg_chisq.values(stilde, htilde, stilde.psd,
                                       snrv, norm, c, d, [l])
            if sgv is not None:
                sg_chisq[i] = sgv[0]

            if self.newsnr_threshold:
                newsnr = ranking.newsnr(results['snr'][i], chisq[i])
                if newsnr >= self.newsnr_threshold:
                    keep.append(i)

        if self.newsnr_threshold:
            keep = numpy.array(keep, dtype=numpy.uint32)
            for key in results:
                results[key] = results[key][keep]

        return results

    def _process_batch(self):
        """Process only a single batch group of data"""
        if self.block_id == len(self.tgroups):
            return None, None

        tgroup = self.tgroups[self.block_id]
        psize = self.chunk_tsamples[self.block_id]
        mid = self.mids[self.block_id]
        stilde = self.data.overwhitened_data(tgroup[0].delta_f)
        psd = stilde.psd

        valid_end = int(psize - self.data.trim_padding)
        valid_start = int(valid_end - self.data.blocksize * self.data.sample_rate)

        seg = slice(valid_start, valid_end)

        self.corr[self.block_id].execute(stilde)
        self.ifts[mid].execute()

        self.block_id += 1

        snr = numpy.zeros(len(tgroup), dtype=numpy.complex64)
        time = numpy.zeros(len(tgroup), dtype=numpy.float64)
        templates = numpy.zeros(len(tgroup), dtype=numpy.uint64)
        sigmasq = numpy.zeros(len(tgroup), dtype=numpy.float32)

        time[:] = self.data.start_time

        result = {}
        tkeys = tgroup[0].params.dtype.names
        for key in tkeys:
            result[key] = []

        veto_info = []

        # Find the peaks in our SNR times series from the various templates
        i = 0
        for htilde in tgroup:
            l = htilde.out[seg].abs_arg_max()

            sgm = htilde.sigmasq(psd)
            norm = 4.0 * htilde.delta_f / (sgm ** 0.5)

            l += valid_start
            snrv = numpy.array([htilde.out[l]])

            # If nothing is above threshold we can exit this template
            s = abs(snrv[0]) * norm
            if s < self.snr_threshold:
                continue

            time[i] += float(l - valid_start) / self.data.sample_rate

            # We have an SNR so high that we will drop the entire analysis
            # of this chunk of time!
            if self.snr_abort_threshold is not None and s > self.snr_abort_threshold:
                logging.info("We are seeing some *really* high SNRs, lets"
                             " assume they aren't signals and just give up")
                return False, []

            veto_info.append((snrv, norm, l, htilde, stilde))

            snr[i] = snrv[0] * norm
            sigmasq[i] = sgm
            templates[i] = htilde.id
            if not hasattr(htilde, 'dict_params'):
                htilde.dict_params = {}
                for key in tkeys:
                    htilde.dict_params[key] = htilde.params[key]

            for key in tkeys:
                result[key].append(htilde.dict_params[key])
            i += 1


        result['snr'] = abs(snr[0:i])
        result['coa_phase'] = numpy.angle(snr[0:i])
        result['end_time'] = time[0:i]
        result['template_id'] = templates[0:i]
        result['sigmasq'] = sigmasq[0:i]

        for key in tkeys:
            result[key] = numpy.array(result[key])

        return result, veto_info

def followup_event_significance(ifo, data_reader, bank,
                                template_id, coinc_times,
                                coinc_threshold=0.005,
                                lookback=150, duration=0.095):
    """ Followup an event in another detector and determine its significance
    """
    from pycbc.waveform import get_waveform_filter_length_in_time
    tmplt = bank.table[template_id]
    length_in_time = get_waveform_filter_length_in_time(tmplt['approximant'],
                                                        tmplt)

    # calculate onsource time range
    from pycbc.detector import Detector
    onsource_start = -numpy.inf
    onsource_end = numpy.inf
    fdet = Detector(ifo)

    for cifo in coinc_times:
        time = coinc_times[cifo]
        dtravel =  Detector(cifo).light_travel_time_to_detector(fdet)
        if time - dtravel > onsource_start:
            onsource_start = time - dtravel
        if time + dtravel < onsource_end:
            onsource_end = time + dtravel

    # Source must be within this time window to be considered a possible
    # coincidence
    onsource_start -= coinc_threshold
    onsource_end += coinc_threshold

    # Calculate how much time needed to calculate significance
    trim_pad = (data_reader.trim_padding * data_reader.strain.delta_t)
    bdur = int(lookback + 2.0 * trim_pad + length_in_time)
    if bdur > data_reader.strain.duration * .75:
        bdur = data_reader.strain.duration * .75

    # Require all strain be valid within lookback time
    if data_reader.state is not None:
        state_start_time = data_reader.strain.end_time \
                - data_reader.reduced_pad * data_reader.strain.delta_t - bdur
        if not data_reader.state.is_extent_valid(state_start_time, bdur):
            return None, None, None, None

    # We won't require that all DQ checks be valid for now, except at
    # onsource time.
    if data_reader.dq is not None:
        dq_start_time = onsource_start - duration / 2.0
        dq_duration = onsource_end - onsource_start + duration
        if not data_reader.dq.is_extent_valid(dq_start_time, dq_duration):
            return None, None, None, None

    # Calculate SNR time series for this duration
    htilde = bank.get_template(template_id, min_buffer=bdur)
    stilde = data_reader.overwhitened_data(htilde.delta_f)

    sigma2 = htilde.sigmasq(stilde.psd)
    snr, _, norm = matched_filter_core(htilde, stilde, h_norm=sigma2)

    # Find peak in on-source and determine p-value
    onsrc = snr.time_slice(onsource_start, onsource_end)
    peak = onsrc.abs_arg_max()
    peak_time = peak * snr.delta_t + onsrc.start_time
    peak_value = abs(onsrc[peak])

    bstart = float(snr.start_time) + length_in_time + trim_pad
    bkg = abs(snr.time_slice(bstart, onsource_start)).numpy()

    window = int((onsource_end - onsource_start) * snr.sample_rate)
    nsamples = int(len(bkg) / window)

    peaks = bkg[:nsamples*window].reshape(nsamples, window).max(axis=1)
    pvalue = (1 + (peaks >= peak_value).sum()) / float(1 + nsamples)

    # Return recentered source SNR for bayestar, along with p-value, and trig
    baysnr = snr.time_slice(peak_time - duration / 2.0,
                            peak_time + duration / 2.0)

    logging.info('Adding %s to candidate, pvalue %s, %s samples', ifo,
                 pvalue, nsamples)

    return baysnr * norm, peak_time, pvalue, sigma2

def compute_followup_snr_series(data_reader, htilde, trig_time,
                                duration=0.095, check_state=True,
                                coinc_window=0.05):
    """Given a StrainBuffer, a template frequency series and a trigger time,
    compute a portion of the SNR time series centered on the trigger for its
    rapid sky localization and followup.

    If the trigger time is too close to the boundary of the valid data segment
    the SNR series is calculated anyway and might be slightly contaminated by
    filter and wrap-around effects. For reasonable durations this will only
    affect a small fraction of the triggers and probably in a negligible way.

    Parameters
    ----------
    data_reader : StrainBuffer
        The StrainBuffer object to read strain data from.

    htilde : FrequencySeries
        The frequency series containing the template waveform.

    trig_time : {float, lal.LIGOTimeGPS}
        The trigger time.

    duration : float (optional)
        Duration of the computed SNR series in seconds. If omitted, it defaults
        to twice the Earth light travel time plus 10 ms of timing uncertainty.

    check_state : boolean
        If True, and the detector was offline or flagged for bad data quality
        at any point during the inspiral, then return (None, None) instead.

    coinc_window : float (optional)
        Maximum possible time between coincident triggers at different
        detectors. This is needed to properly determine data padding.

    Returns
    -------
    snr : TimeSeries
        The portion of SNR around the trigger. None if the detector is offline
        or has bad data quality, and check_state is True.
    """
    if check_state:
        # was the detector observing for the full amount of involved data?
        state_start_time = trig_time - duration / 2 - htilde.length_in_time
        state_end_time = trig_time + duration / 2
        state_duration = state_end_time - state_start_time
        if data_reader.state is not None:
            if not data_reader.state.is_extent_valid(state_start_time,
                                                     state_duration):
                return None

        # was the data quality ok for the full amount of involved data?
        dq_start_time = state_start_time - data_reader.dq_padding
        dq_duration = state_duration + 2 * data_reader.dq_padding
        if data_reader.dq is not None:
            if not data_reader.dq.is_extent_valid(dq_start_time, dq_duration):
                return None

    stilde = data_reader.overwhitened_data(htilde.delta_f)
    snr, _, norm = matched_filter_core(htilde, stilde,
                                          h_norm=htilde.sigmasq(stilde.psd))

    valid_end = int(len(snr) - data_reader.trim_padding)
    valid_start = int(valid_end - data_reader.blocksize * snr.sample_rate)

    half_dur_samples = int(snr.sample_rate * duration / 2)
    coinc_samples = int(snr.sample_rate * coinc_window)
    valid_start -= half_dur_samples + coinc_samples
    valid_end += half_dur_samples
    if valid_start < 0 or valid_end > len(snr)-1:
        raise ValueError(('Requested SNR duration ({0} s)'
                          ' too long').format(duration))

    # Onsource slice for Bayestar followup
    onsource_idx = float(trig_time - snr.start_time) * snr.sample_rate
    onsource_idx = int(round(onsource_idx))
    onsource_slice = slice(onsource_idx - half_dur_samples,
                           onsource_idx + half_dur_samples + 1)
    return snr[onsource_slice] * norm

__all__ = ['match', 'matched_filter', 'sigmasq', 'sigma', 'get_cutoff_indices',
           'sigmasq_series', 'make_frequency_series', 'overlap',
           'overlap_cplx', 'matched_filter_core', 'correlate',
           'MatchedFilterControl', 'LiveBatchMatchedFilter',
           'MatchedFilterSkyMaxControl', 'MatchedFilterSkyMaxControlNoPhase',
           'compute_max_snr_over_sky_loc_stat_no_phase',
           'compute_max_snr_over_sky_loc_stat',
           'compute_followup_snr_series',
           'compute_u_val_for_sky_loc_stat_no_phase',
           'compute_u_val_for_sky_loc_stat',
           'followup_event_significance']
