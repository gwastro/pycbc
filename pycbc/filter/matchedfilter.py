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
from math import log,ceil,sqrt
from pycbc.types import TimeSeries,FrequencySeries,zeros, Array, complex64
from pycbc.types import complex_same_precision_as,real_same_precision_as
from pycbc.fft import fft, ifft
import pycbc.scheme
from pycbc import events
import pycbc
import numpy
from scipy import interpolate

BACKEND_PREFIX="pycbc.filter.matchedfilter_"

@pycbc.scheme.schemed(BACKEND_PREFIX)
def correlate(x, y, z):
    pass


class MatchedFilterControl(object):
    """
    This class holds all the information/settings related to matched-filtering
    and provides a simple interface for pycbc to use when performing
    matched-filtering.
    """
    def __init__(self, lower_frequency_cutoff, upper_frequency_cutoff,
                 snr_threshold, strain_segments,
                 cluster_method=None,
                 cluster_window=0,
                 downsample_factor=None,
                 downsample_precheck_threshold=None,
                 upsample_method='pruned_fft'):
        """
        Initialize instance and set various parameters.

        Parameters
        -----------
        lower_frequency_cutoff : float
            The lower frequency cutoff to use when filtering.
        upper_frequency_cutoff : float
            The upper frequency cutoff to use when filtering.
            WARNING: This currently is not hooked up and will do nothing.
            This will be fixed soon!
        snr_threshold : float
            The lower SNR threshold at which to accept triggers.
        cluster_method : string, optional (default: None)
            The method to use to cluster over individual templates' SNR time
            series. Choices are:
            * Window: Cluster over a fixed time window (see cluster_window)
            * Template: Cluster over the length of the template
            * None: Perform no clustering; *not* recommended.
        cluster_window : float
            If using cluster_method=window, specifies the length of time, in
            seconds, to cluster over.
        downsample_factor : int, optional (default: None)
            If provided the data will be prefiltered at a lower frequency than
            that of the input data by this factor and then upsampled according
            to options below. This can provide a significant speed-up, but
            some of the options may decrease sensitivity. This has only been
            tested with binary values (ie. 2, 4, ...), it is not recommended to
            use any other value.
        downsample_precheck_threshold : float
            This is required if downsample_factor is given. This, multiplied by
            the SNR threshold, sets the threshold in the downsampled data for
            upsampling ... ie. we will not upsample parts of the SNR time
            series where SNR is less than SNR threshold multiplied by this.
            A recommended value is to use a downsample_factor of 4 and a
            downsample_precheck_threshold of 0.98. Other values/settings may
            also be possible/appropriate.
        upsample_method : string, default = "pruned_fft"
            Used only if downsample_factor is given. This specifies the method
            that will be used for upsampling the data. Current options are:
            * pruned_fft: Use Alex's fancy FFT (which he will provide more
              details of), to upsample the data with no loss in SNR. The only
              potential loss in sensitivity here is if we have not identified
              the necessary number of points to be upsampled.
            * interpolation: Use a 2nd order interpolation technique to
              evaluate SNR at the higher sample rate. As this performs an
              interpolation at the downsampled data, any higher-frequency
              content not included in the downsampled filter will be lost.
              Therefore this can reduce sensitivity for systems with power at
              high frequencies. (Alex also doesn't trust the interpolation
              completely).
        """
        # Set values
        self.f_low = lower_frequency_cutoff
        self.f_upper = upper_frequency_cutoff
        self.snr_threshold = snr_threshold
        self.segments = strain_segments.fourier_segments()
        self.cluster_method = cluster_method
        self.cluster_window = cluster_window
        self.downsample_factor=downsample_factor
        self.downsample_precheck_threshold = downsample_precheck_threshold
        self.upsample_method = upsample_method

        # Error checking
        if self.downsample_factor == 1:
            # If factor is 1, just turn it off
            self.downsample_factor = None
        elif self.downsample_factor is not None and self.downsample_factor < 1:
            # Huh?
            err_msg = "Downsample factor must be greater than one."
            err_msg += "Otherwise you are upsampling, which is a bad idea."
            raise ValueError(err_msg)
        if self.downsample_factor and ((not self.cluster_method) or\
                (self.cluster_method == 'window' and not self.cluster_window)):
            err_msg = "Downsampling requires that we are clustering during the"
            err_msg += "matched-filtering process. Please provide a cluster "
            err_msg += "method, and window if appropriate to the "
            err_msg += "MatchedFilterControl."     


        # Set derived values
        self.flen = strain_segments.freq_len
        self.tlen = strain_segments.time_len
        self.delta_f = strain_segments.delta_f
        self.delta_t = strain_segments.delta_t

        # Initialize memory (starting to look like C-code!)
        self.snr_mem = zeros(self.tlen, dtype=complex64)
        self.corr_mem = zeros(self.tlen, dtype=complex64)
        if downsample_factor:
            N2 = self.tlen / downsample_factor
            self.downsampled_snr_mem = zeros(N2, dtype=complex64)
            self.downsampled_corr_mem = zeros(N2, dtype=complex64)
            self.downsampled_pruned_mem = zeros(self.tlen, dtype=complex64)

        # FIXME: Should overwhitening occur within here?

    def matched_filter_and_cluster(self, template, stilde):
        """
        Compute the matched filter between two data streams. We assume that any
        PSD has already been included by (over-)whitening one (or both) of 
        these inputs. Also cluster output to find points above threshold and
        loudest within the cluster window.

        Parameters
        ----------
        template : PyCBC FrequencyArray or TimeArray
            The template, given either as a pycbc FrequencyArray or TimeArray.
        stilde : PyCBC FrequencyArray
            The data, usually over-whitened.

        Returns
        --------
        snr : PyCBC complex TimeArray
            The calculated complex SNR time series. This is *not* normalized,
            to get sane values, multiply the entire array by norm. For some
            methods values are only populated at a small subset of points where
            the SNR is loud/above threshold.
        norm : float
            The normalization factor of the template, which needs to be applied
            to snr as described above.
        corr : PyCBC complex FrequencyArray
            The frequency-domain complex correlation between template and
            stilde. The IFFT of this gives SNR.
        idx : Array of ints
            The locations in the snr array that are above threshold and have
            passed any clustering checks
        snrv : Array of complex
            The values of SNR at the list of idx. Also needs normalizing by the
            norm factor.
        """
        if self.cluster_method == "window":
            cluster_window = int(self.cluster_window * (1 / self.delta_t))
        elif self.cluster_method == "template":
            cluster_window = int(template.length_in_time * (1 / self.delta_t))
        else:
            cluster_window = None

        # Choose a path
        if self.downsample_factor is None:
            # This is the standard filtering engine
            snr, corr, norm = matched_filter_core(template, stilde,
                                               h_norm=template.sigmasq,
                                               low_frequency_cutoff=self.f_low,
                                               out=self.snr_mem,
                                               corr_out=self.corr_mem)
            idx, snrv = events.threshold(snr[stilde.analyze],
                                         self.snr_threshold / norm)
            if len(idx) == 0:
                return [], 0, [], [], []
            logging.info("%s points above threshold" % str(len(idx)))

            if cluster_window:
                idx, snrv = events.cluster_reduce(idx, snrv, cluster_window)
                logging.info("%s points after clustering" % str(len(idx)))
            return snr, norm, corr, idx, snrv
        else:
            snr, norm, corr, idx, snrv = \
                    dynamic_rate_thresholded_matched_filter(template,
                            stilde, template.sigmasq, self.downsample_factor,
                            self.downsample_precheck_threshold,
                            self.snr_threshold, cluster_window,
                            low_frequency_cutoff=self.f_low,
                            upsample_method=self.upsample_method,
                            snr_mem=self.snr_mem, corr_mem=self.corr_mem,
                            downsampled_snr_mem=self.downsampled_snr_mem,
                            downsampled_corr_mem=self.downsampled_corr_mem,
                            downsampled_pruned_mem=self.downsampled_pruned_mem)
            return snr, norm, corr, idx, snrv
       

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
        n = N/2+1    
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
    """Return the power of the waveform. 

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

    if psd and ht.delta_f != psd.delta_f:
        raise ValueError('Waveform does not have same delta_f as psd')

    if psd is None:
        sq = ht.inner(ht)
    else:
        sq = ht.weighted_inner(ht, psd[kmin:kmax])
        
    return sq.real * norm

def sigma(htilde, psd = None, low_frequency_cutoff=None,
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
    else:
        kmin = 1
    if fhigh:
        kmax = int(fhigh / df )
    else:
        # int() truncates towards 0, so this is
        # equivalent to the floor of the float
        kmax = int((N + 1)/2.)
        
    return kmin,kmax
    
# Workspace Memory for the matchedfilter
_qtilde_t = None

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
    if corr_out is not None:
        _qtilde = corr_out
    else:
        global _qtilde_t
        _qtilde = _qtilde_t
  
    htilde = make_frequency_series(template)
    stilde = make_frequency_series(data)

    if len(htilde) != len(stilde):
        raise ValueError("Length of template and data must match")

    N = (len(stilde)-1) * 2   
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, stilde.delta_f, N)   

    if out is None:
        _q = zeros(N, dtype=complex_same_precision_as(data))
    elif (len(out) == N) and type(out) is Array and out.kind =='complex':
        _q = out
    else:
        raise TypeError('Invalid Output Vector: wrong length or dtype')
        
    if corr_out:
        pass
    elif (_qtilde is None) or (len(_qtilde) != N) or _qtilde.dtype != data.dtype:
        _qtilde_t = _qtilde = zeros(N, dtype=complex_same_precision_as(data))
    else:
        _qtilde.clear()         
    
    correlate(htilde[kmin:kmax], stilde[kmin:kmax], _qtilde[kmin:kmax])

    if psd is not None:
        if isinstance(psd, FrequencySeries):
            if psd.delta_f == stilde.delta_f :
                _qtilde[kmin:kmax] /= psd[kmin:kmax]
            else:
                raise TypeError("PSD delta_f does not match data")
        else:
            raise TypeError("PSD must be a FrequencySeries")
            
    ifft(_qtilde, _q)
    
    if h_norm is None:
        h_norm = sigmasq(htilde, psd, low_frequency_cutoff, high_frequency_cutoff)     

    norm = (4.0 * stilde.delta_f) / sqrt( h_norm)
    delta_t = 1.0 / (N * stilde.delta_f)
    
    return (TimeSeries(_q, epoch=stilde._epoch, delta_t=delta_t, copy=False), 
           FrequencySeries(_qtilde, epoch=stilde._epoch, delta_f=htilde.delta_f, copy=False), 
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
           
def dynamic_rate_thresholded_matched_filter(htilde, stilde, h_norm,
                                            downsample_factor,
                                            downsample_threshold,
                                            snr_threshold,
                                            cluster_window,
                                            low_frequency_cutoff=None,
                                            high_frequency_cutoff=None,
                                            upsample_method='pruned_fft',
                                            snr_mem=None,
                                            corr_mem=None,
                                            downsampled_snr_mem=None,
                                            downsampled_corr_mem=None,
                                            downsampled_pruned_mem=None):
    """
    Perform a matched-filter calculation by:
    * Downsampling the data
    * Performing a matched-filter at the downsampled rate,   
    * Thresholding at a reduced threshold
    * Computing the upsampled SNR at the points above the downsampled threshold

    Parameters
    ----------
    htilde : FrequencySeries 
        The template waveform
    stilde : FrequencySeries 
        The strain data to be filtered.
    h_norm : float
        The template normalization.
        FIXME: Should be able to calculate this internally (and cache) if
        not provided.
    downsample_factor : int
        The factor by which the sample rate of the input data should be
        reduced before initial filtering. Binary values (2,4,...) are
        recommended.
    downsample_threshold : float
        The factor, multiplied by the SNR threshold, to use when thresholding
        the downsampled SNR time series.
    snr_threshold : float
        The threshold with which to accept triggers at the full sample rate.
    cluster_window : int
        The length, in sample points (at full sample rate), to use when
        performing clustering.
    psd : {FrequencySeries}, optional
        The noise weighting of the filter.
    low_frequency_cutoff : {None, float}, optional
        The frequency to begin the filter calculation. If None, begin at the
        first frequency after DC.
    high_frequency_cutoff : {None, float}, optional
        The frequency to stop the filter calculation. If None, continue to the 
        the nyquist frequency.
    upsample_method : string, default = "pruned_fft"
        Used only if downsample_factor is given. This specifies the method
        that will be used for upsampling the data. Current options are:
        * pruned_fft: Use Alex's fancy FFT (which he will provide more
          details of), to upsample the data with no loss in SNR. The only
          potential loss in sensitivity here is if we have not identified
          the necessary number of points to be upsampled.
        * interpolation: Use a 2nd order interpolation technique to
          evaluate SNR at the higher sample rate. As this performs an
          interpolation at the downsampled data, any higher-frequency
          content not included in the downsampled filter will be lost.
          Therefore this can reduce sensitivity for systems with power at
          high frequencies. (Alex also doesn't trust the interpolation
          completely).
    snr_mem : {None, Array}, optional
        An array to use as memory for snr storage. If None, memory is allocated 
        internally.
    corr_mem : {None, Array}, optional
        An array to use as memory for correlation storage. If None, memory is
        allocated internally.
    downsampled_snr_mem : {None, Array}, optional
        An array to use as memory for reduced filter SNR storage.
        If None, memory is allocated internally.
    downsampled_corr_mem : {None, Array}, optional
        An array to use as memory for reduced filter correlation storage.
        If None, memory is allocated internally.
    downsampled_pruned_mem : {None, Array}, optional
        An array to use as memory for the intermediate pruned_fft calculation.
        If None, memory is allocated internally.
    """
    valid_methods = ['pruned_fft', 'interpolation']
    if upsample_method not in valid_methods:
        err_msg = "Upsample method %s is not recognized." %(upsample_method)
        err_msg = "Supported methods are %s." %(' '.join(valid_methods))
        raise ValueError(err_msg)

    from pycbc.fft.fftw_pruned import pruned_c2cifft, fft_transpose
    from pycbc.events import threshold, cluster_reduce

    N = (len(stilde)-1) * 2   
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, stilde.delta_f, N)                                     
    N2 = N / downsample_factor
    kmin2, kmax2 = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, stilde.delta_f, N2)   
    norm = (4.0 * stilde.delta_f) / sqrt(h_norm)
    ctype = complex_same_precision_as(htilde)

    # Set up the arrays if not already provided
    if snr_mem is None:
        snr_mem = zeros(N, dtype=ctype)
    delta_t = 1.0 / (N * stilde.delta_f)
    q = TimeSeries(snr_mem, epoch=stilde._epoch, delta_t=delta_t, copy=False)
    if corr_mem is None:
        corr_mem = zeros(N, dtype=ctype)
    qtilde = FrequencySeries(corr_mem, delta_f=stilde.delta_f, copy=False)
    if downsampled_snr_mem is None:
        downsampled_snr_mem = zeros(N2, dtype=ctype)
    q2 = downsampled_snr_mem 
    if downsampled_corr_mem is None:
        downsampled_corr_mem = zeros(N2, dtype=ctype)
    qtilde2 = downsampled_corr_mem
    if downsampled_pruned_mem is None:
        downsampled_pruned_mem = zeros(N2, dtype=ctype)
    tempvec = downsampled_pruned_mem 

    correlate(htilde[kmin2:kmax2], stilde[kmin2:kmax2], qtilde2[kmin2:kmax2])    
    ifft(qtilde2, q2)    
    
    q2s = q2[stilde.analyze.start/downsample_factor:stilde.analyze.stop/downsample_factor]
    idx2, snrv2 = threshold(q2s, snr_threshold / norm * downsample_threshold)
    if len(idx2) == 0:
        return [], None, [], [], []
    # Cluster
    idx2, _ = cluster_reduce(idx2, snrv2, cluster_window / downsample_factor)
    logging.info("%s points above threshold at lower filter resolution"\
                  %(str(len(idx2)),))

    # This is a simple linear interpolation. Any content above the reduced
    # Nyquist frequency is lost, but this does have the desired time resolution
    if upsample_method=='interpolation':
        idx_shift_fac = stilde.analyze.start/downsample_factor
        snr_indexes = []
        snrv = []
        for index2 in idx2:
            interp_idxs = numpy.arange(index2-2, index2+3, dtype='int')
            try:
                interp_snrs = q2[interp_idxs + idx_shift_fac]
            except IndexError:
                # This happens when the point is at the edge of the SNR time
                # series. Here we don't have the necessary +/- 2 points on
                # either side to perform interpolation so we just fall back
                # on the pruned FFT.
                # NOTE: This should now not be possible to happen as we allow
                # points in the SNR time series that are not within
                # stilde.analyze
                upsample_method='pruned_fft'
                logging.info("Something went wrong, falling back to pruned FFT")
                break
            q_idx = index2 * downsample_factor
            interp_rsmpl_idxs = interp_idxs*downsample_factor
            q_idxs = numpy.arange(q_idx-downsample_factor/2,
                                  q_idx+downsample_factor/2 + 1, dtype='int')
            snr_indexes.extend(q_idxs)
            interp_func = interpolate.interp1d(interp_rsmpl_idxs, interp_snrs, 
                                               kind='quadratic')
            upsampled_snrs = interp_func(q_idxs)
            snrv.extend(upsampled_snrs)
            # FIXME: In numpy, this is done with a one line array operation. I
            # couldn't get that to work here th_ough!
            for i in range(len(upsampled_snrs)):
                q[q_idxs[i] + stilde.analyze.start] = upsampled_snrs[i] 
        else:
            # We end up here only if the loop above is successful. Otherwise
            # this is skipped and we fall back on the fancyfft.
            correlate(htilde[kmin:kmax], stilde[kmin:kmax], qtilde[kmin:kmax])
            snr_indexes = numpy.array(snr_indexes)
            snrv = numpy.array(snrv, dtype=numpy.complex64)
            logging.info("%s points above threshold" % str(len(snr_indexes)))
            snr_indexes, snrv = events.cluster_reduce(snr_indexes, snrv,
                                                      cluster_window)
            msg = "%s points at full filter resolution" %str(len(snr_indexes))
            msg += " after interpolation upsampling and clustering."
            logging.info(msg)
            return q, norm, qtilde, snr_indexes, snrv

    # The fancy upsampling is here
    if upsample_method=='pruned_fft':
        idx = (idx2*downsample_factor) + stilde.analyze.start
        idx = smear(idx, downsample_factor)
        correlate(htilde[kmin:kmax], stilde[kmin:kmax], qtilde[kmin:kmax])
        
        # If there are too many points, revert back to IFFT
        # FIXME: What should this value be??

        if len (idx) > 150:
            msg = "Too many points at lower sample rate, reverting to IFFT"
            logging.info(msg)
            ifft(qtilde, q)
            qs = q[stilde.analyze.start:stilde.analyze.stop]
            idx, snrv = threshold(qs, snr_threshold / norm)
            if len(idx) == 0:
                return [], None, [], [], []
            # FIXME: Use proper cluster window!
            idx, snrv = cluster_reduce(idx, snrv, cluster_window)
            msg = "%s points at full filter resolution" %(str(len(idx)),)
            msg += " after full sample rate filter and clustering."
            logging.info(msg)
            return q, norm, qtilde, idx, snrv
        # Or do the fancy upsampling
        else:             
            # cache transposed  versions of htilde and stilde
            if not hasattr(corr_mem, 'transposed'):
                corr_mem.transposed = qtilde * 1
            
            if not hasattr(htilde, 'transposed'):
                htilde.transposed = qtilde * 1
                htilde.transposed[kmin:kmax] = htilde[kmin:kmax]
                htilde.transposed = fft_transpose(htilde.transposed)
            
            if not hasattr(stilde, 'transposed'):
                stilde.transposed = qtilde * 1
                stilde.transposed[kmin:kmax] = stilde[kmin:kmax]
                stilde.transposed = fft_transpose(stilde.transposed)  
            
            correlate(htilde.transposed, stilde.transposed, corr_mem.transposed)      
            snrv = pruned_c2cifft(corr_mem.transposed, tempvec, idx, pretransposed=True)   
            q.data[idx] = snrv
            idx = idx - stilde.analyze.start
            
            msg = "%s points at full filter resolution" %(str(len(idx)),)
            msg += " after pruned FFT upsample and clustering."
            logging.info(msg)
            return q, norm, qtilde, idx, snrv

    # I shouldn't have gotten here
    err_msg = "Something went wrong somewhere. Please contact a developer."
    raise ValueError(err_msg) 
    
           
def matched_filter(template, data, psd, low_frequency_cutoff=None,
                  high_frequency_cutoff=None):
    """ Return the complex snr and normalization. 
    
    Return the complex snr, along with its associated normalization of the template,
    matched filtered against the data. 

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

    Returns
    -------
    snr : TimeSeries
        A time series containing the complex snr. 
    """
    snr, corr, norm = matched_filter_core(template, data, psd, low_frequency_cutoff, high_frequency_cutoff)
    return snr * norm
    
_snr = None 
def match(vec1, vec2, psd=None, low_frequency_cutoff=None,
          high_frequency_cutoff=None, v1_norm=None, v2_norm=None):
    """ Return the match between the two TimeSeries or FrequencySeries.
    
    Return the match between two waveforms. This is equivelant to the overlap 
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
    """

    htilde = make_frequency_series(vec1)
    stilde = make_frequency_series(vec2)

    N = (len(htilde)-1) * 2

    global _snr
    if _snr is None or _snr.dtype != htilde.dtype or len(_snr) != N:
        _snr = zeros(N,dtype=complex_same_precision_as(vec1))
    snr, corr, snr_norm = matched_filter_core(htilde,stilde,psd,low_frequency_cutoff,
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
        

__all__ = ['match', 'matched_filter', 'sigmasq', 'sigma', 'dynamic_rate_thresholded_matched_filter',
           'sigmasq_series', 'make_frequency_series', 'overlap', 'overlap_cplx',
           'matched_filter_core', 'correlate', 'MatchedFilterControl']

