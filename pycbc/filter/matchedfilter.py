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
from math import log, ceil, sqrt
from pycbc.types import TimeSeries, FrequencySeries, zeros, Array, complex64
from pycbc.types import complex_same_precision_as, real_same_precision_as
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
    def __init__(self, low_frequency_cutoff, high_frequency_cutoff, 
                snr_threshold, tlen, delta_f, dtype, downsample_factor=1, 
                upsample_threshold=1, upsample_method='pruned_fft'):
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
        downsample_factor : {1, int}, optional
            The factor by which to reduce the sample rate when doing a heirarchical
            matched filter
        upsample_threshold : {1, float}, optional
            The fraction of the snr_threshold to trigger on the subsampled filter.
        upsample_method : {pruned_fft, str}
            The method to upsample or interpolate the reduced rate filter.
        """

        self.tlen = tlen
        self.delta_f = delta_f
        self.dtype = dtype
        self.snr_threshold = snr_threshold    
        self.flow = low_frequency_cutoff
        self.fhigh = high_frequency_cutoff    
                
        if downsample_factor == 1:
            self.matched_filter_and_cluster = self.full_matched_filter_and_cluster
            self.snr_mem = zeros(self.tlen, dtype=self.dtype)
            self.corr_mem = zeros(self.tlen, dtype=self.dtype)           
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
              
    def full_matched_filter_and_cluster(self, template, template_norm, stilde, window):
        """ Return the complex snr and normalization. 
    
        Calculated the matched filter, threshold, and cluster. 

        Parameters
        ----------
        htilde : FrequencySeries 
            The template waveform. Must come from the FilterBank class.
        template_norm : float
            The htilde, template normalization factor.
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
        corrrelation: FrequencySeries
            A frequency series containing the correlation vector. 
        idx : Array
            List of indices of the triggers.
        snrv : Array
            The snr values at the trigger locations.
        """
        snr, corr, norm = matched_filter_core(template, stilde, 
                                              low_frequency_cutoff=self.flow, 
                                              high_frequency_cutoff=self.fhigh, 
                                              h_norm=template_norm,
                                              out=self.snr_mem,
                                              corr_out=self.corr_mem)
        idx, snrv = events.threshold(snr[stilde.analyze], self.snr_threshold / norm)            

        if len(idx) == 0:
            return [], [], [], [], []            
        logging.info("%s points above threshold" % str(len(idx)))              
  
        idx, snrv = events.cluster_reduce(idx, snrv, window)
        logging.info("%s clustered points" % str(len(idx)))
        
        return snr, norm, corr, idx, snrv   
        
    def heirarchical_matched_filter_and_cluster(self, htilde, template_norm, stilde, window):
        """ Return the complex snr and normalization. 
    
        Calculated the matched filter, threshold, and cluster. 

        Parameters
        ----------
        htilde : FrequencySeries 
            The template waveform. Must come from the FilterBank class.
        template_norm : float
            The htilde, template normalization factor.
        stilde : FrequencySeries 
            The strain data to be filtered.
        window : int
            The size of the cluster window in samples.

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
        

__all__ = ['match', 'matched_filter', 'sigmasq', 'sigma', 'get_cutoff_indices',
           'sigmasq_series', 'make_frequency_series', 'overlap', 'overlap_cplx',
           'matched_filter_core', 'correlate', 'MatchedFilterControl']

