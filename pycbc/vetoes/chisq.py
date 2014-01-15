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
import numpy
import logging
import pycbc.fft

from pycbc.types import Array, zeros, real_same_precision_as, TimeSeries, complex_same_precision_as, FrequencySeries
from pycbc.filter import sigmasq_series, make_frequency_series, sigmasq, matched_filter_core
from pycbc.scheme import schemed

BACKEND_PREFIX="pycbc.vetoes.chisq_"

def power_chisq_bins_from_sigmasq_series(sigmasq_series, num_bins, kmin, kmax):
    """Returns bins of equal power for use with the chisq functions
    
    Parameters
    ----------
    
    sigmasq_series: FrequencySeries
        A frequency series containing the cumulative power of a filter template
        preweighted by a psd.
    num_bins: int
    kmin: int
        The number of chisq bins to calculate.
    kmax: int    
        
    Returns
    -------
    
    bins: List of ints
        A list of the edges of the chisq bins is returned. 
    
    """
    sigmasq = sigmasq_series[kmax - 1]                 
    edge_vec = numpy.arange(0, num_bins) * sigmasq / num_bins
    bins = numpy.searchsorted(sigmasq_series[kmin:kmax], edge_vec, side='right')
    bins += kmin
    return numpy.append(bins, kmax)

def power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff=None, 
                     high_frequency_cutoff=None):
    """Returns bins of equal power for use with the chisq functions
    
    Parameters
    ----------
    
    htilde: FrequencySeries
        A frequency series containing the template waveform
    num_bins: int
        The number of chisq bins to calculate.
    psd: FrequencySeries
        A frequency series containing the psd. Its length must be commensurate
        with the template waveform.
    low_frequency_cutoff: {None, float}, optional
        The low frequency cutoff to apply
    high_frequency_cutoff: {None, float}, optional
        The high frequency cutoff to apply
    
    Returns
    -------
    
    bins: List of ints
        A list of the edges of the chisq bins is returned. 
    """
    sigma_vec = sigmasq_series(htilde, psd, low_frequency_cutoff, 
                               high_frequency_cutoff).numpy() 
    kmin = int(low_frequency_cutoff / htilde.delta_f)
    kmax = len(sigma_vec) - 1
    return power_chisq_bins_from_sigmasq_series(sigma_vec, num_bins, kmin, kmax)
    
@schemed(BACKEND_PREFIX)
def chisq_accum_bin(chisq, q):
    pass
   
@schemed(BACKEND_PREFIX) 
def shift_sum(v1, shifts, slen=None, offset=0):
    """ Calculate the time shifted sum of the FrequencySeries
    """
    pass
 
def power_chisq_at_points_from_precomputed(corr, snr, snr_norm, bins, indices):
    """Calculate the chisq timeseries from precomputed values for only select
    points.
    
    This function calculates the chisq at each point by explicitly time shifting
    and summing each bin. No FFT is involved.
    
    Parameters
    ----------  
    corr: FrequencySeries
        The product of the template and data in the frequency domain.
    snr: Array
        The unnormalized array of snr values at only the selected points in `indices`.
    norm: float
        The normalization of the snr
    bins: List of integers
        The edges of the equal power bins
    indices: Array
        The indices where we will calculate the chisq. These must be relative 
        to the given `corr` series.           
    
    Returns
    -------
    chisq: Array
        An array containing only the chisq at the selected points.
    """
    snr = Array(snr, copy=False)   
    chisq = zeros(len(indices), dtype=real_same_precision_as(corr))     
    num_bins = len(bins) - 1
    
    for j in range(len(bins)-1):
        k_min = int(bins[j])
        k_max = int(bins[j+1])                 
        qi = shift_sum(corr[k_min:k_max], indices, slen=len(corr), offset=k_min)
        chisq_accum_bin(chisq, qi)
        
    return (chisq * num_bins - snr.squared_norm()) * (snr_norm ** 2.0)
    
def power_chisq_from_precomputed(corr, snr, snr_norm, bins):
    """Calculate the chisq timeseries from precomputed values
    
    This function calculates the chisq at all times by performing an 
    inverse FFT of each bin.
    
    Parameters
    ----------
    
    corr: FrequencySeries
        The produce of the template and data in the frequency domain.
    snr: TimeSeries
        The unnormalized snr time series.
    snr_norm:
        The snr normalization factor. (true snr = snr * snr_norm)
    bins: List of integers
        The edges of the chisq bins.   
    
    Returns
    -------
    chisq: TimeSeries
    """      
    q = zeros(len(snr), dtype=complex_same_precision_as(snr))
    qtilde = zeros(len(snr), dtype=complex_same_precision_as(snr))
    chisq = TimeSeries(zeros(len(snr), dtype=real_same_precision_as(snr)), 
                       delta_t=snr.delta_t, epoch=snr.start_time, copy=False)
    
    chisq_norm = snr_norm ** 2.0
    num_bins = len(bins) - 1
    
    for j in range(len(bins)-1): 
        k_min = int(bins[j])
        k_max = int(bins[j+1])
        
        qtilde[k_min:k_max] = corr[k_min:k_max]
        pycbc.fft.ifft(qtilde, q) 
        qtilde[k_min:k_max].clear()
        chisq_accum_bin(chisq, q)
        
    return (chisq * num_bins - snr.squared_norm()) * chisq_norm
    
def fastest_power_chisq_at_points(corr, snr, snr_norm, bins, indices):
    """Calculate the chisq values for only select points.
    
    This function looks at the number of point that need to evaluated and
    selects the fastest method (FFT, or direct time shift and sum). In either
    case, only the selected points are returned.
    
    Parameters
    ----------
    
    corr: FrequencySeries
        The product of the template and data in the frequency domain.
    snr: Array
        The unnormalized snr
    snr_norm: float
        The snr normalization factor
    bins: List of integers
        The edges of the equal power bins
    indices: Array
        The indices where we will calculate the chisq. These must be relative 
        to the given `snr` series.           
    
    Returns
    -------
    chisq: Array
        An array containing only the chisq at the selected points.
    """ 
    # Time per fft (ms)
    FFT_T = 50
    # Time per single point calculation
    P_T = 3
    # This is empirically chosen from tests on SUGAR. It may not be correct
    # into the future. Replace with better estimate or auto-tuning.
    POINT_THRESHOLD = (len(bins)-1)*FFT_T / P_T
    
    # This is a temporary hack, standard gpu support is intended and should be forthcoming
    # Note, the code for gpu support is in place but atm very slow, optimization required
    # before turning it on by default
    import pycbc.scheme
    if (len(indices) < POINT_THRESHOLD and type(pycbc.scheme.mgr.state) is pycbc.scheme.CPUScheme) :
        # We don't have that many points so do the direct time shift.
        return power_chisq_at_points_from_precomputed(corr, snr.take(indices), snr_norm, bins, indices)
    else:
        # We have a lot of points so it is faster to use the fourier transform
        return power_chisq_from_precomputed(corr, snr, snr_norm, bins).take(indices)    

def power_chisq(template, data, num_bins, psd, low_frequency_cutoff=None, high_frequency_cutoff=None):
    """Calculate the chisq timeseries 

    Parameters
    ----------
    template: FrequencySeries or TimeSeries
        A time or frequency series that contains the filter template. The length
        must be commensurate with the data. 
    data: FrequencySeries or TimeSeries
        A time ore frequency series that contains the data to filter. The length
        must be commensurate with the template.
    num_bins: int
        The number of bins in the chisq. Note that the dof goes as 2*num_bin-2. 
    psd: FrequencySeries
        The psd of the data. 
    low_frequency_cutoff: {None, float}, optional
        The low frequency cutoff to apply.
    high_frequency_cutoff: {None, float}, optional
        The high frequency cutoff to apply. 

    Returns
    -------
    chisq: TimeSeries
        TimeSeries containing the chisq values for all times. 
    """
    htilde = make_frequency_series(template)
    stilde = make_frequency_series(data)   
     
    bins = power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff, high_frequency_cutoff)   
    corra = zeros((len(htilde)-1)*2, dtype=htilde.dtype)   
    total_snr, corr, tnorm = matched_filter_core(htilde, stilde, psd,
                           low_frequency_cutoff, high_frequency_cutoff, corr_out=corra)

    return power_chisq_from_precomputed(corr, total_snr, tnorm, bins)

class SingleDetPowerChisq(object):
    """Class that handles precomputation and memory management for efficiently
    running the power chisq in a single detector inspiral analysis.
    """
    def __init__(self, num_bins):
        if num_bins > 0:
            self.do = True
            self.column_name = "chisq"
            self.table_dof_name = "chisq_dof"
            self.dof = num_bins * 2 - 2
            
            self._num_bins = num_bins
            self._bins = None
            self._template = None
        else:
            self.do = False           

    def values(self, corr, snr, snr_norm, psd, indices, template, bank, low_frequency_cutoff):
        if self.do:
            # Compute the chisq bins if we haven't already
            # Only recompute the bins if the template changes
            
            if self._template is None or self._template != template:        
                if bank.sigmasq_vec is not None:
                    logging.info("...Calculating fast power chisq bins")
                    kmin = int(low_frequency_cutoff / corr.delta_f)
                    kmax = template.end_idx
                    bins = power_chisq_bins_from_sigmasq_series(bank.sigmasq_vec, self._num_bins, kmin, kmax)
                else:  
                    logging.info("...Calculating power chisq bins")
                    bins = power_chisq_bins(template, self._num_bins, psd, low_frequency_cutoff)
                self._template = template
                self._bins = bins
                
            logging.info("...Doing power chisq")     
            return fastest_power_chisq_at_points(corr, snr, snr_norm, self._bins, indices)
