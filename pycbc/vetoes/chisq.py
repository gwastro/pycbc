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
from pycbc.types import zeros, real_same_precision_as, TimeSeries, complex_same_precision_as, FrequencySeries
from pycbc.filter import sigmasq_series, make_frequency_series, sigmasq, matched_filter_core
import numpy
from pycbc.scheme import schemed
import pycbc.fft

BACKEND_PREFIX="pycbc.vetoes.chisq_"

def power_chisq_bins_from_sigmasq_series(sigmasq_series, num_bins, kmin, kmax):
    """Returns bins of equal power for use with the chisq functions
    """
    sigmasq = sigmasq_series[kmax - 1]                        
    edge_vec = numpy.arange(0, num_bins) * sigmasq / num_bins
    bins = numpy.searchsorted(sigmasq_series[kmin:kmax], edge_vec, side='right')
    bins += kmin
    return numpy.append(bins, kmax)

def power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff=None, 
                     high_frequency_cutoff=None):
    """Returns bins of equal power for use with the chisq functions
    """
    sigma_vec = sigmasq_series(htilde, psd, low_frequency_cutoff, 
                               high_frequency_cutoff).numpy() 
    kmin = int(low_frequency_cutoff / htilde.delta_f)
    kmax = len(sigma_vec)
    return power_chisq_bins_from_sigmasq_series(sigma_vec, num_bins, kmin, kmax)
    
@schemed(BACKEND_PREFIX)
def chisq_accum_bin(chisq, snrp, q):
    pass
    
def power_chisq_from_precomputed(corr, snr, bins, snr_norm):
    """ Returns the chisq time series
    """
    q = zeros(len(snr), dtype=complex_same_precision_as(snr))
    qtilde = zeros(len(snr), dtype=complex_same_precision_as(snr))
    chisq = TimeSeries(zeros(len(snr), dtype=real_same_precision_as(snr)), 
                       delta_t=snr.delta_t, epoch=snr.start_time, copy=False)
    
    chisq_norm = snr_norm ** 2.0
    num_bins = len(bins) - 1
    
    snrp = snr/num_bins
    
    for j in range(len(bins)-1): 
        k_min = int(bins[j])
        k_max = int(bins[j+1])
        
        qtilde[k_min:k_max] = corr[k_min:k_max]
        pycbc.fft.ifft(qtilde, q) 
        qtilde[k_min:k_max].clear()
        chisq_accum_bin(chisq, snrp, q)
        
    return chisq * (num_bins * chisq_norm)

def power_chisq(template, data, num_bins, psd, low_frequency_cutoff=None, high_frequency_cutoff=None):
    """ Return the chisq time series.
    """  
    htilde = make_frequency_series(template)
    stilde = make_frequency_series(data)
    
    bins = power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff, high_frequency_cutoff)
    
    corra = zeros((len(htilde)-1)*2, dtype=htilde.dtype)
    
    total_snr, corr, tnorm = matched_filter_core(htilde, stilde, psd,
                           low_frequency_cutoff, high_frequency_cutoff, corr_out=corra)

    return power_chisq_from_precomputed(corr, total_snr, bins, tnorm)

    
                
        


