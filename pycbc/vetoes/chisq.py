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
from pycbc.types import zeros, real_same_precision_as, TimeSeries, complex_same_precision_as
from pycbc.filter import sigmasq_series, make_frequency_series, sigmasq, matched_filter
import numpy
import pycbc.fft

def power_chisq_bins_from_sigmasq_series(sigmasq_series, num_bins, kmin, kmax):
    """Returns bins of equal power for use with the chisq functions
    """
    sigmasq = sigmasq_series[kmax]                        
    edge_vec = numpy.arange(0, sigmasq, sigmasq/num_bins)
    bins = numpy.searchsorted(sigmasq_series[kmin:kmax], edge_vec, side='right')
    bins += kmin
    return numpy.append(bins, kmax)

def power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff):
    """Returns bins of equal power for use with the chisq functions
    """
    sigma_vec = sigmasq_series(htilde, psd, low_frequency_cutoff, 
                               high_frequency_cutoff).numpy() 
    kmin = int(low_frequency_cutoff / htilde.delta_f)
    kmax = len(sigma_vec) 
    return power_chisq_bins_from_sigmasq_series(sigma_vec, num_bins, kmin, kmax)
    
def power_chisq_from_precomputed(corr, snr, bins, h_norm):
    """ Returns the chisq time series
    """
    q = zeros(len(snr), dtype=complex_same_precision_as(snr))
    qtilde = zeros(len(snr), dtype=complex_same_precision_as(snr))
    chisq = TimeSeries(zeros(len(snr), dtype=real_same_precision_as(snr)), 
                       delta_t=snr.delta_t, epoch=snr.start_time, copy=False)
    
    chisq_norm = 16 * corr.delta_f * corr.delta_f / h_norm
    num_bins = len(bins) - 1
    
    for j in range(len(bins)-1): 
        k_min = bins[j]
        k_max = bins[j+1]
        qtilde.clear()
        qtilde[k_min:k_max] = corr[k_min:k_max]
        pycbc.fft.ifft(qtilde, q) 
        chisq += (snr / num_bins - q).squared_norm()
        
    return chisq * num_bins * chisq_norm

def power_chisq(template, data, num_bins, psd=None, low_frequency_cutoff=None, 
          high_frequency_cutoff=None):
    """ Return the chisq time series.
    """  
    htilde = make_frequency_series(template)
    stilde = make_frequency_series(data)
    
    bins = power_chisq_bins(htilde, num_bins, psd, low_frequency_cutoff, high_frequency_cutoff)
    
    total_snr, tnorm = matched_filter(htilde, stilde, psd,
                           low_frequency_cutoff, high_frequency_cutoff)
    print tnorm                       
   # print sigmasq(htilde, psd=psd, low_frequency_cutoff=low_frequency_cutoff, high_frequency_cutoff=high_frequency_cutoff)

    chisq_vec = zeros(len(total_snr), dtype=real_same_precision_as(htilde))
    
    for i in range(len(bins)):
        f_min = bins[i] * template.delta_f
        
        if (i+1) < len(bins):
            f_max = bins[i+1] * template.delta_f
        else:
            f_max = high_frequency_cutoff

        snr_part, pnorm = matched_filter(htilde, stilde, psd, f_min, f_max)
       # print sigmasq(htilde, psd=psd, low_frequency_cutoff=f_min, high_frequency_cutoff=f_max)
        
        chisq_vec += (total_snr*tnorm/15 - snr_part*tnorm).squared_norm()*15
        print pnorm, snr_part[100000], f_min, f_max
    
   # print chisq_vec[100000]
    #chisq_vec *= num_bins
   # print chisq_vec[100000], total_snr.squared_norm()[100000], tnorm, tnorm*tnorm
    #chisq_vec = (chisq_vec - total_snr.squared_norm()) * tnorm * tnorm

    return TimeSeries(chisq_vec, delta_t=total_snr.delta_t, epoch=total_snr.start_time, copy=False)

    
                
        


