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
This modules provides matchedfiltering and related match and chisq calculations.
"""

from pycbc.types import *
from pycbc.fft import fft,ifft
from math import log,ceil,sqrt


def get_padded_frequencyseries(vec):
    if not isinstance(vec,TimeSeries):
        raise TypeError("Input must be a Timeseries")
    else:
        power = ceil(log(len(vec),2))+1
        N = 2 ** power
        n = N/2+1
        
        vec_pad = TimeSeries(zeros(N),delta_t=vec.delta_t,dtype=real_same_precision_as(vec))
        vec_pad[0:len(vec)] = vec
        
        vectilde = FrequencySeries(zeros(n),delta_f=1, dtype=complex_same_precision_as(vec))
        
        fft(vec_pad,vectilde)
        
        return vectilde

def get_frequencyseries(vec):
    if isinstance(vec,FrequencySeries):
        return vec
    if isinstance(vec,TimeSeries):
        N = len(vec)
        n = N/2+1    
        delta_f = 1.0 / N / vec.delta_t
        vectilde = FrequencySeries(zeros(n),delta_f=delta_f, dtype=complex_same_precision_as(vec))
        fft(vec,vectilde)   
        return vectilde
    else:
        raise TypeError("Can only convert a TimeSeries to a FrequencySeries")



def sigmasq_series(htilde,psd = None,low_frequency_cutoff=None,high_frequency_cutoff=None):    
    N = (len(htilde)-1) * 2 
    norm = 4.0 / (N * N * htilde.delta_f) 
    kmin,kmax = get_cutoff_indices(low_frequency_cutoff,high_frequency_cutoff,htilde.delta_f,N)   
    
    moment = htilde.squared_norm()
    
    if psd is not None:
        moment[kmin:kmax] /= psd[kmin:kmax]      
        
    return moment[kmin:kmax],norm

def sigmasq(htilde,psd = None,low_frequency_cutoff=None,high_frequency_cutoff=None):
    moment,norm = sigmasq_series(htilde,psd,low_frequency_cutoff,high_frequency_cutoff)
    return moment.sum() * norm
    
def get_cutoff_indices(flow,fhigh,df,N):
    if flow:
        kmin = int(flow / df)
    else:
        kmin = 1
    if fhigh:
        kmax = int(fhigh / df )
    else:
        kmax = N/2 + 1
        
    return kmin,kmax
    

_q = None
_qtilde = None

def matchedfilter(template,data,psd=None,low_frequency_cutoff=None,high_frequency_cutoff=None):
    global _q
    global _qtilde
  
    # Get the Inputs in terms of htilde and stilde
    htilde = get_frequencyseries(template)
    stilde = get_frequencyseries(data)

    # Calculate the length we need for the temporary memory 
    # Check that this is a power of two?
    N = (len(htilde)-1) * 2   
    kmin,kmax = get_cutoff_indices(low_frequency_cutoff,high_frequency_cutoff,stilde.delta_f,N) 
   
    # Create workspace memory
    if (_q is None) or (len(_q) != N):
        _q = zeros(N,dtype=complex_same_precision_as(data))
    else:
        _q.fill(0)
        
    if (_qtilde is None) or (len(_q) != N):
        _qtilde = zeros(N,dtype=complex_same_precision_as(data))
    else:
        _qtilde.fill(0)

    correlate(htilde[kmin:kmax],stilde[kmin:kmax],_qtilde[kmin:kmax])
    if psd is not None:
        _qtilde[kmin:kmax] /= psd[kmin:kmax]

    ifft(_qtilde,_q) 

    norm = sqrt(((4.0 / (N * N * stilde.delta_f)) **2) / sigmasq(htilde,psd,low_frequency_cutoff,high_frequency_cutoff) )
    
    #return complex snr
    return _q,norm
    
    
def match(vec1,vec2,psd=None,low_frequency_cutoff=None,high_frequency_cutoff=None):
    htilde = get_frequencyseries(vec1)
    stilde = get_frequencyseries(vec2)
    snr,norm = matchedfilter(htilde,stilde,psd,low_frequency_cutoff,high_frequency_cutoff)

    maxsnrsq = (snr.squared_norm()).max()
    return sqrt(maxsnrsq/sigmasq(stilde,psd,low_frequency_cutoff,high_frequency_cutoff))*norm
    
def chisq(template, data,num_bins, psd = None , low_frequency_cutoff=None,high_frequency_cutoff=None):
    bins = get_chisq_bin_sizes(num_bins,template,psd,low_frequency_cutoff,high_frequency_cutoff)
 
    total_snr,norm = matchedfilter(template,data,psd,low_frequency_cutoff,high_frequency_cutoff)
    
    bin_snrs = []
    N = (len(htilde)-1) * 2   
    delta_t = 1.0 / N / data.delta_f
    
    chisq_series = TimeSeries(zeros(N),delta_t=delta_t,dtype = real_same_precision_as(data))
    
    for kmin,kmax in bins:
        template_piece = template[kmin:kmax]
        snr,part_norm = matchedfilter(template_piece,data,psd,low_frequency_cutoff,high_frequency_cutoff)
        delta = (snr/num_bins - total_snr)
        chisq_series += (delta.conj()*delta)
        
    chisq_series *= norm * num_bins
    return chisq_series
    
   
def get_chisq_bin_sizes(num_bins,template,psd=None,low_frequency_cutoff=None,high_frequency_cutoff=None):
    pass
    
    




    
