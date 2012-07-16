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

from math import log,ceil,sqrt
from pycbc.types import TimeSeries,FrequencySeries,zeros
from pycbc.types import complex_same_precision_as,real_same_precision_as
from pycbc.fft import fft,ifft
import pycbc.scheme
import pycbc

def get_padded_frequencyseries(vec):
    """Pad a TimeSeries with a length of zeros greater than its length, such
    that the total length is the closest power of 2. This prevents the effects 
    of wraparound.
    """
    if not isinstance(vec,TimeSeries):
        raise TypeError("Input must be a Timeseries")
    else:
        power = ceil(log(len(vec),2))+1
        N = 2 ** power
        n = N/2+1
        
        vec_pad = TimeSeries(zeros(N),delta_t=vec.delta_t,
                             dtype=real_same_precision_as(vec))
        vec_pad[0:len(vec)] = vec
        vectilde = FrequencySeries(zeros(n),delta_f=1, 
                                   dtype=complex_same_precision_as(vec))
        
        fft(vec_pad,vectilde)
        
        return vectilde

def get_frequencyseries(vec):
    """Convenience function that returns immediately if given a FrequencySeries,
    or ffts it and returns a frequency series.
    """
    if isinstance(vec,FrequencySeries):
        return vec
    if isinstance(vec,TimeSeries):
        N = len(vec)
        n = N/2+1    
        delta_f = 1.0 / N / vec.delta_t
        vectilde = FrequencySeries(zeros(n),delta_f=delta_f, 
                                   dtype=complex_same_precision_as(vec))
        fft(vec,vectilde)   
        return vectilde
    else:
        raise TypeError("Can only convert a TimeSeries to a FrequencySeries")

def sigmasq(htilde,psd = None,low_frequency_cutoff=None,
            high_frequency_cutoff=None):
    """
    """
    N = (len(htilde)-1) * 2 
    norm = 4.0 / (N * N * htilde.delta_f) 
    kmin,kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff,htilde.delta_f,N)  
 
    moment = htilde.squared_norm()   
    
    if psd is not None:
        moment[kmin:kmax] /= psd[kmin:kmax] 
        
    sq = moment.sum() 
        
    return sq * norm
    
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
    
# Workspace Memory for the matchedfilter
_q = None
_qtilde = None

def matchedfilter(template,data,psd=None,low_frequency_cutoff=None,
                  high_frequency_cutoff=None,calculate_norm=True):
    """Return the complex SNR and normalization (SNR,norm) of the template 
       filtered against the data, where the normalized SNR is SNR' = SNR * norm.
    """
    global _q
    global _qtilde
  
    htilde = get_frequencyseries(template)
    stilde = get_frequencyseries(data)

    N = (len(htilde)-1) * 2   
    kmin,kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff,stilde.delta_f,N) 
   
    if (_q is None) or (len(_q) != N) or _q.dtype != data.dtype:
        _q = zeros(N,dtype=complex_same_precision_as(data))
    else:
        pass
        
    if (_qtilde is None) or (len(_q) != N) or _qtilde.dtype != data.dtype:
        _qtilde = zeros(N,dtype=complex_same_precision_as(data))
    else:
        _q.fill(0)        
    
    #REPLACE with in place operations once they are fixed in PyCUDA
    _qtilde[kmin:kmax] = htilde[kmin:kmax].conj() * stilde[kmin:kmax]

    # Only weight by the psd if it was explictly given.
    # In most cases, the expectation is to prewhiten the data 
    if psd is not None:
        _qtilde[kmin:kmax] /= psd[kmin:kmax]

    ifft(_qtilde,_q)
    
    # Only calculate the normalization if needed. For SPA waveforms this can be 
    # this can be done ahead of time.
    if calculate_norm:
        s_norm = sigmasq(htilde,psd,low_frequency_cutoff,high_frequency_cutoff)
        norm = (4.0 / (N * N * stilde.delta_f)) / sqrt( s_norm) 
    else:
        norm = None
        
    return _q,norm
    
    
def match(vec1,vec2,psd=None,low_frequency_cutoff=None,high_frequency_cutoff=None):
    """ Return the match between the two TimeSeries or FrequencySeries.
    """
    htilde = get_frequencyseries(vec1)
    stilde = get_frequencyseries(vec2)
    snr,norm = matchedfilter(htilde,stilde,psd,low_frequency_cutoff,
                             high_frequency_cutoff)
    maxsnrsq = (snr.squared_norm()).max()
    vec1_normsq = sigmasq(stilde,psd,low_frequency_cutoff,high_frequency_cutoff)
    return sqrt(maxsnrsq/vec1_normsq)*norm
    




    
