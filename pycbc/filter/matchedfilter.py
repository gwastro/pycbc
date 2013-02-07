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

from math import log,ceil,sqrt
from pycbc.types import TimeSeries,FrequencySeries,zeros, Array
from pycbc.types import complex_same_precision_as,real_same_precision_as
from pycbc.fft import fft,ifft
import pycbc.scheme
import pycbc
import numpy

BACKEND_PREFIX="pycbc.filter.matchedfilter_"

@pycbc.scheme.schemed(BACKEND_PREFIX)
def correlate(x, y, z):
    pass


def make_frequency_series(vec):
    """Returns a FrequencySeries from the input vector
    """
    if isinstance(vec, FrequencySeries):
        return vec
    if isinstance(vec, TimeSeries):
        N = len(vec)
        n = N/2+1    
        delta_f = 1.0 / N / vec.delta_t
        vectilde =  FrequencySeries(zeros(n, dtype=complex_same_precision_as(vec)), delta_f=delta_f, copy=False)
        fft(vec, vectilde)   
        return vectilde
    else:
        raise TypeError("Can only convert a TimeSeries to a FrequencySeries")

def sigmasq_series(htilde, psd=None, low_frequency_cutoff=None,
            high_frequency_cutoff=None):
    """Returns the accumulated loudness of the template over frequency bins. 
    """
    htilde = make_frequency_series(htilde)
    N = (len(htilde)-1) * 2 
    norm = 4.0 * htilde.delta_f
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, htilde.delta_f, N)  
   
    sigma_vec = FrequencySeries(zeros(len(htilde)), delta_f = htilde.delta_f, 
                        dtype=real_same_precision_as(htilde), copy=False)
    
    mag = htilde.squared_norm()

    sigma_vec[kmin:kmax] = mag[kmin:kmax].cumsum()
        
    return sigma_vec*norm
    

def sigmasq(htilde, psd = None, low_frequency_cutoff=None,
            high_frequency_cutoff=None):
    """Returns the loudness squared of the template.
    """
    htilde = make_frequency_series(htilde)
    N = (len(htilde)-1) * 2 
    norm = 4.0 * htilde.delta_f
    kmin, kmax = get_cutoff_indices(low_frequency_cutoff,
                                   high_frequency_cutoff, htilde.delta_f, N)  
    ht = htilde[kmin:kmax] 

    if psd is None:
        sq = ht.inner(ht)
    else:
        sq = ht.weighted_inner(ht, psd[kmin:kmax])
        
    return sq.real * norm

def sigma(htilde, psd = None, low_frequency_cutoff=None,
        high_frequency_cutoff=None):
    """Returns the loudness of the template.
    """
    return sqrt(sigmasq(htilde, psd, low_frequency_cutoff, high_frequency_cutoff))
    

def get_cutoff_indices(flow, fhigh, df, N):
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
_qtilde = None

def matched_filter(template, data, psd=None, low_frequency_cutoff=None,
                  high_frequency_cutoff=None, h_norm=None, out=None):
    """Return the complex SNR and normalization 
    """
    global _qtilde
  
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
        
    if (_qtilde is None) or (len(_qtilde) != N) or _qtilde.dtype != data.dtype:
        _qtilde = zeros(N, dtype=complex_same_precision_as(data))
    else:
        _qtilde.clear()         
    
    correlate(htilde[kmin:kmax], stilde[kmin:kmax], _qtilde[kmin:kmax])

    # Only weight by the psd if it was explictly given.
    # In most cases, the expectation is to overwhiten the data 
    if psd is not None:
        if isinstance(psd, FrequencySeries):
            if psd.delta_f == stilde.delta_f :
                _qtilde[kmin:kmax] /= psd[kmin:kmax]
            else:
                raise TypeError("PSD delta_f does not match data")
        else:
            raise TypeError("PSD must be a FrequencySeries")
            
    ifft(_qtilde,_q)
    
    # Only calculate the normalization if needed. For SPA waveforms
    # this can be done ahead of time.
    if h_norm is None:
        h_norm = sigmasq(htilde, psd, low_frequency_cutoff, high_frequency_cutoff)     

    norm = (4.0 * stilde.delta_f) / sqrt( h_norm)
    delta_t = 1.0 / (N * stilde.delta_f)
    return TimeSeries(_q, epoch=stilde._epoch, delta_t=delta_t, copy=False), norm
    
_snr = None 
def match(vec1, vec2, psd=None, low_frequency_cutoff=None,
          high_frequency_cutoff=None, s_norm=None, h_norm=None):
    """ Return the match between the two TimeSeries or FrequencySeries.
    """

    htilde = make_frequency_series(vec1)
    stilde = make_frequency_series(vec2)

    N = (len(htilde)-1) * 2

    global _snr
    if _snr is None or _snr.dtype != htilde.dtype or len(_snr) != N:
        _snr = zeros(N,dtype=complex_same_precision_as(vec1))
    snr, snr_norm = matched_filter(htilde,stilde,psd,low_frequency_cutoff,
                             high_frequency_cutoff, h_norm, out=_snr)
    maxsnr, max_id = snr.abs_max_loc()
    if s_norm is None:
        s_norm = sigmasq(stilde, psd, low_frequency_cutoff, high_frequency_cutoff)
    return maxsnr * snr_norm / sqrt(s_norm), max_id
    
def overlap(vec1, vec2, psd=None, low_frequency_cutoff=None,
          high_frequency_cutoff=None, normalized=True):
    """ Return the overlap between the two TimeSeries or FrequencySeries.
    """
    htilde = make_frequency_series(vec1)
    stilde = make_frequency_series(vec2)

    if psd:
        inner = htilde.weighted_inner(stilde, psd).real
    else:
        inner = htilde.inner(stilde).real
     
    if normalized:   
        sig1 = sigma(vec1, psd=psd, low_frequency_cutoff=low_frequency_cutoff, high_frequency_cutoff=high_frequency_cutoff)
        sig2 = sigma(vec2, psd=psd, low_frequency_cutoff=low_frequency_cutoff, high_frequency_cutoff=high_frequency_cutoff)
        norm = 1 / sig1 / sig2
    else:
        norm = 1
        
    return 4 * htilde.delta_f * inner * norm
      

__all__ = ['match', 'matched_filter', 'sigmasq', 'sigma', 'sigmasq_series', 'make_frequency_series', 'overlap']

