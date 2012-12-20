# Copyright (C) 2012  Alex Nitz
#
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
"""This module contains functions to generate guassian noise colored with a 
noise spectrum. 
"""
from pycbc.types import TimeSeries, zeros
from lalsimulation import SimNoise
import lal

def noise_from_psd(length, delta_t, psd, seed=0):
    """Produce noise colored by a psd.
    """
    noise_ts = TimeSeries(zeros(length), delta_t=delta_t)
    
    randomness = lal.gsl_rng("ranlux", seed)
    
    N = int (1.0 / delta_t / psd.delta_f)
    n = N/2+1
    stride = N/2
    
    if n > len(psd):
        raise ValueError("PSD not compatible with requested delta_t")
        
    psd = (psd[0:n]).lal()
    psd.data.data[n-1] = 0 
    
    segment = TimeSeries(zeros(N), delta_t=delta_t).lal()
    length_generated = 0 
  
    SimNoise(segment, 0, psd, randomness) 
    while (length_generated < length):
        if (length_generated + stride) < length:
            noise_ts.data[length_generated:length_generated+stride]  = segment.data.data[0:stride]
        else: 
            noise_ts.data[length_generated:length]  = segment.data.data[0:length-length_generated]
    
        length_generated += stride
        SimNoise(segment, stride, psd, randomness)
        
    return noise_ts
    
    
