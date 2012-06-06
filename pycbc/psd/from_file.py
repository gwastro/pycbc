# Copyright (C) 2012  Alex Nitz
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
"""
Helper functions to load psds and asds from ascii files.
"""

import numpy
from pycbc import DYN_RANGE_FAC
from pycbc.types import FrequencySeries
from scipy.interpolate import interp1d
from pycbc.filter import *

def psd_from_file(filename,delta_f,length,low_frequency_cutoff=None,high_frequency_cutoff=None,dtype=None):
    """Returns the psd from an ascii file containing an asd
    """
    fpsd = numpy.loadtxt(filename)          
    freq_data=fpsd[:,0]
    psd_data=fpsd[:,1]*DYN_RANGE_FAC
    psd_data = psd_data **2
    psd_interp= interp1d(freq_data,psd_data) 
    N = (length-1) * 2  

    kmin,kmax = get_cutoff_indices(low_frequency_cutoff,high_frequency_cutoff,delta_f,N) 

    psd = numpy.zeros(length)
    for k in range(0,kmax,1):
        if (k<kmin):
            psd[k]=0
        else:
            psd[k]=float(psd_interp( k* delta_f ) )                   

   
    return FrequencySeries(psd,delta_f=delta_f,dtype=dtype,copy=False)
  
#def asd_from_file(filename,delta_f,length):
#    """Returns the interpolated asd from an ascii file containing an asd 
#    """
#    asd = psd_from_file(filename,delta_f,length)
#    return asd ** 2
    
