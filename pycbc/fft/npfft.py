# Copyright (C) 2012  Josh Willis
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
This module provides the numpy backend of the fast Fourier transform
for the PyCBC package.
"""

import pycbc.array
import numpy.fft

def fft(invec,outvec,prec,itype,otype):
    if itype is 'complex' and otype is 'complex':
        outvec.data() # Just to move output, if necessary
        outvec._data = numpy.fft.fft(invec.data())
    elif itype is 'real' and otype is 'complex':
        outvec.data() # Just to move output, if necessary
        outvec._data = numpy.fft.rfft(invec.data())

def ifft(invec,outvec,backend=None):
    if itype is 'complex' and otype is 'complex':
        outvec.data() # Just to move output, if necessary
        outvec._data = numpy.fft.ifft(invec.data())
        outvec *= len(outvec)
    elif itype is 'complex' and otype is 'real':
        outvec.data() # Just to move output, if necessary
        outvec._data = numpy.fft.irfft(invec.data())
        outvec *= len(outvec)

