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

import numpy.fft

def fft(invec,outvec,prec,itype,otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("numpy backend of pycbc.fft does not support in-place transforms")
    if itype == 'complex' and otype == 'complex':
        outvec.data = numpy.asarray(numpy.fft.fft(invec.data),dtype=outvec.dtype)
    elif itype == 'real' and otype == 'complex':
        outvec.data = numpy.asarray(numpy.fft.rfft(invec.data),dtype=outvec.dtype)

def ifft(invec,outvec,prec,itype,otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("numpy backend of pycbc.fft does not support in-place transforms")
    if itype == 'complex' and otype == 'complex':
        outvec.data = numpy.asarray(numpy.fft.ifft(invec.data),dtype=outvec.dtype)
        outvec *= len(outvec)
    elif itype == 'complex' and otype == 'real':
        outvec.data = numpy.asarray(numpy.fft.irfft(invec.data,len(outvec)),
                                    dtype=outvec.dtype)
        outvec *= len(outvec)

