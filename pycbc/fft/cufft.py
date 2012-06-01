# Copyright (C) 2012  Josh Willis, Andrew Miller
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
This module provides the cufft backend of the fast Fourier transform
for the PyCBC package.
"""

import pycbc.array
import scikits.cuda.fft as cu_fft

def fft(invec,outvec,prec,itype,otype):
    outvec.data #Move output if necessary
    invec.data #Move input if necessary
    if itype =='complex' and otype == 'complex':
        cuplan = cu_fft.Plan((len(invec),),invec.dtype,outvec.dtype)
        cu_fft.fft(invec.data,outvec.data,cuplan)

    elif itype=='real' and otype=='complex':
        #The cufft algorithm doesn't return exact zeros for imaginary parts of this transform
        #it returns imaginary components on the order of 10^-16. Because of this, the forward
        #real to complex tests do not currently pass the unit tests.
        cuplan = cu_fft.Plan((len(invec),),invec.dtype,outvec.dtype)
        cu_fft.fft(invec.data,outvec.data,cuplan)

def ifft(invec,outvec,prec,itype,otype):
    outvec.data #Move output if necessary
    invec.data #Move input if necessary
    if itype =='complex' and otype == 'complex':
        cuplan = cu_fft.Plan((len(invec),),invec.dtype,outvec.dtype)
        cu_fft.ifft(invec.data,outvec.data,cuplan)

    elif itype=='complex' and otype=='real':
        #This plan wants the actual size, not the N/2+1
        #That's why this plan uses the outvec length, instead of the invec
        cuplan = cu_fft.Plan((len(outvec),),invec.dtype,outvec.dtype)
        cu_fft.ifft(invec.data,outvec.data,cuplan)

