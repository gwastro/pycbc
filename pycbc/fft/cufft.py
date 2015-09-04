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

import pycbc.scheme
import skcuda.fft as cu_fft
from .core import _BaseFFT, _BaseIFFT

_forward_plans = {}
_reverse_plans = {}

#These dicts need to be cleared before the cuda context is destroyed
def _clear_plan_dicts():
    _forward_plans.clear()
    _reverse_plans.clear()

pycbc.scheme.register_clean_cuda(_clear_plan_dicts)

#itype and otype are actual dtypes here, not strings
def _get_fwd_plan(itype, otype, inlen, batch=1):
    try:
        theplan = _forward_plans[(itype, otype, inlen, batch)]
    except KeyError:
        theplan = cu_fft.Plan((inlen,), itype, otype, batch=batch)
        _forward_plans.update({(itype, otype, inlen) : theplan })

    return theplan

#The complex to real plan wants the actual size, not the N/2+1
#That's why the inverse plans use the outvec length, instead of the invec
def _get_inv_plan(itype, otype, outlen, batch=1):
    try:
        theplan = _reverse_plans[(itype, otype, outlen, batch)]
    except KeyError:
        theplan = cu_fft.Plan((outlen,), itype, otype, batch=batch)
        _reverse_plans.update({(itype, otype, outlen) : theplan })

    return theplan


def fft(invec, outvec, prec, itype, otype):
    cuplan = _get_fwd_plan(invec.dtype, outvec.dtype, len(invec))
    cu_fft.fft(invec.data, outvec.data, cuplan)

def ifft(invec, outvec, prec, itype, otype):
    cuplan = _get_inv_plan(invec.dtype, outvec.dtype, len(outvec))
    cu_fft.ifft(invec.data, outvec.data, cuplan)
    
class FFT(_BaseFFT):
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(FFT, self).__init__(invec, outvec, nbatch, size)
        self.plan = _get_fwd_plan(invec.dtype, outvec.dtype, len(invec), batch=nbatch)
        self.invec = invec.data
        self.outvec = outvec.data

    def execute(self):
        cu_fft.fft(self.invec, self.outvec, self.plan)

class IFFT(_BaseIFFT):
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(IFFT, self).__init__(invec, outvec, nbatch, size)
        self.plan = _get_inv_plan(invec.dtype, outvec.dtype, len(outvec), batch=nbatch)

        self.invec = invec.data
        self.outvec = outvec.data

    def execute(self):
        cu_fft.ifft(self.invec, self.outvec, self.plan)

