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
This module provides the SwigLAL backend of the fast Fourier transform
for the PyCBC package.
"""

import pycbc.types
import swiglal
from swiglal import XLALCreateForwardCOMPLEX16FFTPlan as _CreateForwardCOMPLEX16FFTPlan
from swiglal import XLALCreateForwardCOMPLEX8FFTPlan as _CreateForwardCOMPLEX8FFTPlan
from swiglal import XLALCreateForwardREAL8FFTPlan as _CreateForwardREAL8FFTPlan
from swiglal import XLALCreateForwardREAL4FFTPlan as _CreateForwardREAL4FFTPlan
from swiglal import XLALCreateReverseCOMPLEX16FFTPlan as _CreateReverseCOMPLEX16FFTPlan
from swiglal import XLALCreateReverseCOMPLEX8FFTPlan as _CreateReverseCOMPLEX8FFTPlan
from swiglal import XLALCreateReverseREAL8FFTPlan as _CreateReverseREAL8FFTPlan
from swiglal import XLALCreateReverseREAL4FFTPlan as _CreateReverseREAL4FFTPlan

from swiglal import XLALCOMPLEX16VectorFFT as _COMPLEX16VectorFFT
from swiglal import XLALCOMPLEX8VectorFFT as _COMPLEX8VectorFFT
from swiglal import XLALREAL8ForwardFFT as _REAL8ForwardFFT
from swiglal import XLALREAL4ForwardFFT as _REAL4ForwardFFT
from swiglal import XLALCOMPLEX16VectorFFT as _COMPLEX16VectorFFT
from swiglal import XLALCOMPLEX8VectorFFT as _COMPLEX8VectorFFT
from swiglal import XLALREAL8ReverseFFT as _REAL8ReverseFFT
from swiglal import XLALREAL4ReverseFFT as _REAL4ReverseFFT

_default_measurelvl = 1
_forward_plans = {}
_reverse_plans = {}
_forward_plan_fn_dict = {('double','complex','complex') : _CreateForwardCOMPLEX16FFTPlan,
                    ('single','complex','complex') : _CreateForwardCOMPLEX8FFTPlan,
                    ('double','real','complex') : _CreateForwardREAL8FFTPlan,
                    ('single','real','complex') : _CreateForwardREAL4FFTPlan}
_reverse_plan_fn_dict = {('double','complex','complex') : _CreateReverseCOMPLEX16FFTPlan,
                    ('single','complex','complex') : _CreateReverseCOMPLEX8FFTPlan,
                    ('double','complex','real') : _CreateReverseREAL8FFTPlan,
                    ('single','complex','real') : _CreateReverseREAL4FFTPlan}
_forward_fft_fn_dict = {('double','complex','complex') : _COMPLEX16VectorFFT,
                    ('single','complex','complex') : _COMPLEX8VectorFFT,
                    ('double','real','complex') : _REAL8ForwardFFT,
                    ('single','real','complex') : _REAL4ForwardFFT}
_reverse_fft_fn_dict = {('double','complex','complex') : _COMPLEX16VectorFFT,
                    ('single','complex','complex') : _COMPLEX8VectorFFT,
                    ('double','complex','real') : _REAL8ReverseFFT,
                    ('single','complex','real') : _REAL4ReverseFFT}

def _get_fwd_plan(prec,itype,otype,inlen):
    try:
        theplan = _forward_plans[(prec,itype,otype,inlen)]
    except KeyError:
        theplan = _forward_plan_fn_dict[(prec,itype,otype)](inlen,_default_measurelvl)
        _forward_plans.update({(prec,itype,otype,inlen) : theplan})
    return theplan

def _get_inv_plan(prec,itype,otype,outlen):
    try:
        theplan = _reverse_plans[(prec,itype,otype,outlen)]
    except KeyError:
        theplan = _reverse_plan_fn_dict[(prec,itype,otype)](outlen,_default_measurelvl)
        _reverse_plans.update({(prec,itype,otype,outlen) : theplan})
    return theplan


def fft(invec,outvec,prec,itype,otype):
    theplan = _get_fwd_plan(prec,itype,otype,len(invec))
    _forward_fft_fn_dict[(prec,itype,otype)](outvec.lal,invec.lal,theplan)

def ifft(invec,outvec,prec,itype,otype):
    theplan = _get_inv_plan(prec,itype,otype,len(outvec))
    _reverse_fft_fn_dict[(prec,itype,otype)](outvec.lal,invec.lal,theplan)
