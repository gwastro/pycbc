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

import pycbc.array
import swiglal

_default_measurelvl = 1
_forward_plans = {}
_reverse_plans = {}
_forward_fn_dict = {('double','complex','complex') : swiglal.XLALCreateForwardCOMPLEX16FFTPlan,
                    ('single','complex','complex') : swiglal.XLALCreateForwardCOMPLEX8FFTPlan,
                    ('double','real','complex') : swiglal.XLALCreateForwardREAL8FFTPlan,
                    ('single','real','complex') : swiglal.XLALCreateForwardREAL4FFTPlan}
_reverse_fn_dict = {('double','complex','complex') : swiglal.XLALCreateReverseCOMPLEX16FFTPlan,
                    ('single','complex','complex') : swiglal.XLALCreateReverseCOMPLEX8FFTPlan,
                    ('double','complex','real') : swiglal.XLALCreateReverseREAL8FFTPlan,
                    ('single','complex','real') : swiglal.XLALCreateReverseREAL4FFTPlan}

def _get_fwd_plan(prec,itype,otype,inlen):
    try:
        theplan = _forward_plans[(prec,itype,otype,inlen)]
    except KeyError:
        theplan = _forward_fn_dict[(prec,itype,otype)](inlen,_default_measurelvl)
        _forward_plans.update({(prec,itype,otype,inlen) : theplan})
    return theplan

def _get_inv_plan(prec,itype,otype,outlen):
    try:
        theplan = _reverse_plans[(prec,itype,otype,outlen)]
    except KeyError:
        theplan = _reverse_fn_dict[(prec,itype,otype)](outlen,_default_measurelvl)
        _reverse_plans.update({(prec,itype,otype,outlen) : theplan})
    return theplan


def fft(invec,outvec,prec,itype,otype):
    theplan = _get_fwd_plan(prec,itype,otype,len(invec))
    if itype is 'complex' and otype is 'complex':
        if prec is 'single':
            swiglal.XLALCOMPLEX8VectorFFT(outvec.lal,invec.lal,theplan)
        else:
            swiglal.XLALCOMPLEX16VectorFFT(outvec.lal,invec.lal,theplan)
    elif itype is 'real' and otype is 'complex':
        if prec is 'single':
            swiglal.XLALREAL4ForwardFFT(outvec.lal,invec.lal,theplan)
        else:
            swiglal.XLALREAL8ForwardFFT(outvec.lal,invec.lal,theplan)

def ifft(invec,outvec,prec,itype,otype):
    theplan = _get_inv_plan(prec,itype,otype,len(outvec))
    if itype is 'complex' and otype is 'complex':
        if prec is 'single':
            swiglal.XLALCOMPLEX8VectorFFT(outvec.lal,invec.lal,theplan)
        else:
            swiglal.XLALCOMPLEX16VectorFFT(outvec.lal,invec.lal,theplan)
    elif itype is 'complex' and otype is 'real':
        if prec is 'single':
            swiglal.XLALREAL4ReverseFFT(outvec.lal,invec.lal,theplan)
        else:
            swiglal.XLALREAL8ReverseFFT(outvec.lal,invec.lal,theplan)

