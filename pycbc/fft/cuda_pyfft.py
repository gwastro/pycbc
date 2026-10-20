# Copyright (C) 2012  Andrew Miller
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
This module provides the pyfft backend of the fast Fourier transform
for the PyCBC package.
"""

import pycbc.scheme
from pyfft.cuda import Plan

_plans = {}

#These dicts need to be cleared before the cuda context is destroyed
def _clear_plan_dict():
    _plans.clear()

pycbc.scheme.register_clean_cuda(_clear_plan_dict)


#itype and otype are actual dtypes here, not strings
def _get_plan(itype,otype,inlen):
    try:
        theplan = _plans[(itype,otype,inlen)]
    except KeyError:
        theplan = Plan(inlen,dtype = itype,normalize=False,fast_math=True)
        _plans.update({(itype,otype,inlen) : theplan })

    return theplan

def fft(invec,outvec,prec,itype,otype):
    if itype =='complex' and otype == 'complex':
        pyplan=_get_plan(invec.dtype, outvec.dtype, len(invec))
        pyplan.execute(invec.data,outvec.data)

    elif itype=='real' and otype=='complex':
        raise NotImplementedError("Only Complex to Complex FFTs for pyfft currently.")

def ifft(invec,outvec,prec,itype,otype):
    if itype =='complex' and otype == 'complex':
        pyplan=_get_plan(invec.dtype,outvec.dtype,len(invec))
        pyplan.execute(invec.data,outvec.data,inverse=True)

    elif itype=='complex' and otype=='real':
        raise NotImplementedError("Only Complex to Complex IFFTs for pyfft currently.")

