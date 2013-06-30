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
import pycbc
from pyfft.cl import Plan
import numpy
from pycbc.types import zeros

def near_two(length):
    n = numpy.log2(length)
    n = numpy.ceil(n)
    return 2**n

_plans = {}


#itype and otype are actual dtypes here, not strings
def _get_plan(itype,otype,inlen):
    try:
        theplan = _plans[(itype,otype,inlen)]
    except KeyError:
        theplan = Plan(inlen,dtype=itype,queue=pycbc.scheme.mgr.state.queue,normalize=False,fast_math=True)
        _plans.update({(itype,otype,inlen) : theplan })

    return theplan

def fft(invec,outvec,prec,itype,otype):
# This has been hacked horrible to make it work more generally
    N = int(near_two(len(outvec)))
    N2 = int (near_two(len(invec)))
    if N2 > N:
        N = N2

    i = zeros(N, dtype = outvec.dtype)
    o = zeros(N, dtype = outvec.dtype)    
    
    i[0:len(invec)] = invec

    pyplan=_get_plan(i.dtype, o.dtype, N)
    pyplan.execute(i.data.data, o.data.data)
    
    outvec.data = o[0:len(outvec)]

def ifft(invec,outvec,prec,itype,otype):
    if itype =='complex' and otype == 'complex':
        pyplan=_get_plan(invec.dtype,outvec.dtype,len(invec))
        pyplan.execute(invec.data.data,outvec.data.data,inverse=True)

    elif itype=='complex' and otype=='real':
        raise NotImplementedError("Only Complex to Complex IFFTs for pyfft currently.")

