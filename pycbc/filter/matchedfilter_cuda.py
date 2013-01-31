# Copyright (C) 2012  Alex Nitz
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

from pycuda.elementwise import get_elwise_kernel,ElementwiseKernel
from pytools import match_precision
from pycuda.reduction import ReductionKernel
from pycuda.tools import context_dependent_memoize
from pycuda.tools import dtype_to_ctype
from pycbc.types import real_same_precision_as
from pycuda.gpuarray import _get_common_dtype
import numpy

@context_dependent_memoize
def get_correlate_kernel(dtype_x, dtype_y,dtype_out):
    return ElementwiseKernel(
            "%(tp_x)s *x, %(tp_y)s *y, %(tp_z)s *z" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_y": dtype_to_ctype(dtype_y),
                "tp_z": dtype_to_ctype(dtype_out),
                },
            "z[i] = conj(x[i]) * y[i]",
            "correlate")

def correlate(a, b, out, stream=None):
    dtype_out = _get_common_dtype(a,b)
    krnl = get_correlate_kernel(a.dtype, b.dtype, dtype_out)
    krnl(a.data,b.data,out.data)




