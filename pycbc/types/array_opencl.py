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
"""PyOpenCL based backend for pycbc array
"""
from pyopencl.tools import context_dependent_memoize
from pyopencl.scan import InclusiveScanKernel
from pyopencl.reduction import ReductionKernel
from pytools import match_precision, memoize_method
from pyopencl.array import _get_common_dtype
import pyopencl.array
from pyopencl.elementwise import ElementwiseKernel, complex_dtype_to_name
from pyopencl.tools import dtype_to_ctype
from pycbc.scheme import mgr 
import numpy as np

complex_headers = """
#define PYOPENCL_DEFINE_CDOUBLE
#include <pyopencl-complex.h>
"""

@context_dependent_memoize
def get_cumsum_kernel(dtype):
    return InclusiveScanKernel(mgr.state.context, dtype, "a+b", 
                                          neutral="0", preamble=complex_headers)

def icumsum(vec):
    krnl = get_cumsum_kernel(vec.dtype)
    return krnl(vec)

 
@context_dependent_memoize
def get_weighted_inner_kernel(dtype_x, dtype_y, dtype_w, dtype_out):
    if (dtype_x == np.complex64) or (dtype_x == np.complex128):
        if (dtype_y == np.float64) or (dtype_y == np.float32):
            ys = "%s_fromreal(y[i])" % complex_dtype_to_name(dtype_x)
        else:
            ys = "y[i]"
        inner_map="%s_mul(%s_conj(x[i]), %s)" % (complex_dtype_to_name(dtype_x), complex_dtype_to_name(dtype_x), ys)
    else:
        inner_map="x[i]*y[i]"  
        
    if (dtype_w == np.float64) or (dtype_w == np.float32):
        inner_map = inner_map + "/w[i]"  
    else:
        inner_map = "%s_divide(%s, %s)" % (complex_dtype_to_name(dtype_x), inner_map, "w[i]")
               
    return ReductionKernel(mgr.state.context, dtype_out,
            neutral="0",
            arguments="__global const %(tp_x)s *x, __global const %(tp_y)s *y, __global const %(tp_w)s *w" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_y": dtype_to_ctype(dtype_y),
                "tp_w": dtype_to_ctype(dtype_w),
                },
            reduce_expr="a+b",
            map_expr=inner_map,
            name="weighted_inner")

@context_dependent_memoize
def get_inner_kernel(dtype_x, dtype_y, dtype_out):
    if (dtype_x == np.complex64) or (dtype_x == np.complex128):
        if (dtype_y == np.float64) or (dtype_y == np.float32):
            ys = "%s_fromreal(y[i])" % complex_dtype_to_name(dtype_x)
        else:
            ys = "y[i]"
        inner_map="%s_mul(%s_conj(x[i]), %s)" % (complex_dtype_to_name(dtype_x), complex_dtype_to_name(dtype_x), ys)
    else:
        inner_map="x[i]*y[i]"    
               
    return ReductionKernel(mgr.state.context, dtype_out,
            neutral="0",
            arguments="__global const %(tp_x)s *x, __global const %(tp_y)s *y" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_y": dtype_to_ctype(dtype_y),
                },
            reduce_expr="a+b",
            map_expr=inner_map,
            name="inner")

def inner(a, b):
    dtype_out = _get_common_dtype(a, b, mgr.state.queue)
    krnl = get_inner_kernel(a.dtype, b.dtype, dtype_out)
    return krnl(a, b)

def weighted_inner(self, b, w):
    if weight is None:
        return self.inner(other) 
    a = self.data
    dtype_out = _get_common_dtype(a, b, mgr.state.queue)
    krnl = get_weighted_inner_kernel(a.dtype, b.dtype, w.dtype, dtype_out)
    return krnl(a, b, w).get().max()
 

@context_dependent_memoize
def get_norm_kernel(dtype_x, dtype_out):
    return ElementwiseKernel(mgr.state.context, 
            "%(tp_x)s *x, %(tp_z)s *z" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_z": dtype_to_ctype(dtype_out),
                },
            "z[i] = norm(x[i])",
            "norm")

def squared_norm(a):
#    dtype_out = match_precision(np.dtype('float64'), a.dtype)
#    out = a._new_like_me(dtype=dtype_out)
#    krnl = get_norm_kernel(a.dtype,dtype_out)
#    krnl(a,out)
#    return out 
    return a.conj() * a

def zeros(length, dtype=np.float64):
    return pyopencl.array.zeros(mgr.state.queue, length, dtype)

def ptr(self):
    return self._data.data

def dot(self, other):
    return pyopencl.array.dot(self._data,other).get().max()
    
def min(self):
    return pyopencl.array.min(self._data).get().max()  

def cumsum(self):
    tmp = self.data*1
    return icumsum(tmp)

def max(self):
    return pyopencl.array.max(self._data).get().max()
                   


