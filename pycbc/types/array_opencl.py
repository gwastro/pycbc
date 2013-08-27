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
from pyopencl.array import vdot as pvdot
import pyopencl.array
from pyopencl.elementwise import ElementwiseKernel, complex_dtype_to_name
from pyopencl.tools import dtype_to_ctype, get_or_register_dtype
from pycbc.scheme import mgr 
import numpy as np

complex_headers = """
#define PYOPENCL_DEFINE_CDOUBLE
#include <pyopencl-complex.h>
"""

# Define pyopencl MAXLOC for both single and double precission ################## 
       
maxloc_preamble = """

    typedef struct {
        TTYPE max;
        LTYPE   loc;
    }MAXLOCN;
    
    __global
    MAXLOCN maxloc_red(MAXLOCN a, MAXLOCN b){
        if (a.max > b.max)
            return a;
        else 
            return b;  
    }
    
    __global
    MAXLOCN maxloc_start(){
        MAXLOCN t;
        t.max=0;
        t.loc=0;
        return t;
    }
    
    __global
    MAXLOCN maxloc_map(TTYPE val, LTYPE loc){
        MAXLOCN t;
        t.max = val;
        t.loc = loc;
        return t;
    }

    """
    
maxloc_preamble_single = """
    #define MAXLOCN maxlocs
    #define TTYPE float
    #define LTYPE int
""" + maxloc_preamble

maxloc_preamble_double = """
    #define MAXLOCN maxlocd
    #define TTYPE double
    #define LTYPE long
""" + maxloc_preamble
    
maxloc_dtype_double = np.dtype([("max", np.float64), ("loc", np.int64)])
maxloc_dtype_single = np.dtype([("max", np.float32), ("loc", np.int32)])

maxloc_dtype_single = get_or_register_dtype("maxlocs", dtype=maxloc_dtype_single)
maxloc_dtype_double = get_or_register_dtype("maxlocd", dtype=maxloc_dtype_double)

mls = ReductionKernel(mgr.state.context, maxloc_dtype_single, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(x[i], i)",
        arguments="float *x", preamble=maxloc_preamble_single)

mld = ReductionKernel(mgr.state.context, maxloc_dtype_double, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(x[i], i)",
        arguments="double *x", preamble=maxloc_preamble_double)
        
max_loc_map = {'single':mls,'double':mld}

cfloat  = complex_dtype_to_name(np.complex64)
cdouble = complex_dtype_to_name(np.complex128)

get_or_register_dtype('cfloat', np.complex64)
get_or_register_dtype('cdouble', np.complex128)

amls = ReductionKernel(mgr.state.context, maxloc_dtype_single, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(fabs(x[i]), i)",
        arguments="float *x", preamble=maxloc_preamble_single)

amld = ReductionKernel(mgr.state.context, maxloc_dtype_double, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(fabs(x[i]), i)",
        arguments="double *x", preamble=maxloc_preamble_double)

amlsc = ReductionKernel(mgr.state.context, maxloc_dtype_single, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(%s_abs(x[i]), i)" % cfloat,
       arguments="%s *x" % cfloat, preamble=maxloc_preamble_single)

amldc = ReductionKernel(mgr.state.context, maxloc_dtype_double, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(%s_abs(x[i]), i)" % cdouble,
        arguments="%s *x" % cdouble, preamble=maxloc_preamble_double)

abs_max_loc_map = {'single':{ 'real':amls, 'complex':amlsc }, 'double':{ 'real':amld, 'complex':amldc }}


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

def inner(self, b):
   # a = self.data
   # dtype_out = _get_common_dtype(a, b, mgr.state.queue)
   # krnl = get_inner_kernel(a.dtype, b.dtype, dtype_out)
   # return krnl(a, b).get().max()
   return pvdot(self.data, b).get().max()

vdot = inner   

def weighted_inner(self, b, w):
  #  if w is None:
  #      return self.inner(other) 
  #  a = self.data
  #  dtype_out = _get_common_dtype(a, b, mgr.state.queue)
  #  krnl = get_weighted_inner_kernel(a.dtype, b.dtype, w.dtype, dtype_out)
  #  return krnl(a, b, w).get().max()
    return pvdot(self.data, b/w).get().max()
 

@context_dependent_memoize
def get_norm_kernel(dtype_x, dtype_out):
    if dtype_x == np.float32 or dtype_x == np.float64:
        op = "z[i] = x[i] * x[i]"
    if dtype_x == np.complex64 or dtype_x == np.complex128:
        op = "z[i] = x[i].x*x[i].x + x[i].y*x[i].y"
    return ElementwiseKernel(mgr.state.context, 
            "%(tp_x)s *x, %(tp_z)s *z" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_z": dtype_to_ctype(dtype_out),
                },
            op,
            "normsq")

def squared_norm(self):
    dtype_out = match_precision(np.dtype('float64'), self.dtype)
    out = self.data._new_like_me(dtype=dtype_out)
    krnl = get_norm_kernel(self.dtype, dtype_out)
    krnl(self.data, out)
    return out 

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
                   
def numpy(self):
    return self._data.get()     
 
def take(self, indices):
    indices = pyopencl.array.to_device(mgr.state.queue, indices)
    return pyopencl.array.take(self.data, indices)
    
def abs_max_loc(self):
    maxloc = abs_max_loc_map[self.precision][self.kind](self._data)
    maxloc = maxloc.get()
    return float(maxloc['max']),int(maxloc['loc'])

def max_loc(self):
    maxloc = max_loc_map[self.precision](self._data)
    maxloc = maxloc.get()
    return float(maxloc['max']),int(maxloc['loc'])

def _copy(self, self_ref, other_ref):
    if (len(other_ref) <= len(self_ref)) :
        self_ref._copy(self_ref, other_ref)
    else:
        raise RuntimeError("The arrays must the same length")

def _getvalue(self, index):
    return self._data.get()[index]
    
def sum(self):
    return pyopencl.array.sum(self._data).get().max() 

def clear(self):
    self.fill(0)
    
def _scheme_matches_base_array(array):
    if isinstance(array, pyopencl.array.Array):
        return True
    else:
        return False
        
def _copy_base_array(array):
    #Unsure how to copy memory directly with OpenCL, these two lines are 
    #a temporary workaround, still probably better than going to the host and back though
    # This function doesn't behave nicely when the array is empty
    # because it tries to fill it with zeros
    if len(array) > 0:
        data = pyopencl.array.zeros_like(array)
        data += array
    else:
        data = pyopencl.array.Array(mgr.state.queue,(0,), array.dtype)
    return data
    
def _to_device(array):
    return pyopencl.array.to_device(mgr.state.queue, array)
