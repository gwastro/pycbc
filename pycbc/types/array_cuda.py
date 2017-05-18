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
"""Pycuda based 
"""
import pycuda.driver
from pycuda.elementwise import ElementwiseKernel
from pycuda.reduction import ReductionKernel
from pycuda.tools import get_or_register_dtype
from pycuda.tools import context_dependent_memoize
from pycuda.tools import dtype_to_ctype
from pytools import match_precision
from pycuda.gpuarray import _get_common_dtype, empty, GPUArray
import pycuda.gpuarray
from pycuda.scan import InclusiveScanKernel
import numpy as np

include_complex = """
#include <pycuda-complex.hpp>
"""

@context_dependent_memoize
def get_cumsum_kernel(dtype):
    return InclusiveScanKernel(dtype, "a+b", preamble=include_complex)

def icumsum(vec):
    krnl = get_cumsum_kernel(vec.dtype)
    return krnl(vec)

@context_dependent_memoize
def call_prepare(self, sz, allocator):
    MAX_BLOCK_COUNT = 1024
    SMALL_SEQ_COUNT = 4        

    if sz <= self.block_size*SMALL_SEQ_COUNT*MAX_BLOCK_COUNT:
        total_block_size = SMALL_SEQ_COUNT*self.block_size
        block_count = (sz + total_block_size - 1) // total_block_size
        seq_count = SMALL_SEQ_COUNT
    else:
        block_count = MAX_BLOCK_COUNT
        macroblock_size = block_count*self.block_size
        seq_count = (sz + macroblock_size - 1) // macroblock_size

    if block_count == 1:
        result = empty((), self.dtype_out, allocator)
    else:
        result = empty((block_count,), self.dtype_out, allocator)

    grid_size = (block_count, 1)
    block_size =  (self.block_size, 1, 1)

    return result, block_count, seq_count, grid_size, block_size

class LowerLatencyReductionKernel(ReductionKernel):
    def __init__(self, dtype_out,
            neutral, reduce_expr, map_expr=None, arguments=None,
            name="reduce_kernel", keep=False, options=None, preamble=""):
            ReductionKernel.__init__(self, dtype_out,
                neutral, reduce_expr, map_expr, arguments,
                name, keep, options, preamble)

            self.shared_size=self.block_size*self.dtype_out.itemsize


    def __call__(self, *args, **kwargs):
        f = self.stage1_func
        s1_invocation_args = [] 
        for arg in args:
            s1_invocation_args.append(arg.gpudata)
        sz = args[0].size

        result, block_count, seq_count, grid_size, block_size = call_prepare(self, sz, args[0].allocator)

        f(grid_size, block_size, None,
                *([result.gpudata]+s1_invocation_args+[seq_count, sz]),
                shared_size=self.shared_size)

        while True:
            f = self.stage2_func
            sz = result.size
            result2 = result
            result, block_count, seq_count, grid_size, block_size = call_prepare(self, sz, args[0].allocator)

            f(grid_size, block_size, None,
                    *([result.gpudata, result2.gpudata]+s1_invocation_args+[seq_count, sz]),
                    shared_size=self.shared_size)

            if block_count == 1:
                return result



@context_dependent_memoize
def get_norm_kernel(dtype_x, dtype_out):
    return ElementwiseKernel(
            "%(tp_x)s *x, %(tp_z)s *z" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_z": dtype_to_ctype(dtype_out),
                },
            "z[i] = norm(x[i])",
            "normalize")

def squared_norm(self):
    a = self.data
    dtype_out = match_precision(np.dtype('float64'), a.dtype)
    out = a._new_like_me(dtype=dtype_out)
    krnl = get_norm_kernel(a.dtype, dtype_out)
    krnl(a, out)
    return out     

# FIXME: Write me!
#def multiply_and_add(self, other, mult_fac):
#    """
#    Return other multiplied by mult_fac and with self added.
#    Self will be modified in place. This requires all inputs to be of the same
#    precision.
#    """
 
@context_dependent_memoize
def get_weighted_inner_kernel(dtype_x, dtype_y, dtype_w, dtype_out):
    if (dtype_x == np.complex64) or (dtype_x == np.complex128):
        inner_map="conj(x[i])*y[i]/w[i]"
    else:
        inner_map="x[i]*y[i]/w[i]"       
    return LowerLatencyReductionKernel(dtype_out,
            neutral="0",
            arguments="%(tp_x)s *x, %(tp_y)s *y,  %(tp_w)s *w" % {
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
        inner_map="conj(x[i])*y[i]"
    else:
        inner_map="x[i]*y[i]"           
    return LowerLatencyReductionKernel(dtype_out,
            neutral="0",
            arguments="%(tp_x)s *x, %(tp_y)s *y" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_y": dtype_to_ctype(dtype_y),
                },
            reduce_expr="a+b",
            map_expr=inner_map,
            name="inner")

def inner(self, b):
    a = self.data
    dtype_out = _get_common_dtype(a,b)
    krnl = get_inner_kernel(a.dtype, b.dtype, dtype_out)
    return krnl(a, b).get().max()
    
vdot = inner

def weighted_inner(self, b, w):
    if w is None:
        return self.inner(b)  
    a = self.data
    dtype_out = _get_common_dtype(a, b)
    krnl = get_weighted_inner_kernel(a.dtype, b.dtype, w.dtype, dtype_out)
    return krnl(a, b, w).get().max()

# Define PYCUDA MAXLOC for both single and double precission ################## 
       
maxloc_preamble = """

    struct MAXLOCN{
        TTYPE max;
        LTYPE   loc;
        
        __device__
        MAXLOCN(){}
        
        __device__
        MAXLOCN(MAXLOCN const &src): max(src.max), loc(src.loc){}
        __device__
        MAXLOCN(MAXLOCN const volatile &src): max(src.max), loc(src.loc){}
        
        __device__
        MAXLOCN volatile &operator=( MAXLOCN const &src) volatile{
            max = src.max;
            loc = src.loc;
            return *this;
        }
    };
    
    __device__
    MAXLOCN maxloc_red(MAXLOCN a, MAXLOCN b){
        if (a.max > b.max)
            return a;
        else 
            return b;  
    }
    
    __device__
    MAXLOCN maxloc_start(){
        MAXLOCN t;
        t.max=0;
        t.loc=0;
        return t;
    }
    
    __device__
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

mls = LowerLatencyReductionKernel(maxloc_dtype_single, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(x[i], i)",
        arguments="float *x", preamble=maxloc_preamble_single)

mld = LowerLatencyReductionKernel(maxloc_dtype_double, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(x[i], i)",
        arguments="double *x", preamble=maxloc_preamble_double)
        
max_loc_map = {'single':mls,'double':mld}


amls = LowerLatencyReductionKernel(maxloc_dtype_single, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(abs(x[i]), i)",
        arguments="float *x", preamble=maxloc_preamble_single)

amld = LowerLatencyReductionKernel(maxloc_dtype_double, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(abs(x[i]), i)",
        arguments="double *x", preamble=maxloc_preamble_double)

amlsc = LowerLatencyReductionKernel(maxloc_dtype_single, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(abs(x[i]), i)",
        arguments="pycuda::complex<float> *x", preamble=maxloc_preamble_single)

amldc = LowerLatencyReductionKernel(maxloc_dtype_double, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(abs(x[i]), i)",
        arguments="pycuda::complex<double> *x", preamble=maxloc_preamble_double)

abs_max_loc_map = {'single':{ 'real':amls, 'complex':amlsc }, 'double':{ 'real':amld, 'complex':amldc }}

def zeros(length, dtype=np.float64):
    result = GPUArray(length, dtype=dtype)
    nwords = result.nbytes / 4
    pycuda.driver.memset_d32(result.gpudata, 0, nwords)
    return result

def ptr(self):
    return self._data.ptr

def dot(self, other):
    return pycuda.gpuarray.dot(self._data,other).get().max()

def min(self):
    return pycuda.gpuarray.min(self._data).get().max()

def abs_max_loc(self):
    maxloc = abs_max_loc_map[self.precision][self.kind](self._data)
    maxloc = maxloc.get()
    return float(maxloc['max']),int(maxloc['loc'])

def cumsum(self):
    tmp = self.data*1
    return icumsum(tmp)

def max(self):
    return pycuda.gpuarray.max(self._data).get().max()

def max_loc(self):
    maxloc = max_loc_map[self.precision](self._data)
    maxloc = maxloc.get()
    return float(maxloc['max']),int(maxloc['loc'])
    
def take(self, indices):
    if not isinstance(indices, pycuda.gpuarray.GPUArray):
        indices = pycuda.gpuarray.to_gpu(indices)
    return pycuda.gpuarray.take(self.data, indices)
    
def numpy(self):
    return self._data.get()
     
def _copy(self, self_ref, other_ref):
    if (len(other_ref) <= len(self_ref)) :
        from pycuda.elementwise import get_copy_kernel
        func = get_copy_kernel(self.dtype, other_ref.dtype)
        func.prepared_async_call(self_ref._grid, self_ref._block, None,
                self_ref.gpudata, other_ref.gpudata,
                self_ref.mem_size)
    else:
        raise RuntimeError("The arrays must the same length")

def _getvalue(self, index):
    return self._data.get()[index]
    
def sum(self):
    return pycuda.gpuarray.sum(self._data).get().max()
    
def clear(self):
    n32 = self.data.nbytes / 4
    pycuda.driver.memset_d32(self.data.gpudata, 0, n32)
    
def _scheme_matches_base_array(array):
    if isinstance(array, pycuda.gpuarray.GPUArray):
        return True
    else:
        return False

def _copy_base_array(array):
    data = pycuda.gpuarray.GPUArray((array.size), array.dtype)
    if len(array) > 0:
        pycuda.driver.memcpy_dtod(data.gpudata, array.gpudata, array.nbytes)
    return data

def _to_device(array):
    return pycuda.gpuarray.to_gpu(array)
    
   
   
