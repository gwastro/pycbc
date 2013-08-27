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
import pycuda.driver
from pycuda.elementwise import ElementwiseKernel
from pycuda.reduction import ReductionKernel
from pycuda.tools import get_or_register_dtype
from pycuda.tools import context_dependent_memoize
from pycuda.tools import dtype_to_ctype
from pytools import match_precision, memoize_method
from pycuda.gpuarray import _get_common_dtype, empty, GPUArray
import pycuda.gpuarray
from pycuda.scan import InclusiveScanKernel
import numpy as np

@context_dependent_memoize
def get_accum_diff_sq_kernel(dtype_x, dtype_z):
    return ElementwiseKernel(
            "%(tp_a)s *x,  %(tp_c)s *z" % {
                "tp_a": dtype_to_ctype(dtype_x),
                "tp_c": dtype_to_ctype(dtype_z),
                },
            "x[i] += norm(z[i]) ",
            "chisq_accum")    
 
def chisq_accum_bin(chisq, q):
    krnl = get_accum_diff_sq_kernel(chisq.dtype, q.dtype)
    krnl(chisq.data, q.data)
    
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
        arg_types = self.stage1_arg_types
        stage1_args = args
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
            arg_types = self.stage2_arg_types
            sz = result.size
            result2 = result
            result, block_count, seq_count, grid_size, block_size = call_prepare(self, sz, args[0].allocator)

            f(grid_size, block_size, None,
                    *([result.gpudata, result2.gpudata]+s1_invocation_args+[seq_count, sz]),
                    shared_size=self.shared_size)

            if block_count == 1:
                return result
   
    
    
    
