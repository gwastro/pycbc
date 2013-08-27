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
import pycbc.types
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
from mako.template import Template
from pycbc.types import Array

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
    
        final_result = kwargs.pop("result", None)
        
        f = self.stage1_func
        arg_types = self.stage1_arg_types
        stage1_args = args
        s1_invocation_args = [] 
        for arg in args:
            if isinstance(arg, GPUArray):
                s1_invocation_args.append(arg.gpudata)
            else:
                s1_invocation_args.append(arg)
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

            if block_count == 1 and final_result is not None:
                result = final_result
                
            f(grid_size, block_size, None,
                    *([result.gpudata, result2.gpudata]+s1_invocation_args+[seq_count, sz]),
                    shared_size=self.shared_size)

            if block_count == 1:
                return result
                
   
@context_dependent_memoize               
def get_shift_kernel(num_shifts):
    shift_preamble = Template("""
    struct shift_t${n}{
        % for i in range(n):
             float vr${i};
             float vi${i};
        % endfor
        
        __device__
        shift_t${n}(){}
        
        __device__
        shift_t${n}(shift_t${n} const &src): 
                    % for i in range(n-1):
                        vr${i}(src.vr${i}), 
                        vi${i}(src.vi${i}),
                    % endfor
                    vr${n-1}(src.vr${n-1}), 
                    vi${n-1}(src.vi${n-1})
                    {}
                   
        __device__
        shift_t${n}(shift_t${n} const volatile &src): 
                    % for i in range(n-1):
                        vr${i}(src.vr${i}), 
                        vi${i}(src.vi${i}),
                    % endfor
                    vr${n-1}(src.vr${n-1}), 
                    vi${n-1}(src.vi${n-1})
                    {}
        
        __device__
        shift_t${n} volatile &operator=( shift_t${n} const &src) volatile{
            % for i in range(n):
                 vr${i} = src.vr${i};
                 vi${i} = src.vi${i};
            % endfor
            return *this;
        }
    };
    


    __device__ shift_t${n} shift_red(shift_t${n} a, shift_t${n} b){
        % for i in range(n):
             a.vr${i} += b.vr${i};
             a.vi${i} += b.vi${i};
        % endfor
        return a;
    }

    __device__ shift_t${n} shift_start(){
        shift_t${n} t;
        % for i in range(n):
             t.vr${i}=0;
             t.vi${i}=0;
        % endfor
        return t;
    }

    __device__ shift_t${n} shift_map(pycuda::complex<float> x, 
                               % for i in range(n):
                                    float shift${i},
                               % endfor
                               float offset, float slen){
        shift_t${n} t; 
        float pphase = offset * 2 * 3.141592653  / slen;
        float  pr, pi;
        
        % for i in range(n):
            __sincosf(pphase * shift${i}, &pi, &pr);
           
            // Phase shift the input data (x) to correspond to a time shift
            t.vr${i} = x._M_re * pr - x._M_im * pi;
            t.vi${i} = x._M_re * pi + x._M_im * pr;  
        % endfor
        return t;
    }
    """).render(n = num_shifts)
    
    shift_map_args = ""
    shift_krnl_args = ""
    for i in range(num_shifts):
        shift_map_args += " shift%s," % i
        shift_krnl_args += " float shift%s, " % i        

    sd = np.dtype([("v1", np.complex64, num_shifts)])
    shift_t = get_or_register_dtype('shift_t%s' % num_shifts, sd)

    shift_krnl = LowerLatencyReductionKernel(shift_t, neutral="shift_start()",
                reduce_expr="shift_red(a, b)", map_expr="shift_map(x[i], " + shift_map_args + " offset+i, slen)",
                arguments="pycuda::complex<float> *x," + shift_krnl_args + "float offset, float slen ",
                preamble=shift_preamble)
                
    return shift_krnl

chisq_buf = pycbc.types.zeros(4096*256, dtype=np.complex64)

def shift_sum(v1, shifts, slen=None, offset=0):
    global chisq_buf
    vlen = len(v1)
    shifts = list(shifts)
    if slen is None:
        slen = vlen
    
    n = len(shifts)
    group_size = 10
    
    num_full_pass =  n // group_size
    remainder = n - group_size * num_full_pass
    result = chisq_buf

    for i in range(num_full_pass):
        f = result[i*group_size:(i+1)*group_size]
        shift_krnl = get_shift_kernel(group_size) 
        args = [v1.data] + shifts[i*group_size:(i+1)*group_size] + [offset, slen]
        shift_krnl(*args, result=f.data)
        
    if remainder > 0:
        f = result[group_size*num_full_pass:n]
        shift_krnl = get_shift_kernel(remainder) 
        args = [v1.data] + shifts[group_size*num_full_pass:n] + [offset, slen]
        shift_krnl(*args, result=f.data)
        
    return Array(result[0:n], copy=False)
   
    
    
