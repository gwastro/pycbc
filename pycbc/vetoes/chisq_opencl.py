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
import numpy
from pycbc.types import Array
import pycbc.scheme
from pyopencl import array as clarray
from pyopencl.reduction import ReductionKernel
from pyopencl.tools import dtype_to_ctype, get_or_register_dtype
from pyopencl.elementwise import ElementwiseKernel, complex_dtype_to_name
from pyopencl.tools import context_dependent_memoize
from mako.template import Template

def chisq_accum_bin(chisq, q):
    chisq += q.squared_norm()
    
cfloat  = complex_dtype_to_name(numpy.complex64)
cdouble = complex_dtype_to_name(numpy.complex128)

get_or_register_dtype('cfloat', numpy.complex64)
get_or_register_dtype('cdouble', numpy.complex128)

class LowerLatencyReductionKernel(ReductionKernel):
    def __call__(self, *args, **kwargs):
        MAX_GROUP_COUNT = 1024
        SMALL_SEQ_COUNT = 4

        from pyopencl.array import empty

        stage_inf = self.stage_1_inf
        final_result = kwargs.pop("result", None)
        queue = kwargs.pop("queue", None)
        wait_for = kwargs.pop("wait_for", None)
        return_event = kwargs.pop("return_event", False)

        if kwargs:
            raise TypeError("invalid keyword argument to reduction kernel")

        stage1_args = args

        while True:
            invocation_args = []
            vectors = []

            from pyopencl.tools import VectorArg
            for arg, arg_tp in zip(args, stage_inf.arg_types):
                if isinstance(arg_tp, VectorArg):
                    if not arg.flags.forc:
                        raise RuntimeError("ReductionKernel cannot "
                                "deal with non-contiguous arrays")

                    vectors.append(arg)
                    invocation_args.append(arg.base_data)
                    if arg_tp.with_offset:
                        invocation_args.append(arg.offset)
                else:
                    invocation_args.append(arg)

            repr_vec = vectors[0]
            sz = repr_vec.size

            if queue is not None:
                use_queue = queue
            else:
                use_queue = repr_vec.queue

            if sz <= stage_inf.group_size*SMALL_SEQ_COUNT*MAX_GROUP_COUNT:
                total_group_size = SMALL_SEQ_COUNT*stage_inf.group_size
                group_count = (sz + total_group_size - 1) // total_group_size
                seq_count = SMALL_SEQ_COUNT
            else:
                group_count = MAX_GROUP_COUNT
                macrogroup_size = group_count*stage_inf.group_size
                seq_count = (sz + macrogroup_size - 1) // macrogroup_size

            if group_count == 1:
                if final_result is None:
                    result = empty(use_queue,
                        (), self.dtype_out,
                        allocator=repr_vec.allocator)
                else:
                    result = final_result
            else:
                result = empty(use_queue, (group_count,), self.dtype_out)

            last_evt = stage_inf.kernel(
                    use_queue,
                    (group_count*stage_inf.group_size,),
                    (stage_inf.group_size,),
                    *([result.data]+invocation_args+[seq_count, sz]),
                    **dict(wait_for=wait_for))

            wait_for = [last_evt]

            if group_count == 1:
                if return_event:
                    return result, last_evt
                else:
                    return result
            else:
                stage_inf = self.stage_2_inf
                args = (result,) + stage1_args
 
@context_dependent_memoize               
def get_shift_kernel(num_shifts):
    shift_preamble = Template("""
    typedef struct {
        % for i in range(n):
             cfloat_t v${i};
        % endfor
    }shift_t${n};

    __global shift_t${n} shift_red(shift_t${n} a, shift_t${n} b){
        % for i in range(n):
             a.v${i} += b.v${i};
        % endfor
        return a;
    }

    __global shift_t${n} shift_start(){
        shift_t${n} t;
        % for i in range(n):
             t.v${i}=0;
        % endfor
        return t;
    }

    __global shift_t${n} shift_map(cfloat_t x, 
                               % for i in range(n):
                                    float shift${i},
                               % endfor
                               float offset, float slen){
        shift_t${n} t; 
        float pphase = offset * 2 * 3.141592653  / slen;
        float phase, pr, pi;
        
        % for i in range(n):
            phase = pphase * shift${i};
            pi = sincos(phase, &pr);
           
            // Phase shift the input data (x) to correspond to a time shift
            t.v${i}.x = x.x * pr - x.y * pi;
            t.v${i}.y = x.x * pi + x.y * pr;  
        % endfor
        return t;
    }
    """).render(n = num_shifts)
    
    shift_map_args = ""
    shift_krnl_args = ""
    for i in range(num_shifts):
        shift_map_args += " shift%s," % i
        shift_krnl_args += " float shift%s, " % i        

    sd = numpy.dtype([("v1", numpy.complex64, num_shifts)])
    shift_t = get_or_register_dtype('shift_t%s' % num_shifts, sd)

    shift_krnl = LowerLatencyReductionKernel(pycbc.scheme.mgr.state.context,
                shift_t, neutral="shift_start()",
                reduce_expr="shift_red(a, b)", map_expr="shift_map(x[i], " + shift_map_args + " offset+i, slen)",
                arguments="__global cfloat_t *x," + shift_krnl_args + "float offset, float slen ",
                preamble=shift_preamble)
                
    return shift_krnl


def shift_sum(v1, shifts, slen=None, offset=0):
    vlen = len(v1)
    shifts = list(shifts)
    if slen is None:
        slen = vlen

    v = []
    
    n = len(shifts)
    group_size = 4
    
    num_full_pass =  n // group_size
    remainder = n - group_size * num_full_pass

    for i in range(num_full_pass):
        f = pycbc.types.zeros(group_size, dtype=numpy.complex64)
        shift_krnl = get_shift_kernel(group_size) 
        args = [v1.data] + shifts[i*group_size:(i+1)*group_size] + [offset, slen]
        v.append(shift_krnl(*args, result=f.data))
        
    if remainder > 0:
        f = pycbc.types.zeros(remainder, dtype=numpy.complex64)
        shift_krnl = get_shift_kernel(remainder) 
        args = [v1.data] + shifts[group_size*num_full_pass:n] + [offset, slen]
        v.append(shift_krnl(*args, result=f.data))
     
    result = clarray.concatenate(v)   
        
    return Array(result, copy=False)

