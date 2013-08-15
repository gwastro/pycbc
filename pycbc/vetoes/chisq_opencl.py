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
import pycbc.scheme
from pyopencl.reduction import ReductionKernel
from pyopencl.tools import dtype_to_ctype, get_or_register_dtype
from pyopencl.elementwise import ElementwiseKernel, complex_dtype_to_name

def chisq_accum_bin(chisq, q):
    chisq += q.squared_norm()


cfloat  = complex_dtype_to_name(numpy.complex64)
cdouble = complex_dtype_to_name(numpy.complex128)

get_or_register_dtype('cfloat', numpy.complex64)
get_or_register_dtype('cdouble', numpy.complex128)

shift_preamble = """
typedef struct {
    cfloat_t vals[2];
}shift_t;

__global shift_t shift_red(shift_t a, shift_t b, int slen){
    for (int i=0; i < slen; i++){
        a.vals[i] += b.vals[i];
    }
    return a;
}

__global shift_t shift_start(){
    shift_t t;
    return t;
}

__global
shift_t shift_map(cfloat_t x, __global float* s, __global float* o, int slen){
    shift_t t;
    for (int i=0; i< slen; i++){
        t.vals[i] = x;
    }
    return t;
}
"""
with pycbc.scheme.OpenCLScheme():

    sd = numpy.dtype([("v1", numpy.complex64), ("v2", numpy.complex64)])
    shift_t = get_or_register_dtype('shift_t', sd)

    shift_krnl = ReductionKernel(pycbc.scheme.mgr.state.context,
                shift_t, neutral="shift_start()",
                reduce_expr="shift_red(a, b, slen)", map_expr="shift_map(x[i], s, o, slen)",
                arguments="__global cfloat_t *x, __global float* s, __global float* o, int slen",
                preamble=shift_preamble)

    import pycbc.types
    a = pycbc.types.zeros(2**20, dtype=numpy.complex64) + 3 + 3j
    s = pycbc.types.zeros(10, dtype=numpy.float32)
    o = pycbc.types.zeros(10, dtype=numpy.float32)
    slen = float(2)
    shift_krnl(a.data, s.data, o.data, slen)





