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
import pycuda.driver, pycbc.types, numpy
from pycuda.elementwise import ElementwiseKernel
from pycuda.tools import get_or_register_dtype, context_dependent_memoize, dtype_to_ctype
from pytools import match_precision, memoize_method
from pycuda.gpuarray import _get_common_dtype, empty, GPUArray
import pycuda.gpuarray
from mako.template import Template
from pycbc.types import Array
from pycuda.compiler import SourceModule

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


chisqkernel = Template("""
#include <stdio.h>
__global__ void power_chisq_at_points_${NP}(float2* corr, float2* outc, unsigned int N,
                                      %for p in range(NP):
                                        float phase${p},
                                      %endfor
                                      unsigned int* kmin, 
                                      unsigned int* kmax,
                                      unsigned int* bv, 
                                      unsigned int nbins){
    __shared__ unsigned int s;
    __shared__ unsigned int e;
    __shared__ float2 chisq[${NT} * ${NP}];
    
    // load integration boundaries (might not be bin boundaries if bin is large)
    if (threadIdx.x == 0){
        s = kmin[blockIdx.x];
        e = kmax[blockIdx.x];
    }
    
    % for p in range(NP):
        chisq[threadIdx.x + ${NT*p}].x = 0;
        chisq[threadIdx.x + ${NT*p}].y = 0;
    % endfor   
    __syncthreads();

    // calculate the chisq integral for each thread
    // sliding reduction for each thread from s, e
    for (int i = threadIdx.x + s; i < e; i += blockDim.x){
        float re, im;
        float2 qt = corr[i];
        
        %for p in range(NP):
            __sincosf(phase${p} * i, &im, &re);
            chisq[threadIdx.x + ${NT*p}].x += re * qt.x - im * qt.y;
            chisq[threadIdx.x + ${NT*p}].y += im * qt.x + re * qt.y;
        %endfor   
    }

    float x, y, x2, y2;
    // logarithmic reduction within thread block
    for (int j=${NT} / 2; j>=1; j/=2){
        if (threadIdx.x <j){
            %for p in range(NP):
                __syncthreads();
                x = chisq[threadIdx.x + ${NT*p}].x;
                y = chisq[threadIdx.x + ${NT*p}].y;
                x2 = chisq[threadIdx.x + j + ${NT*p}].x;
                y2 = chisq[threadIdx.x + j + ${NT*p}].y;
                 __syncthreads();
                chisq[threadIdx.x + ${NT*p}].x = x + x2;
                chisq[threadIdx.x + ${NT*p}].y = y + y2;
            %endfor
        }            
    }
  
    if (threadIdx.x == 0){
        % for p in range(NP):
            atomicAdd(&outc[bv[blockIdx.x] + nbins * ${p}].x, chisq[0 + ${NT*p}].x);
            atomicAdd(&outc[bv[blockIdx.x] + nbins * ${p}].y, chisq[0 + ${NT*p}].y);
        % endfor
    }

}
""")

_pchisq_cache = {}
def get_pchisq_fn(np):
    if np not in _pchisq_cache:
        nt = 256
        mod = SourceModule(chisqkernel.render(NT=nt, NP=np))
        fn = mod.get_function("power_chisq_at_points_%s" % (np))
        _pchisq_cache[np] = (fn, nt)
    return _pchisq_cache[np]

_bcache = {}
def get_cached_bin_layout(bins):
    key = id(bins)
    if key not in _bcache:
        bv, kmin, kmax = [], [], []
        for i in range(len(bins)-1):
            s, e = bins[i], bins[i+1]
            BS = 4096
            if (e - s) < BS:
                bv.append(i)
                kmin.append(s)
                kmax.append(e)
            else:
                k = list(numpy.arange(s, e, BS/2))
                kmin += k
                kmax += k[1:] + [e]
                bv += [i]*len(k)
        bv = pycuda.gpuarray.to_gpu_async(numpy.array(bv, dtype=numpy.uint32)) 
        kmin = pycuda.gpuarray.to_gpu_async(numpy.array(kmin, dtype=numpy.uint32))
        kmax = pycuda.gpuarray.to_gpu_async(numpy.array(kmax, dtype=numpy.uint32))  
        _bcache[key] = (kmin, kmax, bv) 
    return _bcache[key]

def shift_sum(corr, points, bins):
    corr = corr.data
    kmin, kmax, bv = get_cached_bin_layout(bins)
    nb = len(kmin)
    N = numpy.uint32(len(corr))
    nbins = numpy.uint32(len(bins) - 1)
    outc = pycuda.gpuarray.zeros((len(points), nbins), dtype=numpy.complex64)
    outp = outc.reshape(nbins * len(points))
    phase = [numpy.float32(p * -2.0 * numpy.pi / N) for p in points]

    np = len(points)
    while np > 0:
        if np >= 4:
            fn, nt = get_pchisq_fn(4)
            fn(corr, outp, N, phase[0], phase[1], phase[2], phase[3], 
               kmin, kmax, bv, nbins, block=(nt, 1, 1), grid=(nb, 1))
            outp = outp[4*nb:]
            phase = phase[4:]
            np -= 4    
            continue
        elif np >=3:
            fn, nt = get_pchisq_fn(3)
            fn(corr, outp, N, phase[0], phase[1], phase[2], 
                kmin, kmax, bv, nbins, block=(nt, 1, 1), grid=(nb, 1))
            np -= 3
            outp = outp[3*nb:]
            phase = phase[3:]
            continue
        elif np >=2:
            fn, nt = get_pchisq_fn(2)
            fn(corr, outp, N, phase[0], phase[1],
                kmin, kmax, bv, nbins, block=(nt, 1, 1), grid=(nb, 1))
            np -= 2
            outp = outp[2*nb:]
            phase=phase[2:]
            continue
        elif np == 1:
            fn, nt = get_pchisq_fn(1)
            fn(corr, outp, N, phase[0], 
                kmin, kmax, bv, nbins, block=(nt, 1, 1), grid=(nb, 1))
            np -= 1
            continue
    o = outc.get()
    return (o.conj() * o).sum(axis=1).real
    
    
