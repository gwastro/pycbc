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
from pycbc.types import zeros, Array
from pycuda.gpuarray import to_gpu
from pycuda.tools import get_or_register_dtype, dtype_to_ctype
from pycuda.elementwise import ElementwiseKernel
from events import complex64_subset
from pycuda.compiler import SourceModule

complex64_subset = get_or_register_dtype("event", dtype=complex64_subset)

preamble = """
    #include <stdio.h>
    struct event{
        pycuda::complex<float> val;
        long loc;
    };
    
    """

threshold_op = """
    if (i == 0)
        bn[0] = 0;

    pycuda::complex<float> val = in[i];
    event nv;
    if ( abs(val) > threshold){
        nv.val = val;
        nv.loc = i;
        int n_w = atomicAdd(bn, 1) ;
        out[n_w] = nv;
    }

"""

threshold_kernel = ElementwiseKernel(
            " %(tp_in)s *in, %(tp_out)s *out, %(tp_th)s threshold, %(tp_n)s *bn" % {
                "tp_in": dtype_to_ctype(numpy.complex64),
                "tp_out": dtype_to_ctype(complex64_subset),
                "tp_th": dtype_to_ctype(numpy.float32),
                "tp_n": dtype_to_ctype(numpy.int32),
                },
            threshold_op,
            "getstuff", preamble=preamble)
            
n_events = numpy.zeros(1, dtype=numpy.int64)
n_events = to_gpu(n_events)
buffer_vec = numpy.zeros(4096*2048, dtype=complex64_subset)
buffer_vec = to_gpu(buffer_vec)
            
def threshold(series, value):
    threshold_kernel(series.data, buffer_vec, value, n_events)
    n = n_events.get()[0]
    return numpy.sort(buffer_vec[0:n].get(), order='loc')
  


threshold_cluster_mod = SourceModule("""
#include<pycuda-complex.hpp>
#include<stdio.h>

__global__ void threshold_and_cluster(pycuda::complex<float>* series, float threshold, unsigned int window, unsigned int blen, unsigned int tlen){
    unsigned int start = blockIdx.x * blen;
    unsigned int end = start + blen;
    
    //printf("start: %i  end  : %i block: %i\\n", start, end, blockIdx.x);
    
    if (end >= tlen)
        end = tlen;
        
    float consider_val = -1;
    unsigned int consider_index=start;
    
    float test_val;
        
    for (unsigned int i=start; i < start + blen; i++){
        test_val = abs(series[i]);
        
        if (test_val < threshold)
            continue;

        if (i > consider_index + window){
            series[consider_index] = -1;   
       }else 
            if (consider_val > test_val)
                continue;  
                
        consider_val = test_val;
        consider_index = i;
    }
}
""")

threshold_cluster_krnl = threshold_cluster_mod.get_function("threshold_and_cluster")
  

def threshold_and_cluster(series, threshold, window):
     
    threshold = numpy.float32(threshold)
    window = numpy.uint32(window)
    tlen = numpy.uint32(len(series))
    
    blocklen = numpy.uint32(tlen / 32)
    numblocks = int(tlen / blocklen)
    
    print numblocks, blocklen, tlen
    
    threshold_cluster_krnl(series.data.gpudata, threshold, window, blocklen, tlen, 
                           block=(1, 1, 1), grid=(numblocks, 1))

