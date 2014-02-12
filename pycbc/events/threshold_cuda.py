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
from pycuda.compiler import SourceModule
import pycbc.scheme

threshold_op = """
    if (i == 0)
        bn[0] = 0;

    pycuda::complex<float> val = in[i];
    if ( abs(val) > threshold){
        int n_w = atomicAdd(bn, 1);
        outv[n_w] = val;
        outl[n_w] = i;
    }

"""

threshold_kernel = ElementwiseKernel(
            " %(tp_in)s *in, %(tp_out1)s *outv, %(tp_out2)s *outl, %(tp_th)s threshold, %(tp_n)s *bn" % {
                "tp_in": dtype_to_ctype(numpy.complex64),
                "tp_out1": dtype_to_ctype(numpy.complex64),
                "tp_out2": dtype_to_ctype(numpy.uint32),
                "tp_th": dtype_to_ctype(numpy.float32),
                "tp_n": dtype_to_ctype(numpy.uint32),
                },
            threshold_op,
            "getstuff")
            
import pycuda.driver as drv
n = drv.pagelocked_empty((1), numpy.uint32, mem_flags=drv.host_alloc_flags.DEVICEMAP)
nptr = numpy.intp(n.base.get_device_pointer())

val = drv.pagelocked_empty((4096*256), numpy.complex64, mem_flags=drv.host_alloc_flags.DEVICEMAP)
vptr = numpy.intp(val.base.get_device_pointer())

loc = drv.pagelocked_empty((4096*256), numpy.int32, mem_flags=drv.host_alloc_flags.DEVICEMAP)
lptr = numpy.intp(loc.base.get_device_pointer())
            
class T():
    pass

tn = T()
tv = T()
tl = T()
tn.gpudata = nptr
tv.gpudata = vptr
tl.gpudata = lptr
tn.flags = tv.flags = tl.flags = n.flags

def threshold(series, value):
    threshold_kernel(series.data, tv, tl, value, tn)
    pycbc.scheme.mgr.state.context.synchronize()
    n0 = n[0]
    srt = numpy.argsort(loc[0:n0])
    return loc[srt], val[srt]
  

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

