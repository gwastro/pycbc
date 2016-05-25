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
import numpy, mako.template
from pycuda import gpuarray
from pycuda.tools import dtype_to_ctype
from pycuda.elementwise import ElementwiseKernel
from pycuda.compiler import SourceModule
from .events import _BaseThresholdCluster
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

tkernel1 = mako.template.Template("""
#include <stdio.h>

__global__ void threshold_and_cluster(float2* in, float2* outv, int* outl, int window, float threshold){
    int s = window * blockIdx.x;
    int e = s + window;
    
    // shared memory for chuck size candidates
    __shared__ float svr[${chunk}];
    __shared__ float svi[${chunk}];
    __shared__ int sl[${chunk}];
    
    // shared memory for the warp size candidates
    __shared__ float svv[32];
    __shared__ int idx[32];
    
    int ml = -1;
    float mvr = 0;
    float mvi = 0;
    float re;
    float im;
    
    // Iterate trought the entire window size chunk and find blockDim.x number
    // of candidates
    for (int i = s + threadIdx.x; i < e; i += blockDim.x){
        re = in[i].x;
        im = in[i].y;
        if ((re * re + im * im) > (mvr * mvr + mvi * mvi)){
            mvr = re;
            mvi = im;
            ml = i;
        }
    }
    
    // Save the candidate from this thread to shared memory
    svr[threadIdx.x] = mvr;
    svi[threadIdx.x] = mvi;
    sl[threadIdx.x] = ml;
    
    syncthreads();
    if (threadIdx.x < 32){
        int tl = threadIdx.x;
        
        // Now that we have all the candiates for this chunk in shared memory
        // Iterate through in the warp size to reduce to 32 candidates
        for (int i = threadIdx.x; i < ${chunk}; i += 32){
            re = svr[i];
            im = svi[i];
            if ((re * re + im * im) > (mvr * mvr + mvi * mvi)){
                tl = i;
                mvr = re;
                mvi = im;
            }
        }
        
        // Store the 32 candidates into shared memory
        svv[threadIdx.x] = svr[tl] * svr[tl] + svi[tl] * svi[tl];
        idx[threadIdx.x] = tl;
       
        // Find the 1 candidate we are looking for using a manual log algorithm
        if ((threadIdx.x < 16) && (svv[threadIdx.x] < svv[threadIdx.x + 16])){
            svv[threadIdx.x] = svv[threadIdx.x + 16];
            idx[threadIdx.x] = idx[threadIdx.x + 16];
        }
              
        if ((threadIdx.x < 8) && (svv[threadIdx.x] < svv[threadIdx.x + 8])){
            svv[threadIdx.x] = svv[threadIdx.x + 8];
            idx[threadIdx.x] = idx[threadIdx.x + 8];
        }
    
        if ((threadIdx.x < 4) && (svv[threadIdx.x] < svv[threadIdx.x + 4])){
            svv[threadIdx.x] = svv[threadIdx.x + 4];
            idx[threadIdx.x] = idx[threadIdx.x + 4];
        }
        
        if ((threadIdx.x < 2) && (svv[threadIdx.x] < svv[threadIdx.x + 2])){
            svv[threadIdx.x] = svv[threadIdx.x + 2];
            idx[threadIdx.x] = idx[threadIdx.x + 2];
        }
        
        
        // Save the 1 candidate maximum and location to the output vectors
        if (threadIdx.x == 0){
            if (svv[threadIdx.x] < svv[threadIdx.x + 1]){
                idx[0] = idx[1];
                svv[0] = svv[1];
            }
            
            if (svv[0] > threshold){   
                tl = idx[0];
                outv[blockIdx.x].x = svr[tl];
                outv[blockIdx.x].y = svi[tl];
                outl[blockIdx.x] = sl[tl];
            } else{
                outl[blockIdx.x] = -1;
            }
        }        
    }
}
""")

tkernel2 = mako.template.Template("""
#include <stdio.h>
__global__ void threshold_and_cluster2(float2* outv, int* outl, float threshold, int window){
    __shared__ int loc[${blocks}];
    __shared__ float val[${blocks}];
    
    int i = threadIdx.x;

    int l = outl[i];
    loc[i] = l;
    
    if (l == -1)
        return;
        
    val[i] = outv[i].x * outv[i].x + outv[i].y * outv[i].y;

    
    // Check right
    if ( (i < (${blocks} - 1)) && (val[i + 1] > val[i]) ){
        outl[i] = -1;
        return;
    }
    
    // Check left
    if ( (i > 0) && (val[i - 1] > val[i]) ){
        outl[i] = -1;
        return;
    }
}
""")

tfn_cache = {}
def get_tkernel(slen, window):
    if window < 32:
        raise ValueError("GPU threshold kernel does not support a window smaller than 32 samples")
    
    elif window <= 4096:
        nt = 128
    elif window <= 16384:
        nt = 256
    elif window <= 32768:
        nt = 512
    else:
        nt = 1024
    
    nb = int(numpy.ceil(slen / float(window)))

    if nb > 1024:
        raise ValueError("More than 1024 blocks not supported yet")

    try:
        return tfn_cache[(nt, nb)], nt, nb
    except KeyError:
        mod = SourceModule(tkernel1.render(chunk=nt))
        mod2 = SourceModule(tkernel2.render(blocks=nb))
        fn = mod.get_function("threshold_and_cluster")
        fn.prepare("PPPif")
        fn2 = mod2.get_function("threshold_and_cluster2")
        fn2.prepare("PPfi")
        tfn_cache[(nt, nb)] = (fn, fn2)
        return tfn_cache[(nt, nb)], nt, nb
    
def threshold_and_cluster(series, threshold, window):
    outl = tl.gpudata
    outv = tv.gpudata
    slen = len(series)
    series = series.data.gpudata
    (fn, fn2), nt, nb = get_tkernel(slen, window)
    threshold = numpy.float32(threshold * threshold)
    window = numpy.int32(window)
    
    cl = loc[0:nb]
    cv = val[0:nb]
    
    fn.prepared_call((nb, 1), (nt, 1, 1), series, outv, outl, window, threshold,)
    fn2.prepared_call((1, 1), (nb, 1, 1), outv, outl, threshold, window)   
    pycbc.scheme.mgr.state.context.synchronize()
    w = (cl != -1)
    return cv[w], cl[w]
    
class CUDAThresholdCluster(_BaseThresholdCluster):
    def __init__(self, series):
        self.series = series.data.gpudata

        self.outl = tl.gpudata
        self.outv = tv.gpudata
        self.slen = len(series)

    def threshold_and_cluster(self, threshold, window):
        threshold = numpy.float32(threshold * threshold)
        window = numpy.int32(window)

        (fn, fn2), nt, nb = get_tkernel(self.slen, window)     
        fn = fn.prepared_call
        fn2 = fn2.prepared_call
        cl = loc[0:nb]
        cv = val[0:nb]

        fn((nb, 1), (nt, 1, 1), self.series, self.outv, self.outl, window, threshold,)
        fn2((1, 1), (nb, 1, 1), self.outv, self.outl, threshold, window)   
        pycbc.scheme.mgr.state.context.synchronize()
        w = (cl != -1)
        return cv[w], cl[w]
    
def _threshold_cluster_factory(series):
    return CUDAThresholdCluster

