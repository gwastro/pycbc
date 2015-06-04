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
from pycuda import gpuarray, driver
from pycuda.gpuarray import to_gpu, empty
from pycuda.tools import get_or_register_dtype, dtype_to_ctype
from pycuda.elementwise import ElementwiseKernel
from pycuda.scan import ExclusiveScanKernel
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

def standard_threshold(series, value):
    threshold_kernel(series.data, tv, tl, value, tn)
    pycbc.scheme.mgr.state.context.synchronize()
    n0 = n[0]
    srt = numpy.argsort(loc[0:n0])
    return loc[srt], val[srt]
  

# An attempt at faster thresholding
mark_above = ElementwiseKernel(
        " %(tp_in)s *in, %(tp_mark)s *mark, %(tp_th)s threshold " % {
            "tp_in": dtype_to_ctype(numpy.complex64),
            "tp_mark": dtype_to_ctype(numpy.uint32),
            "tp_th": dtype_to_ctype(numpy.float32) },
        """
            pycuda::complex<float> val = in[i];
            if ( abs(val) > threshold )
                { mark[i] = 1; }
            else
                { mark[i] = 0; }
        """,
        "mark_above_threshold")

scan_sum = ExclusiveScanKernel(numpy.uint32, "a+b", 0)

stream_compact = ElementwiseKernel(
        " %(tp_in)s *in, %(tp_out1)s *out1, %(tp_out2)s *out2, %(tp_wl)s *write_loc, %(tp_th)s threshold, %(tp_n)s *bn " % {
        "tp_in": dtype_to_ctype(numpy.complex64),
        "tp_out1": dtype_to_ctype(numpy.complex64),
        "tp_out2": dtype_to_ctype(numpy.uint32),
        "tp_wl": dtype_to_ctype(numpy.uint32),
        "tp_th": dtype_to_ctype(numpy.float32),
        "tp_n": dtype_to_ctype(numpy.uint32) },
        """
            pycuda::complex<float> val = in[i];
            if ( abs(val) > threshold )
            {
                out1[write_loc[i]] = in[i];
                out2[write_loc[i]] = i;
            }
            if ( i == n-1 ) { bn[0] = write_loc[i]; }
        """,
        "stream_compact")

write_loc = empty((4096*256), numpy.uint32)
def threshold_by_scan(series, value):
    mark_above(series.data, write_loc, value)
    scan_sum(write_loc)
    stream_compact(series.data, tv, tl, write_loc, value, tn)
    pycbc.scheme.mgr.state.context.synchronize()
    n0 = n[0]
    return loc[0:n0], val[0:n0]


# Do thresholding using the GPU only for the abs function
threshold_cpu_op = """ 
    pycuda::complex<float> val = in[i];
    outv[i] = val;
    if ( abs(val) > threshold)
        { outl[i] = 1; }
    else
        { outl[i] = 0; }
"""

threshold_cpu_kernel = ElementwiseKernel(
            " %(tp_in)s *in, %(tp_out1)s *outv, %(tp_out2)s *outl, %(tp_th)s threshold " % {
                "tp_in": dtype_to_ctype(numpy.complex64),
                "tp_out1": dtype_to_ctype(numpy.complex64),
                "tp_out2": dtype_to_ctype(numpy.int32),
                "tp_th": dtype_to_ctype(numpy.float32)
                },
            threshold_cpu_op,
            "getstuff")

def threshold_with_gpu_abs(series, value):
    threshold_cpu_kernel(series.data, tv, tl, value)
    pycbc.scheme.mgr.state.context.synchronize()
    tgt = numpy.where(loc == 1)[0]
    return tgt, val[tgt]
    
# Do the thresholding using a block-level scan operation within the kernel
nt = 128
gs = 32

kernel2 = """
#include <stdio.h>
#define GS %s
#define NT %s
 
__global__ void threshold_comb(float2* in, unsigned int* outl, unsigned int* outl2, float2* outv,  unsigned int* bn, unsigned int* n){ 

    // should be number of blocks in the threshold_seg kernel, we guess here instead
    // which is ok as long as it is greater
    __shared__ int loc[2048];
    int li = threadIdx.x;  
        
    int tt = 0;
    loc[li] = bn[li];
    __syncthreads();   
       
    for (int j=0; j<blockIdx.x; j++)
        tt += loc[j];
    
    int rn = loc[blockIdx.x];
    for (int r=li; r<rn; r+=blockDim.x){
        int loc = outl[blockIdx.x*GS*NT + r];
        outl2[tt + r] =  loc;
        outv[tt + r] = in[loc];
    }

    if (blockIdx.x == (gridDim.x-1) && threadIdx.x == 0)
        n[0] = tt + rn;
}           
"""  % (gs, nt)
mod = SourceModule(kernel2)
stuff2 = mod.get_function("threshold_comb")

kernel = """
#include <stdio.h>

#define GS %s
#define NT %s

__global__ void threshold_seg(float2 *in, unsigned int* outl, unsigned int* bn, float threshold, unsigned int length){  
    float2 val;
    int above, tmp;
    int s = blockIdx.x * blockDim.x * GS;
    int i = s + threadIdx.x;
    int region = s + GS * blockDim.x;
    int li = threadIdx.x;
    int lio = li + NT;
    
    volatile __shared__ int c;
    volatile __shared__ unsigned int loc[NT * 2];
    
    if (li == 0)
        c = 0;
   
    for (; i<region; i+=blockDim.x){
    
        // Do the thresholding for a NUM_THREAD section of memory
        // Above threshold is stored as a 1 in shared memory       
        if (i < length){
            val = in[i];
            above =  (val.x*val.x + val.y*val.y > threshold)?1:0;
        } else
            above =0;
        
        loc[lio] = above;
        loc[li]=0;
        __syncthreads();
              
        // Do a cumulative sum of the shared memory above threshold flags
        // The cumsum is used to find where to store the locations of the
        // above threshold regions within the scope of this block
        // Be caeful ever touching this section, the syncthreads are really
        // important!!
        
        // Note that this is hardcoded to the block size (number of threads)
              
        tmp = loc[lio] +  loc[lio - 1 ];
        __syncthreads();
        loc[lio] = tmp;
        __syncthreads();
        
        tmp = loc[lio] +  loc[lio - 2 ];
        __syncthreads();
        loc[lio] = tmp;
        __syncthreads();
        
        tmp = loc[lio] +  loc[lio - 4 ];
        __syncthreads();
        loc[lio] = tmp;
        __syncthreads();
        
        tmp = loc[lio] +  loc[lio - 8 ];
        __syncthreads();
        loc[lio] = tmp;
        __syncthreads();
        
        tmp = loc[lio] +  loc[lio - 16 ];
        __syncthreads();
        loc[lio] = tmp;
        __syncthreads();
        
        tmp = loc[lio] +  loc[lio - 32 ];
        __syncthreads();
        loc[lio] = tmp;
        __syncthreads();
        
        tmp = loc[lio] +  loc[lio - 64];
        __syncthreads();
        loc[lio] = tmp;
        __syncthreads();
        
        
        //tmp = loc[lio] +  loc[lio - 128 ];
        //__syncthreads();
        //loc[lio] = tmp;
        //__syncthreads();

        // Store the index value of a trigger now that we know where to put it
        if (above)
            outl[c + s + loc[lio] - 1] = i;   
        
        if (threadIdx.x == 0)
            c += loc[NT*2-1];
    }
    if (threadIdx.x == 0)
        bn[blockIdx.x] = c;

}
""" % (gs, nt)
mod = SourceModule(kernel)
stuff = mod.get_function("threshold_seg")

# Workspace memory needed for this kernel
_dn = None
_n = None
_bn = None
_loc_tmp = None
_loc_out = None
_val_out = None
_val = None
_loc = None

def threshold_integrated(series, value):
    global _dn, _n, _bn, _loc_tmp, _loc_out, _val_out, _loc, _val
        
    t = numpy.float32(value**2)
    nb = int(numpy.ceil(float(len(series))/nt/gs))
    
    if _bn is None or len(_bn) < nb:
        _bn = gpuarray.zeros(nb, dtype=numpy.uint32)
        
    if _n is None:
        _n = driver.pagelocked_empty((1), numpy.uint32, mem_flags=drv.host_alloc_flags.DEVICEMAP)
        ptr = numpy.intp(_n.base.get_device_pointer())
        class T():
            pass
        _dn = T()
        _dn.gpudata = ptr
        _dn.flags = _n.flags
        
    if _loc_tmp is None or len(series) > len(_loc_tmp):
        _loc_tmp = gpuarray.zeros(len(series), dtype=numpy.uint32)
        _loc_out = gpuarray.zeros(len(series), dtype=numpy.uint32)
        _val_out = gpuarray.zeros(len(series), dtype=series.dtype)
        _val = driver.pagelocked_empty((4096*256), numpy.complex64)
        _loc = driver.pagelocked_empty((4096*256), numpy.uint32)
    
    #Do the thresholding by block
    stuff(series.data, _loc_tmp, _bn, t, numpy.uint32(len(series)), block=(nt, 1, 1), grid=(nb, 1))
    
    # Recombine the blocks into a final output
    stuff2(series.data, _loc_tmp, _loc_out, _val_out, _bn, _dn, block=(nb, 1, 1), grid=(nb, 1))
    
    # We need to get the data back now
    pycbc.scheme.mgr.state.context.synchronize()
    if _n != 0: 
        driver.memcpy_dtoh_async(_val[0:_n], _val_out.gpudata)
        driver.memcpy_dtoh_async(_loc[0:_n], _loc_out.gpudata)
        pycbc.scheme.mgr.state.context.synchronize()
    return _loc[0:_n], _val[0:_n]

# Do thresholding entirely on the CPU
def threshold_on_cpu(series, value):
    arr = series.data.get()
    locs = numpy.where(arr.real**2 + arr.imag**2 > value**2)[0]
    return locs, arr[locs]


# Select which of the thresholding methods we want to use out of the above.
threshold = threshold_by_scan

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
    assert window > 0, 'Clustering window length is not positive'
    threshold = numpy.float32(threshold)
    window = numpy.uint32(window)
    tlen = numpy.uint32(len(series))
    
    blocklen = numpy.uint32(tlen / 32)
    numblocks = int(tlen / blocklen)
    
    print numblocks, blocklen, tlen
    
    threshold_cluster_krnl(series.data.gpudata, threshold, window, blocklen, tlen, 
                           block=(1, 1, 1), grid=(numblocks, 1))

