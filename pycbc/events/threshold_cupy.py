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

import cupy as cp
import functools
import mako.template
from .eventmgr import _BaseThresholdCluster
import pycbc.scheme

val = None
loc = None

# https://stackoverflow.com/questions/77798014/cupy-rawkernel-cuda-error-not-found-named-symbol-not-found-cupy

tkernel1 = mako.template.Template("""
extern "C" __global__ void threshold_and_cluster(float2* in, float2* outv, int* outl, int window, float* thresholds, int series_length, int analyse_start) {
    int batch_idx = blockIdx.y;  // Batch index
    int s = batch_idx * series_length + window * blockIdx.x + analyse_start;
    int e = s + window;
    float threshold = thresholds[batch_idx];

    // Shared memory remains unchanged, but it now processes series per batch index
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
                                  
    __syncthreads();

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
                outv[batch_idx*${blockmemsize} + blockIdx.x].x = svr[tl];
                outv[batch_idx*${blockmemsize} + blockIdx.x].y = svi[tl];
                outl[batch_idx*${blockmemsize} + blockIdx.x] = sl[tl] % ${slen} - analyse_start;
                } else {
                outl[batch_idx*${blockmemsize} + blockIdx.x] = -1;
            }
        }
    }
}
""")

tkernel2 = mako.template.Template("""
extern "C" __global__ void threshold_and_cluster2(float2* outv, int* outl, float* thresholds, int window){
    __shared__ int loc[${blocks}];
    __shared__ float val[${blocks}];
                                  
    int i = threadIdx.x;
    int posi = i % ${blockmemsize};
    float threshold = thresholds[i / ${blockmemsize}];
                                  
    int l = outl[i];
    loc[i] = l;

    if (l == -1)
        return;

    val[i] = outv[i].x * outv[i].x + outv[i].y * outv[i].y;


    // Check right
    if ( (posi < (${blocksize} - 1)) && (val[i + 1] > val[i]) ){
        outl[i] = -1;
        return;
    }

    // Check left
    if ( (posi > 0) && (val[i - 1] > val[i]) ){
        outl[i] = -1;
        return;
    }
}
""")

@functools.lru_cache(maxsize=None)
def get_tkernel(slen, alen, window, block_mem_size=None, batch_size=None):
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

    nb = int(cp.ceil(alen / float(window)))

    if nb > 1024:
        raise ValueError("More than 1024 blocks not supported yet")

    if block_mem_size is not None:
        blocks = block_mem_size * batch_size
    else:
        blocks = nb
        block_mem_size = nb

    fn = cp.RawKernel(
        tkernel1.render(chunk=nt, slen=slen, blockmemsize=block_mem_size),
        'threshold_and_cluster',
        backend='nvcc'
    )
    fn2 = cp.RawKernel(
        tkernel2.render(blocks=nb, blocksize=nb, blockmemsize=block_mem_size),
        'threshold_and_cluster2',
        backend='nvcc'
    )
    return (fn, fn2), nt, nb

def threshold_and_cluster(series_batch, threshold, window):
    raise NotImplementedError("Needs writing properly")
    # Not sure this function is accessed easily. However, right now, it has not
    # been properly written. A starting point is here though!
    global val
    global loc
    batch_size, series_length = series_batch.shape
    
    if val is None or val.size < batch_size * 4096 * 256:
        val = cp.zeros((batch_size, 4096 * 256), dtype=cp.complex64)
    if loc is None or loc.size < batch_size * 4096 * 256:
        loc = cp.zeros((batch_size, 4096 * 256), dtype=cp.int32)
    
    outl = loc[:, :series_length]
    outv = val[:, :series_length]
    (fn, fn2), nt, nb = get_tkernel(series_length, window)
    threshold = cp.float32(threshold * threshold)
    window = cp.int32(window)

    # Launch kernel with batch dimension
    grid = (nb, batch_size, 1)
    block = (nt, 1, 1)

    fn(grid, block, (series_batch.data, outv, outl, window, threshold, series_length))
    fn2(grid, block, (outv, outl, threshold, window))
    
    results = []
    for batch_idx in range(batch_size):
        w = (outl[batch_idx] != -1)
        results.append((outv[batch_idx][w], outl[batch_idx][w]))
    return results

class FastFilter:
    # Compile kernels only once as class variables
    size_kernel = cp.RawKernel(r'''
    extern "C" __global__
    void compute_sizes(const float2* cv, const int* cl, 
                      int* sizes, const int nb,
                      const int snum) {
        int row = blockDim.x * blockIdx.x + threadIdx.x;
        if (row >= snum) return;
        
        int count = 0;
        int row_offset = row * nb;
        for (int col = 0; col < nb; col++) {
            if (cl[row_offset + col] != -1) {
                count++;
            }
        }
        sizes[row] = count;
    }
    ''', 'compute_sizes')

    filter_kernel = cp.RawKernel(r'''
    extern "C" __global__
    void filter_arrays(const float2* cv, const int* cl,
                      float2* cv_out, int* cl_out,
                      const int* sizes, const int* offsets,
                      const int nb, const int snum) {
        int row = blockDim.x * blockIdx.x + threadIdx.x;
        if (row >= snum) return;
        
        int out_idx = offsets[row];
        int row_offset = row * nb;
        for (int col = 0; col < nb; col++) {
            if (cl[row_offset + col] != -1) {
                cv_out[out_idx].x = cv[row_offset + col].x;
                cv_out[out_idx].y = cv[row_offset + col].y;
                cl_out[out_idx] = cl[row_offset + col];
                out_idx++;
            }
        }
    }
    ''', 'filter_arrays')

    def __init__(self, snum, nb):
        self.snum = snum
        self.nb = nb
        self.sizes = cp.empty(snum, dtype=cp.int32)
        self.offsets = cp.empty(snum, dtype=cp.int32)
        
    def filter_arrays(self, cv, cl):
        # Input handling
        cv = cp.asarray(cv, dtype=cp.complex64)
        cl = cp.asarray(cl, dtype=cp.int32)
        
        if cv.ndim == 2:
            if cv.shape != (self.snum, self.nb):
                raise ValueError(f"Expected shape ({self.snum}, {self.nb}), got {cv.shape}")
            cv = cv.reshape(-1)
        
        if cl.ndim == 2:
            if cl.shape != (self.snum, self.nb):
                raise ValueError(f"Expected shape ({self.snum}, {self.nb}), got {cl.shape}")
            cl = cl.reshape(-1)
            
        # Compute sizes
        n_threads = (self.snum + 1023) // 1024
        blocks = (self.snum + n_threads - 1) // n_threads
        
        self.size_kernel((blocks,), (n_threads,),
                        (cv, cl, self.sizes, self.nb, self.snum))

        # Compute offsets
        cp.cumsum(self.sizes[:-1], out=self.offsets[1:])
        self.offsets[0] = 0
        
        # Allocate output arrays
        total_size = int(self.sizes.sum())
        cv_out = cp.empty(total_size, dtype=cp.complex64)
        cl_out = cp.empty(total_size, dtype=cp.int32)
        

        # Run filter kernel
        # This is the slowest part of the function, everything else is negligible
        self.filter_kernel((blocks,), (n_threads,),
                          (cv, cl, cv_out, cl_out, 
                           self.sizes, self.offsets,
                           self.nb, self.snum))
        # Create template mapping array
        template_map = cp.empty(total_size, dtype=cp.int32)
        
        # Split into separate arrays
        cv_filtered = []
        cl_filtered = []
        start = 0
        for idx, size in enumerate(self.sizes):
            size = int(size)
            if size > 0:
                cv_filtered.append(cv_out[start:start + size])
                cl_filtered.append(cl_out[start:start + size])
                template_map[start:start + size] = idx
            else:
                cv_filtered.append(cp.empty(0, dtype=cp.complex64))
                cl_filtered.append(cp.empty(0, dtype=cp.int32))
            start += size
        return cv_filtered, cl_filtered, cv_out, cl_out, self.sizes, template_map

class CUDAThresholdCluster(_BaseThresholdCluster):
    def __init__(self, series_batch, analyse_slice):
        self.series_batch = series_batch
        # This value is hardcoded as it is the longest length currently
        # supported. Memory usage for this is tiny anyway so no need to be
        # shorter.
        self.batch_mem_size = 1024
        self.batch_size, self.series_length = series_batch.shape
        self.analyse_start = cp.int32(analyse_slice.start)
        self.analyse_end = cp.int32(analyse_slice.stop)
        self.analyse_len = self.analyse_end = self.analyse_start

        global val
        global loc
        if val is None or val.size < self.batch_size * self.batch_mem_size:
            val = cp.zeros((self.batch_size, self.batch_mem_size), dtype=cp.complex64)
        if loc is None or loc.size < self.batch_size * self.batch_mem_size:
            loc = cp.zeros((self.batch_size, self.batch_mem_size), cp.int32)

        self.outl = loc
        self.outv = val
        # This is kind of hardcoded here, sorry. We maybe should pass window in here.
        nb = int(cp.ceil(self.analyse_len / float(2048)))
        self.fast_filter = FastFilter(snum=self.batch_size, nb=nb)

    def threshold_and_cluster(self, threshold, window):
        threshold = threshold * threshold
        threshold = cp.asarray(threshold, dtype=cp.float32)
        window = cp.int32(window)

        (fn, fn2), nt, nb = get_tkernel(
            self.series_length,
            self.analyse_len,
            window,
            block_mem_size=self.batch_mem_size,
            batch_size=self.batch_size
        )
        grid = (nb, self.batch_size, 1)
        block = (nt, 1, 1)

        fn(grid, block, (self.series_batch.data, self.outv, self.outl, window, threshold, self.series_length, self.analyse_start))
        fn2(grid, block, (self.outv, self.outl, threshold, window))

        results = []
        # for batch_idx in range(self.batch_size):
        #     cl = self.outl[batch_idx][:nb]  # Clustered locations for this batch
        #     cv = self.outv[batch_idx][:nb]  # Clustered values for this batch
        #     w = (cl != -1)  # Valid locations
        #     results.append((cv[w], cl[w]))

        # slightly faster version
        cv_filtered, cl_filtered, cv_out, cl_out, sizes, template_map = \
            self.fast_filter.filter_arrays(self.outv[:, :nb], self.outl[:, :nb])
        results = list(zip(cv_filtered, cl_filtered))
        return results, cv_out, cl_out, sizes, template_map

def _threshold_cluster_factory(*args, **kwargs):
    return CUDAThresholdCluster

