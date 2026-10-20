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

val = None
loc = None

# https://stackoverflow.com/questions/77798014/cupy-rawkernel-cuda-error-not-found-named-symbol-not-found-cupy

tkernel1 = mako.template.Template("""
extern "C" __global__ void threshold_and_cluster(float2* in, float2* outv, int* outl, int window, float threshold){
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
extern "C" __global__ void threshold_and_cluster2(float2* outv, int* outl, float threshold, int window){
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

@functools.lru_cache(maxsize=None)
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

    nb = int(cp.ceil(slen / float(window)))

    if nb > 1024:
        raise ValueError("More than 1024 blocks not supported yet")

    fn = cp.RawKernel(
        tkernel1.render(chunk=nt),
        'threshold_and_cluster',
        backend='nvcc'
    )
    fn2 = cp.RawKernel(
        tkernel2.render(blocks=nb),
        'threshold_and_cluster2',
        backend='nvcc'
    )
    return (fn, fn2), nt, nb

def threshold_and_cluster(series, threshold, window):
    global val
    global loc
    if val is None:
        val = cp.zeros(4096*256, dtype=cp.complex64)
    if loc is None:
        loc = cp.zeros(4096*256, cp.int32)

    outl = loc
    outv = val
    slen = len(series)
    series = series.data
    (fn, fn2), nt, nb = get_tkernel(slen, window)
    threshold = cp.float32(threshold * threshold)
    window = cp.int32(window)

    cl = loc[0:nb]
    cv = val[0:nb]

    fn((nb,), (nt,), (series.data, outv, outl, window, threshold))
    fn2((1,), (nb,), (outv, outl, threshold, window))
    w = (cl != -1)
    return cv[w], cl[w]

class CUDAThresholdCluster(_BaseThresholdCluster):
    def __init__(self, series):
        self.series = series

        global val
        global loc
        if val is None:
            val = cp.zeros(4096*256, dtype=cp.complex64)
        if loc is None:
            loc = cp.zeros(4096*256, cp.int32)

        self.outl = loc
        self.outv = val
        self.slen = len(series)

    def threshold_and_cluster(self, threshold, window):
        threshold = cp.float32(threshold * threshold)
        window = cp.int32(window)

        (fn, fn2), nt, nb = get_tkernel(self.slen, window)
        cl = loc[0:nb]
        cv = val[0:nb]

        fn(
            (nt, 1, 1),
            (nb, 1),
            (self.series.data, self.outv, self.outl, window, threshold)
        )
        fn2(
            (nb, 1, 1),
            (1, 1),
            (self.outv, self.outl, threshold, window)
        )
        w = (cl != -1)
        return cp.asnumpy(cv[w]), cp.asnumpy(cl[w])

def _threshold_cluster_factory(series):
    return CUDAThresholdCluster

