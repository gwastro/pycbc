# Copyright (C) 2015  Alex Nitz, Josh Willis
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

import functools
import numpy
import cupy as cp
import lal
from mako.template import Template

LALARGS = {
    'TWOPI': lal.TWOPI,
}

accum_diff_sq_kernel = cp.ElementwiseKernel(
    "X input",
    "raw Y output",
    "output[i] += norm(input)",
    "accum_diff_sq_kernel"
)

def chisq_accum_bin(chisq, q):
    accum_diff_sq_kernel(q.data, chisq.data)


chisqkernel = Template("""
#include <cstdint>
extern "C" __global__ void power_chisq_at_points_${NP}(
                                      %if fuse:
                                          float2* htilde,
                                          float2* stilde,
                                      %else:
                                          float2* corr,
                                      %endif
                                      float2* outc, unsigned int N,
                                      %for p in range(NP):
                                        float phase${p},
                                      %endfor
                                      uint32_t* kmin,
                                      uint32_t* kmax,
                                      uint32_t* bv,
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

        %if fuse:
            float2 qt, st, ht;
            st = stilde[i];
            ht = htilde[i];
            qt.x = ht.x * st.x + ht.y * st.y;
            qt.y = ht.x * st.y - ht.y * st.x;
        %else:
            float2 qt = corr[i];
        %endif

        %for p in range(NP):
            sincosf(phase${p} * i, &im, &re);
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

chisqkernel_pow2 = Template("""
#include <cstdint>
extern "C" __global__ void power_chisq_at_points_${NP}_pow2(
                                      %if fuse:
                                          float2* htilde,
                                          float2* stilde,
                                      %else:
                                          float2* corr,
                                      %endif
                                      float2* outc, unsigned int N,
                                      %for p in range(NP):
                                        unsigned int points${p},
                                      %endfor
                                      uint32_t* kmin,
                                      uint32_t* kmax,
                                      uint32_t* bv,
                                      unsigned int nbins){
    __shared__ unsigned int s;
    __shared__ unsigned int e;
    __shared__ float2 chisq[${NT} * ${NP}];
    float twopi = ${TWOPI};
    unsigned long long NN;

    NN = (unsigned long long) N;

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

        %if fuse:
            float2 qt, st, ht;
            st = stilde[i];
            ht = htilde[i];
            qt.x = ht.x * st.x + ht.y * st.y;
            qt.y = ht.x * st.y - ht.y * st.x;
        %else:
            float2 qt = corr[i];
        %endif

        %for p in range(NP):
            unsigned long long prod${p} = points${p} * i;
            unsigned int k${p} = (unsigned int) (prod${p}&(NN-1));
            float phase${p} = twopi * k${p}/((float) N);
            __sincosf(phase${p}, &im, &re);
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

@functools.lru_cache(maxsize=None)
def get_pchisq_fn(np, fuse_correlate=False):
    nt = 256
    fn = cp.RawKernel(
        chisqkernel.render(NT=nt, NP=np, fuse=fuse_correlate, **LALARGS),
        f'power_chisq_at_points_{np}',
        backend='nvcc'
    )
    return fn, nt


@functools.lru_cache(maxsize=None)
def get_pchisq_fn_pow2(np, fuse_correlate=False):
    nt = 256
    fn = cp.RawKernel(
        chisqkernel_pow2.render(NT=nt, NP=np, fuse=fuse_correlate, **LALARGS),
        f'power_chisq_at_points_{np}_pow2',
        backend='nvcc'
    )
    return fn, nt

def get_cached_bin_layout(bins):
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
    bv = cp.array(bv, dtype=cp.uint32)
    kmin = cp.array(kmin, dtype=cp.uint32)
    kmax = cp.array(kmax, dtype=cp.uint32)
    return kmin, kmax, bv

def shift_sum_points(num, N, arg_tuple):
    #fuse = 'fuse' in corr.gpu_callback_method
    fuse = False

    fn, nt = get_pchisq_fn(num, fuse_correlate = fuse)
    corr, outp, phase, np, nb, N, kmin, kmax, bv, nbins = arg_tuple
    if fuse:
        args = [corr.htilde.data, corr.stilde.data]
    else:
        args = [corr.data]
    args += [outp, N] + phase[0:num]
    args += [kmin, kmax, bv, nbins]
    fn(
        (nb,),
        (nt,),
        *args,
    )
    
    outp = outp[num*nbins:]
    phase = phase[num:]
    np -= num
    return outp, phase, np

def shift_sum_points_pow2(num, arg_tuple):
    #fuse = 'fuse' in corr.gpu_callback_method
    fuse = False

    fn, nt = get_pchisq_fn_pow2(num, fuse_correlate = fuse)

    corr, outp, points, np, nb, N, kmin, kmax, bv, nbins = arg_tuple
    if fuse:
        args = [corr.htilde.data, corr.stilde.data]
    else:
        args = [corr.data]
    args += [outp, N] + points[0:num] + [kmin, kmax, bv, nbins]
    fn(
        (nb,),
        (nt,),
        tuple(args)
    )
            
    outp = outp[num*nbins:]
    points = points[num:]
    np -= num
    return outp, points, np

@functools.lru_cache(maxsize=None)
def get_cached_pow2(N):
    return not(N & (N-1))

def shift_sum(corr, points, bins):
    kmin, kmax, bv = get_cached_bin_layout(bins)
    nb = len(kmin)
    N = cp.uint32(len(corr))
    is_pow2 = get_cached_pow2(N)
    nbins = cp.uint32(len(bins) - 1)
    outc = cp.zeros((len(points), nbins), dtype=cp.complex64)
    outp = outc.reshape(nbins * len(points))
    np = len(points)

    if is_pow2:
        lpoints = points.tolist()
        while np > 0:
            cargs = (corr, outp, lpoints, np, nb, N, kmin, kmax, bv, nbins)

            if np >= 4:
                outp, lpoints, np = shift_sum_points_pow2(4, cargs)
            elif np >= 3:
                outp, lpoints, np = shift_sum_points_pow2(3, cargs)
            elif np >= 2:
                outp, lpoints, np = shift_sum_points_pow2(2, cargs)
            elif np == 1:
                outp, lpoints, np = shift_sum_points_pow2(1, cargs)
    else:
        phase = [cp.float32(p * 2.0 * cp.pi / N) for p in points]
        while np > 0:
            cargs = (corr, outp, phase, np, nb, N, kmin, kmax, bv, nbins)

            if np >= 4:
                outp, phase, np = shift_sum_points(4, cargs) # pylint:disable=no-value-for-parameter
            elif np >= 3:
                outp, phase, np = shift_sum_points(3, cargs) # pylint:disable=no-value-for-parameter
            elif np >= 2:
                outp, phase, np = shift_sum_points(2, cargs) # pylint:disable=no-value-for-parameter
            elif np == 1:
                outp, phase, np = shift_sum_points(1, cargs) # pylint:disable=no-value-for-parameter

    return (outc.conj() * outc).sum(axis=1).real


# Batched chisq implementation - OPTIMIZED for coalescing and reduced atomics
chisqkernel_pow2_batch = Template("""
#include <cstdint>
extern "C" __global__ void power_chisq_at_points_pow2_batch(
    float2* corr, // 2D: num_templates X N
    float2* outc, // 1D: num_points X max_bin_length 
    unsigned int N, // Scalar
    uint32_t* points, // 1D: num_points
    uint32_t* kmin, // 2D: num_templates X max_bin_size
    uint32_t* kmax, // 2D: num_templates X max_bin_size
    uint32_t* bv, // 2D: num_templates X max_bin_size
    uint32_t* nbins, // 1D: num_templates (per-template bin counts)
    uint32_t* mapping, // 1D: num_points
    uint32_t* nb, // 1D: num_templates
    unsigned int max_nbins, // scalar: max number of bins across templates
    unsigned int num_points, // scalar
    unsigned int num_templates, // scalar
    unsigned int max_bin_size // scalar
)
{
    const unsigned int pnum = blockIdx.y;
    const unsigned int binnum = blockIdx.x;
    const float twopi = ${TWOPI};
    const unsigned long long NN = (unsigned long long) N;
    
    // Early exit check
    if (pnum >= num_points) return;
    
    const unsigned int tempnum = mapping[pnum];
    if (binnum >= nb[tempnum]) return;
    
    // Load bin parameters (all threads load same values, gets cached in L1)
    const unsigned int idx_base = tempnum * max_bin_size + binnum;
    const unsigned int s = kmin[idx_base];
    const unsigned int e = kmax[idx_base];
    
    if (s >= e) return;  // Empty bin
    
    const unsigned int bin_idx = bv[idx_base];
    const unsigned long long point = (unsigned long long) points[pnum];
    const unsigned int corr_base = tempnum * N;
    
    // Each thread accumulates independently (no shared memory needed initially)
    float accum_x = 0.0f;
    float accum_y = 0.0f;
    
    // Main loop - fully coalesced memory access
    const float phase_mult = twopi / ((float) N);
    #pragma unroll 4
    for (unsigned int i = s + threadIdx.x; i < e; i += blockDim.x){
        // Coalesced load
        const float2 qt = corr[corr_base + i];
        
        // Compute phase
        const unsigned int k = (unsigned int)((point * (unsigned long long)i) & (NN-1));
        float re, im;
        __sincosf(phase_mult * k, &im, &re);
        
        // Accumulate
        accum_x += re * qt.x - im * qt.y;
        accum_y += im * qt.x + re * qt.y;
    }
    
    // Use warp shuffle reduction for first stage (no shared memory)
    #pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1) {
        accum_x += __shfl_down_sync(0xffffffff, accum_x, offset);
        accum_y += __shfl_down_sync(0xffffffff, accum_y, offset);
    }
    
    // Shared memory only for cross-warp reduction
    __shared__ float2 warp_sums[16];  // Up to 512 threads = 16 warps
    const unsigned int warp_id = threadIdx.x / 32;
    const unsigned int lane_id = threadIdx.x % 32;
    
    if (lane_id == 0) {
        warp_sums[warp_id].x = accum_x;
        warp_sums[warp_id].y = accum_y;
    }
    __syncthreads();
    
    // Final reduction by first warp
    if (warp_id == 0) {
        float2 sum;
        if (lane_id < (${NT} / 32)) {
            sum = warp_sums[lane_id];
        } else {
            sum.x = 0.0f;
            sum.y = 0.0f;
        }
        
        #pragma unroll
        for (int offset = 8; offset > 0; offset >>= 1) {
            sum.x += __shfl_down_sync(0xffffffff, sum.x, offset);
            sum.y += __shfl_down_sync(0xffffffff, sum.y, offset);
        }
        
        if (lane_id == 0) {
            const unsigned int out_idx = pnum * max_nbins + bin_idx;
            atomicAdd(&outc[out_idx].x, sum.x);
            atomicAdd(&outc[out_idx].y, sum.y);
        }
    }
}
""")


@functools.lru_cache(maxsize=None)
def get_pchisq_fn_pow2_batch():
    nt = 512
    fn = cp.RawKernel(
        chisqkernel_pow2_batch.render(NT=nt, **LALARGS),
        f'power_chisq_at_points_pow2_batch',
        backend='nvcc',
        options=('--use_fast_math',)  # Enable fast math for better performance
    )
    return fn, nt


# Batched chisq for non-power-of-2 FFT lengths
chisqkernel_batch = Template("""
#include <cstdint>
extern "C" __global__ void power_chisq_at_points_batch(
    float2* corr, // 2D: num_templates X N
    float2* outc, // 1D: num_points X max_bin_length 
    unsigned int N, // Scalar
    float* phases, // 1D: num_points (phase multiplier for each point)
    uint32_t* kmin, // 2D: num_templates X max_bin_size
    uint32_t* kmax, // 2D: num_templates X max_bin_size
    uint32_t* bv, // 2D: num_templates X max_bin_size
    uint32_t* nbins, // 1D: num_templates (per-template bin counts)
    uint32_t* mapping, // 1D: num_points
    uint32_t* nb, // 1D: num_templates
    unsigned int max_nbins, // scalar: max number of bins across templates
    unsigned int num_points, // scalar
    unsigned int num_templates, // scalar
    unsigned int max_bin_size // scalar
)
{
    const unsigned int pnum = blockIdx.y;
    const unsigned int binnum = blockIdx.x;
    
    // Early exit check
    if (pnum >= num_points) return;
    
    const unsigned int tempnum = mapping[pnum];
    if (binnum >= nb[tempnum]) return;
    
    // Load bin parameters
    const unsigned int idx_base = tempnum * max_bin_size + binnum;
    const unsigned int s = kmin[idx_base];
    const unsigned int e = kmax[idx_base];
    
    if (s >= e) return;  // Empty bin
    
    const unsigned int bin_idx = bv[idx_base];
    const float phase = phases[pnum];
    const unsigned int corr_base = tempnum * N;
    
    // Each thread accumulates independently
    float accum_x = 0.0f;
    float accum_y = 0.0f;
    
    // Main loop - coalesced memory access
    #pragma unroll 4
    for (unsigned int i = s + threadIdx.x; i < e; i += blockDim.x){
        const float2 qt = corr[corr_base + i];
        
        // Compute phase using sincosf for non-power-of-2
        float re, im;
        sincosf(phase * i, &im, &re);
        
        // Accumulate
        accum_x += re * qt.x - im * qt.y;
        accum_y += im * qt.x + re * qt.y;
    }
    
    // Warp shuffle reduction
    #pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1) {
        accum_x += __shfl_down_sync(0xffffffff, accum_x, offset);
        accum_y += __shfl_down_sync(0xffffffff, accum_y, offset);
    }
    
    // Shared memory for cross-warp reduction
    __shared__ float2 warp_sums[16];
    const unsigned int warp_id = threadIdx.x / 32;
    const unsigned int lane_id = threadIdx.x % 32;
    
    if (lane_id == 0) {
        warp_sums[warp_id].x = accum_x;
        warp_sums[warp_id].y = accum_y;
    }
    __syncthreads();
    
    // Final reduction by first warp
    if (warp_id == 0) {
        float2 sum;
        if (lane_id < (${NT} / 32)) {
            sum = warp_sums[lane_id];
        } else {
            sum.x = 0.0f;
            sum.y = 0.0f;
        }
        
        #pragma unroll
        for (int offset = 8; offset > 0; offset >>= 1) {
            sum.x += __shfl_down_sync(0xffffffff, sum.x, offset);
            sum.y += __shfl_down_sync(0xffffffff, sum.y, offset);
        }
        
        if (lane_id == 0) {
            const unsigned int out_idx = pnum * max_nbins + bin_idx;
            atomicAdd(&outc[out_idx].x, sum.x);
            atomicAdd(&outc[out_idx].y, sum.y);
        }
    }
}
""")


@functools.lru_cache(maxsize=None)
def get_pchisq_fn_batch():
    nt = 512
    fn = cp.RawKernel(
        chisqkernel_batch.render(NT=nt, **LALARGS),
        f'power_chisq_at_points_batch',
        backend='nvcc',
        options=('--use_fast_math',)
    )
    return fn, nt


_bin_layout_cache = {}

def get_cached_bin_layout_batch(bins, bin_lengths):
    """Get or compute cached bin layout arrays.
    
    Returns pre-allocated GPU arrays bv, kmin, kmax, nb that are cached
    based on the bin configuration.
    """
    import time
    t_cache_key_start = time.time()
    # Create cache key from bins - OPTIMIZED
    # Instead of converting to nested tuples, use array hashes which is much faster
    # Hash each bin array and combine them
    import cupy as cp
    cache_key = tuple(
        (len(b), int(b[0]) if len(b) > 0 else 0, int(b[-1]) if len(b) > 0 else 0, 
         hash(b.data.tobytes()) if hasattr(b, 'data') else hash(bytes(b)))
        for b in bins
    )
    t_cache_key = time.time() - t_cache_key_start
    
    if cache_key in _bin_layout_cache:
        print(f"CACHE HIT: Reusing bin layout for {len(bins)} templates (cache_key: {t_cache_key:.4f}s)", flush=True)
        return _bin_layout_cache[cache_key]
    
    print(f"CACHE MISS: Computing bin layout for {len(bins)} templates (cache_key: {t_cache_key:.4f}s)", flush=True)
    
    t_alloc_start = time.time()
    # Convert bins list to batch-friendly format
    # OPTIMIZED: Single pass - build ranges and track max values simultaneously
    BS = 4096
    
    # Pre-compute all bin_ranges for all templates
    all_bin_ranges = []
    max_num_ranges = 0
    max_nbins = 0
    
    for cbins in bins:
        bin_ranges = []
        for i in range(len(cbins)-1):
            s, e = int(cbins[i]), int(cbins[i+1])
            if (e - s) < BS:
                bin_ranges.append((i, s, e))
            else:
                # Calculate chunks without creating arrays
                num_chunks = ((e - s) + BS//2 - 1) // (BS//2)
                for j in range(num_chunks):
                    chunk_s = s + j * (BS//2)
                    chunk_e = min(s + (j + 1) * (BS//2), e)
                    bin_ranges.append((i, chunk_s, chunk_e))
        all_bin_ranges.append(bin_ranges)
        max_num_ranges = max(max_num_ranges, len(bin_ranges))
        max_nbins = max(max_nbins, len(cbins) - 1)
    
    max_bin_size = max_num_ranges
    
    bv = cp.zeros([len(bins), max_bin_size], dtype=cp.uint32)
    kmin = cp.zeros([len(bins), max_bin_size], dtype=cp.uint32)
    kmax = cp.zeros([len(bins), max_bin_size], dtype=cp.uint32)
    nb = cp.zeros(len(bins), dtype=cp.uint32)
    t_alloc = time.time() - t_alloc_start
    
    t_loop_start = time.time()
    for idx1, bin_ranges in enumerate(all_bin_ranges):
        # Bulk copy to GPU (much faster than individual assignments)
        if bin_ranges:
            bin_arr = cp.array(bin_ranges, dtype=cp.uint32)
            bv[idx1, :len(bin_ranges)] = bin_arr[:, 0]
            kmin[idx1, :len(bin_ranges)] = bin_arr[:, 1]
            kmax[idx1, :len(bin_ranges)] = bin_arr[:, 2]
        
        nb[idx1] = len(bin_ranges)
    t_loop = time.time() - t_loop_start
    
    print(f"  Bin layout breakdown - alloc: {t_alloc:.4f}s, loop: {t_loop:.4f}s", flush=True)
    
    # No need to compute max values - we already have them from the first pass
    max_nb_val = max_num_ranges
    max_nbins_val = max_nbins
    
    print(f"  Max values from CPU: max_nb={max_nb_val}, max_nbins={max_nbins_val}", flush=True)
    
    # Cache the result with precomputed max values
    result = (bv, kmin, kmax, nb, max_nb_val, max_nbins_val)
    _bin_layout_cache[cache_key] = result
    return result


def shift_sum_batch(corr, points, bins, bin_lengths, mapping):
    """Compute chisq shift-sum for batched templates.
    
    Parameters
    ----------
    corr : cupy array
        2D array of correlation data (num_templates x frequency_length)
    points : cupy array
        1D array of time indices for triggers
    bins : list of cupy arrays
        List of bin edges for each template
    bin_lengths : cupy array
        Number of bins for each template
    mapping : cupy array
        Maps each point index to its template index
    
    Returns
    -------
    cupy array
        Chisq values for each point
    """
    import time
    t_layout_start = time.time()
    # Get cached bin layout (avoids recomputation) - now includes precomputed max values
    bv, kmin, kmax, nb, max_nb_val, max_nbins_val = get_cached_bin_layout_batch(bins, bin_lengths)
    t_layout = time.time() - t_layout_start
    
    t_setup_start = time.time()
    N = cp.uint32(len(corr[0]))
    is_pow2 = get_cached_pow2(int(N))
    nbins = bin_lengths - 1
    
    # Allocate output with maximum possible size (uniform row width)
    # Use precomputed max_nbins_val to avoid GPU->CPU sync
    outc = cp.zeros((len(points), max_nbins_val), dtype=cp.complex64)
    t_setup = time.time() - t_setup_start

    t_kernel_start = time.time()
    if is_pow2:
        fn, nt = get_pchisq_fn_pow2_batch()
        # Flatten corr array if needed
        if corr.ndim == 1:
            # Single template case - reshape
            corr_flat = corr
        else:
            corr_flat = corr.reshape(-1)
        
        max_bin_size = bv.shape[1]
        args = (corr_flat, outc.reshape(-1), N, points, kmin.reshape(-1), 
            kmax.reshape(-1), bv.reshape(-1), nbins, mapping, nb,
            cp.uint32(max_nbins_val), cp.uint32(len(points)), cp.uint32(len(bins)), cp.uint32(max_bin_size))
        
        # Use precomputed max_nb_val to avoid GPU->CPU sync
        grid_size = (max_nb_val, len(points))
        fn(grid_size, (nt,), args)
    else:
        # Non-power-of-2 case using sincosf
        fn, nt = get_pchisq_fn_batch()
        
        # Flatten corr array if needed
        if corr.ndim == 1:
            corr_flat = corr
        else:
            corr_flat = corr.reshape(-1)
        
        # Calculate phase multipliers for each point
        phases = cp.float32(2 * cp.pi / float(N)) * points.astype(cp.float32)
        
        max_bin_size = bv.shape[1]
        args = (corr_flat, outc.reshape(-1), N, phases, kmin.reshape(-1),
            kmax.reshape(-1), bv.reshape(-1), nbins, mapping, nb,
            cp.uint32(max_nbins_val), cp.uint32(len(points)), cp.uint32(len(bins)), cp.uint32(max_bin_size))
        
        # Use precomputed max_nb_val to avoid GPU->CPU sync
        grid_size = (max_nb_val, len(points))
        fn(grid_size, (nt,), args)
    
    cp.cuda.Stream.null.synchronize()
    t_kernel = time.time() - t_kernel_start
    
    t_reduce_start = time.time()
    result = (outc.conj() * outc).sum(axis=1).real
    cp.cuda.Stream.null.synchronize()
    t_reduce = time.time() - t_reduce_start
    
    print(f"    shift_sum_batch - layout: {t_layout:.4f}s, setup: {t_setup:.4f}s, kernel: {t_kernel:.4f}s, reduce: {t_reduce:.4f}s, triggers: {len(points)}", flush=True)
    
    return result
