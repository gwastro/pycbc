# Copyright (C) 2024  The PyCBC team

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


import cupy as cp
import numpy as np
from mako.template import Template

# The interpolation is the result of the call of two kernels.
#
# The first, find_block_indices(), will find the correct upper
# and lower indices into the frequency texture for each thread
# block in the second kernel. These are placed into global memory,
# as that is the only way to communicate between kernels. The
# indices are found by binary search on the sample frequencies
# texture.
#
# The second kernel, linear_interp, takes these upper and lower
# bounds, the texture of freqency samples, and textures containing
# values of the amplitude and phase at those frequencies, and fills
# an array with the (complex) value of the interpolated waveform.
#
# The three interpolation arrays (node locations, amplitude values,
# and phase values) are stored as 1D textures on the GPU, because many
# threads will need to read them concurrently but never write them, and
# the access pattern of a binary search precludes guaranteeing that
# sequential threads will access sequential memory locations.

kernel_sources = Template("""
#include <cuda_runtime.h>
#include <device_launch_parameters.h>     

__device__ int binary_search(float freq, int lower, int upper, const float* freq_tex){

    /*

       Input parameters:
       =================

       freq:  The target frequency

       lower: The index into the frequency texture at which
              to start the search

       upper: The index into the frequency texture at which
              to end the search
       freq_tex: The frequency values array

       Return value:
       =============
       The largest index into the frequency texture for
       which the value of the texture at that index is less
       than or equal to the target frequency 'freq'.

     */

    int begin = lower;
    int end = upper;

    while (begin != end){
        int mid = (begin + end)/2;
        float fcomp = freq_tex[mid];
        if (fcomp <= freq){
          begin = mid+1;
        } else {
          end = mid;
        }
    }

    return begin - 1;
}


extern "C" __global__ void find_block_indices(
    int *lower, int *upper, int texlen, float df, float flow, const float *freq_tex){

    /*

      Input parameters:
      =================

      texlen: The length of the sample frequency texture

      df:     The difference between successive frequencies in the
              output array

      flow:   The minimum frequency at which to generate an interpolated
              waveform

      Global variable:
      ===================

      freq_tex: Texture of sample frequencies (its length is texlen)

      Output parameters:
      ==================

      lower: array of indices, one per thread block, of the lower
             limit for each block within the frequency arrays.

      upper: array of indices, one per thread block, of the upper
             limit for each block within the frequency arrays.

    */

    // This kernel is launched with only one block; the number of
    // threads will equal the number of blocks in the next kernel.
    int i = threadIdx.x;

    // We want to find the index of the smallest freqency in our
    // texture which is greater than the freqency fmatch below:

    float ffirst = i*df*${ntpb};
    float flast = (i+1)*df*${ntpb}-df;
    if (ffirst < flow){
       ffirst = flow;
    }

    lower[i] = binary_search(ffirst, 0, texlen, freq_tex);
    upper[i] = binary_search(flast, 0, texlen, freq_tex) + 1;

    return;
}


extern "C" __global__ void linear_interp(
    float2 *h, float df, int hlen, float flow, float fmax, int texlen,
    const float *freq_tex, const float *amp_tex, const float *phase_tex,                    
    const int *lower, const int *upper){

    /*

      Input parameters:
      =================

      df:     The difference between successive frequencies in the
              output array

      hlen:   The length of the output array

      flow:   The minimum frequency at which to generate an interpolated
              waveform

      fmax:   The maximum frequency in the sample frequency texture; i.e.,
              freq_tex[texlen-1]

      texlen: The common length of the three sample textures

      lower:  Array that for each thread block stores the index into the
              sample frequency array of the largest sample frequency that
              is less than or equal to the smallest frequency considered
              by that thread block.

      upper:  Array that for each thread block stores the index into the
              sample frequency array of the smallest sample frequency that
              is greater than the next frequency considered *after* that
              thread block.

      freq_tex:  Array of sample frequencies (its length is texlen)

      amp_tex:   Array of amplitudes corresponding to sample frequencies

      phase_tex: Array of phases corresponding to sample frequencies


      Output parameters:
      ==================

      h: array of complex

    */

    __shared__ int low[1];
    __shared__ int high[1];

    if (threadIdx.x == 0) {
        low[0] = lower[blockIdx.x];
        high[0] = upper[blockIdx.x];
    }
    __syncthreads();

    int i = ${ntpb} * blockIdx.x + threadIdx.x;

    if (i < hlen){

        float freq = df*i;
        float2 tmp;

        if ( (freq<flow) || (freq>fmax) ){
          tmp.x = 0.0;
          tmp.y = 0.0;
        } else {
          int idx = binary_search(freq, low[0], high[0], freq_tex);
          float amp, phase, inv_df, x, y;
          float a0, a1, f0, f1, p0, p1;
          
          if (idx < texlen - 1) {
              f0 = freq_tex[idx];
              f1 = freq_tex[idx+1];
              inv_df = 1.0/(f1-f0);
              a0 = amp_tex[idx];
              a1 = amp_tex[idx+1];
              p0 = phase_tex[idx];
              p1 = phase_tex[idx+1];
                          
              amp = a0*inv_df*(f1-freq) + a1*inv_df*(freq-f0);
              phase = p0*inv_df*(f1-freq) + p1*inv_df*(freq-f0);
          } else {
             // We must have idx = texlen-1, so this frequency
             // exactly equals fmax
             amp = amp_tex[idx];
             phase = phase_tex[idx];
          }
          __sincosf(phase, &y, &x);
          tmp.x = amp*x;
          tmp.y = amp*y;
        }

       h[i] = tmp;
    }

}
""")

dckernel_cache = {}
def get_dckernel(slen):
    # Right now, hardcoding the number of threads per block
    nt = 1024
    nb = int(np.ceil(slen / 1024.0))

    if nb > 1024:
        raise ValueError("More than 1024 blocks not supported yet")

    if nb not in dckernel_cache:
        mod = cp.RawModule(code=kernel_sources.render(ntpb=nt))
        fn1 = mod.get_function("find_block_indices")
        fn2 = mod.get_function("linear_interp")
        dckernel_cache[nb] = (fn1, fn2, nt, nb)

    return dckernel_cache[nb]

class CUPYLinearInterpolate(object):
    def __init__(self, output):
        self.output = output.data
        self.df = np.float32(output.delta_f)
        self.hlen = np.int32(len(output))
        lookups = get_dckernel(self.hlen)
        self.fn1 = lookups[0]
        self.fn2 = lookups[1]
        self.nt = lookups[2]
        self.nb = lookups[3]
        self.lower = cp.zeros(self.nb, dtype=np.int32)
        self.upper = cp.zeros(self.nb, dtype=np.int32)

    def interpolate(self, flow, freqs, amps, phases):
        flow = np.float32(flow)
        texlen = np.int32(len(freqs))
        fmax = np.float32(freqs[texlen-1])
        freqs_gpu = cp.asarray(freqs)
        amps_gpu = cp.asarray(amps)
        phases_gpu = cp.asarray(phases)
        self.fn1(
            (1,) , (self.nb,),
            (self.lower, self.upper, texlen, self.df, flow, freqs_gpu))
        self.fn2(
            (self.nb,), (self.nt,),
            (self.output, self.df, self.hlen, flow, fmax, texlen, freqs_gpu, amps_gpu, phases_gpu, self.lower, self.upper)
        )
        return

def inline_linear_interp(amps, phases, freqs, output, df, flow, imin, start_index):
    # Note that imin and start_index are ignored in the GPU code; they are only
    # needed for CPU.
    if output.precision == 'double':
        raise NotImplementedError("Double precision linear interpolation not currently supported on CUDA scheme")
    flow = np.float32(flow)
    texlen = np.int32(len(freqs))
    fmax = np.float32(freqs[texlen-1])
    hlen = np.int32(len(output))
    (fn1, fn2, nt, nb) = get_dckernel(hlen)

    freqs_gpu = cp.asarray(freqs)
    amps_gpu = cp.asarray(amps)
    phases_gpu = cp.asarray(phases)

    df = np.float32(df)
    g_out = output.data
    lower = cp.zeros(nb, dtype=np.int32)
    upper = cp.zeros(nb, dtype=np.int32)
    fn1(
        (1,), (nb,),
        (lower, upper, texlen, df, flow, freqs_gpu)
    )
    fn2(
        (nb,), (nt,),
        (g_out, df, hlen, flow, fmax, texlen, freqs_gpu, amps_gpu, phases_gpu, lower, upper)
    )
    return output
