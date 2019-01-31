# Copyright (C) 2016  Josh Willis
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
import mako.template
from pycuda import gpuarray
from pycuda.compiler import SourceModule
import pycbc.scheme
from pycbc.types import zeros

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

kernel_sources = mako.template.Template("""
texture<float, 1> freq_tex;
texture<float, 1> amp_tex;
texture<float, 1> phase_tex;

__device__ int binary_search(float freq, int lower, int upper){

    /*

       Input parameters:
       =================

       freq:  The target frequency

       lower: The index into the frequency texture at which
              to start the search

       upper: The index into the frequency texture at which
              to end the search

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
        float fcomp = tex1Dfetch(freq_tex, (float) mid);
        if (fcomp <= freq){
          begin = mid+1;
        } else {
          end = mid;
        }
    }

    return begin-1;
}


__global__ void find_block_indices(int *lower, int *upper, int texlen,
                                   float df, float flow, float fmax){

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

    lower[i] = binary_search(ffirst, 0, texlen);
    upper[i] = binary_search(flast, 0, texlen) + 1;

    return;
}


__global__ void linear_interp(float2 *h, float df, int hlen,
                              float flow, float fmax, int texlen,
                              int *lower, int *upper){

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

      Global variables:
      ===================

      freq_tex:  Texture of sample frequencies (its length is texlen)

      amp_tex:   Texture of amplitudes corresponding to sample frequencies

      phase_tex: Texture of phases corresponding to sample frequencies


      Output parameters:
      ==================

      h: array of complex

    */

    __shared__ int low[1];
    __shared__ int high[1];
    int idx;
    float2 tmp;
    float amp, freq, phase, inv_df, x, y;
    float a0, a1, f0, f1, p0, p1;

    // Load values in global memory into shared memory that
    // all threads in this block will use:

    if (threadIdx.x == 0) {
        low[0] = lower[blockIdx.x];
        high[0] = upper[blockIdx.x];
    }
    __syncthreads();

    int i = ${ntpb}*blockIdx.x + threadIdx.x;

    if (i < hlen){

        freq = df*i;

        if ( (freq<flow) || (freq>fmax) ){
          tmp.x = 0.0;
          tmp.y = 0.0;
        } else {
          idx = binary_search(freq, low[0], high[0]);
          if (idx < texlen-1) {
              f0 = tex1Dfetch(freq_tex, idx);
              f1 = tex1Dfetch(freq_tex, idx+1);
              inv_df = 1.0/(f1-f0);
              a0 = tex1Dfetch(amp_tex, idx);
              a1 = tex1Dfetch(amp_tex, idx+1);
              p0 = tex1Dfetch(phase_tex, idx);
              p1 = tex1Dfetch(phase_tex, idx+1);
              amp = a0*inv_df*(f1-freq) + a1*inv_df*(freq-f0);
              phase = p0*inv_df*(f1-freq) + p1*inv_df*(freq-f0);
          } else {
             // We must have idx = texlen-1, so this frequency
             // exactly equals fmax
             amp = tex1Dfetch(amp_tex, idx);
             phase = tex1Dfetch(phase_tex, idx);
          }
          __sincosf(phase, &y, &x);
          tmp.x = amp*x;
          tmp.y = amp*y;
        }

       h[i] = tmp;
    }

    return;

}
""")

dckernel_cache = {}
def get_dckernel(slen):
    # Right now, hardcoding the number of threads per block
    nt = 1024
    nb = int(numpy.ceil(slen / 1024.0))

    if nb > 1024:
        raise ValueError("More than 1024 blocks not supported yet")

    try:
        return dckernel_cache[nb]
    except KeyError:
        mod = SourceModule(kernel_sources.render(ntpb=nt, nblocks=nb))
        freq_tex = mod.get_texref("freq_tex")
        amp_tex = mod.get_texref("amp_tex")
        phase_tex = mod.get_texref("phase_tex")
        fn1 = mod.get_function("find_block_indices")
        fn1.prepare("PPifff", texrefs=[freq_tex])
        fn2 = mod.get_function("linear_interp")
        fn2.prepare("PfiffiPP", texrefs=[freq_tex, amp_tex, phase_tex])
        dckernel_cache[nb] = (fn1, fn2, freq_tex, amp_tex, phase_tex, nt, nb)
        return dckernel_cache[nb]

class CUDALinearInterpolate(object):
    def __init__(self, output):
        self.output = output.data.gpudata
        self.df = numpy.float32(output.delta_f)
        self.hlen = numpy.int32(len(output))
        lookups = get_dckernel(self.hlen)
        self.fn1 = lookups[0]
        self.fn2 = lookups[1]
        self.freq_tex = lookups[2]
        self.amp_tex = lookups[3]
        self.phase_tex = lookups[4]
        self.nt = lookups[5]
        self.nb = lookups[6]
        self.lower = zeros(self.nb, dtype=numpy.int32).data.gpudata
        self.upper = zeros(self.nb, dtype=numpy.int32).data.gpudata

    def interpolate(self, flow, freqs, amps, phases):
        flow = numpy.float32(flow)
        texlen = numpy.int32(len(freqs))
        fmax = numpy.float32(freqs[texlen-1])
        freqs_gpu = gpuarray.to_gpu(freqs)
        freqs_gpu.bind_to_texref_ext(self.freq_tex, allow_offset=False)
        amps_gpu = gpuarray.to_gpu(amps)
        amps_gpu.bind_to_texref_ext(self.amp_tex, allow_offset=False)
        phases_gpu = gpuarray.to_gpu(phases)
        phases_gpu.bind_to_texref_ext(self.phase_tex, allow_offset=False)
        fn1 = self.fn1.prepared_call
        fn2 = self.fn2.prepared_call
        fn1((1, 1), (self.nb, 1, 1), self.lower, self.upper, texlen, self.df, flow, fmax)
        fn2((self.nb, 1), (self.nt, 1, 1), self.output, self.df, self.hlen, flow, fmax, texlen, self.lower, self.upper)
        pycbc.scheme.mgr.state.context.synchronize()
        return

def inline_linear_interp(amps, phases, freqs, output, df, flow, imin, start_index):
    # Note that imin and start_index are ignored in the GPU code; they are only
    # needed for CPU.
    if output.precision == 'double':
        raise NotImplementedError("Double precision linear interpolation not currently supported on CUDA scheme")
    flow = numpy.float32(flow)
    texlen = numpy.int32(len(freqs))
    fmax = numpy.float32(freqs[texlen-1])
    hlen = numpy.int32(len(output))
    (fn1, fn2, ftex, atex, ptex, nt, nb) = get_dckernel(hlen)
    freqs_gpu = gpuarray.to_gpu(freqs)
    freqs_gpu.bind_to_texref_ext(ftex, allow_offset=False)
    amps_gpu = gpuarray.to_gpu(amps)
    amps_gpu.bind_to_texref_ext(atex, allow_offset=False)
    phases_gpu = gpuarray.to_gpu(phases)
    phases_gpu.bind_to_texref_ext(ptex, allow_offset=False)
    fn1 = fn1.prepared_call
    fn2 = fn2.prepared_call
    df = numpy.float32(df)
    g_out = output.data.gpudata
    lower = zeros(nb, dtype=numpy.int32).data.gpudata
    upper = zeros(nb, dtype=numpy.int32).data.gpudata
    fn1((1, 1), (nb, 1, 1), lower, upper, texlen, df, flow, fmax)
    fn2((nb, 1), (nt, 1, 1), g_out, df, hlen, flow, fmax, texlen, lower, upper)
    pycbc.scheme.mgr.state.context.synchronize()
    return output
