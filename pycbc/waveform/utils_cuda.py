# Copyright (C) 2018 Josh Willis
#
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
"""This module contains the CUDA-specific code for
   convenience utilities for manipulating waveforms
"""
from __future__ import absolute_import
from pycbc.types import FrequencySeries
import pycuda.gpuarray
from mako.template import Template
from pycuda.compiler import SourceModule
import numpy

time_shift_kernel = Template("""
__global__ void fseries_ts(float2 *h, float phi,
                           int kmin, int kmax,
                           float2 *out){
    /*
      Input parameters:
      =================

      h:    float2 pointer
            The input frequency series to shift
 
      phi:  float
            Equals -2*pi*delta_f*time_shift

      kmin: int
            minimum index to examine or write

      kmax: int
            maximum index to examine or write

      Output parameters:
      ==================

      out:  float2 pointer
            The output array, may be the same as h

    */

    float x, y;
    unsigned int i;
    float2 tmp;

    i = ${ntpb}*blockIdx.x + threadIdx.x;

    if ((i >= kmin) && (i < kmax)){
       __sincosf(phi*i, &y, &x);
       tmp.x = x*h.x-y*h.y;
       tmp.y = x*h.y+y*h.x;
       h[i] = tmp;
    }
    
    return;
    }
""")

ts_kernel_cache = {}
def get_ts_kernel(nb):
    if nb > 1024:
        raise ValueError("More than 1024 blocks not supported yet")

    try:
        return ts_kernel_cache[nb]
    except KeyError:
        mod = SourceModule(time_shift_kernel.render(ntpb=nt))
        fn = mod.get_function("fseries_ts")
        fn.prepare("PfiiP")
        ts_kernel_cache[hlen] = fn
        return ts_kernel_cache[hlen]

def apply_fseries_time_shift(htilde, dt, kmin=0, copy=True):
    """Shifts a frequency domain waveform in time. The waveform is assumed to
    be sampled at equal frequency intervals.
    """
    if htilde.precision != 'single':
        raise NotImplementedError("CUDA version of apply_fseries_time_shift only supports single precision")

    if copy:
        out = htilde.copy().data.gpudata
    else:
        out = htilde.data.gpudata

    # Right now, hardcoding the number of threads per block
    nt = numpy.int32(1024)
    nb = numpy.int32(numpy.ceil(hlen / 1024.0))
    phi = numpy.float32(-2 * numpy.pi * dt * htilde.delta_f)
    kmax = numpy.int32(len(htilde))
    kmin = numpy.int32(kmin)
    fn = get_ts_kernel(nb).prepared_call
    fn((nb,1), (nt,1,1), htilde.data.gpudata, phi, kmin, kmax, out)
    if copy:
        htilde = FrequencySeries(out, delta_f=htilde.delta_f, epoch=htilde.epoch,
                                 copy=False)
    return htilde
