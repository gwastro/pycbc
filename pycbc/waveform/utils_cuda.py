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
from pycbc.types import FrequencySeries
from mako.template import Template
from pycuda.compiler import SourceModule
import numpy

time_shift_kernel = Template("""
__global__ void fseries_ts(float2 *out, float phi,
                           int kmin, int kmax){
    /*
      Input parameters:
      =================

      out:  float2 pointer
            The input frequency series to shift;
            will be shifted in-place

      phi:  float
            Equals -2*pi*delta_f*time_shift

      kmin: int
            minimum index to examine or write

      kmax: int
            maximum index to examine or write

    */

    float x, y;
    int i;
    float2 tmp, htmp;

    i = ${ntpb}*blockIdx.x + threadIdx.x;

    if ((i >= kmin) && (i < kmax)){
       htmp = out[i];
       __sincosf(phi*i, &y, &x);
       tmp.x = x*htmp.x-y*htmp.y;
       tmp.y = x*htmp.y+y*htmp.x;
       out[i] = tmp;
    }

    return;
    }
""")

# Right now, hardcoding the number of threads per block
nt = 1024
nt_float = numpy.float32(nt)
mod = SourceModule(time_shift_kernel.render(ntpb=nt))
fseries_ts_fn = mod.get_function("fseries_ts")
fseries_ts_fn.prepare("Pfii")

def apply_fseries_time_shift(htilde, dt, kmin=0, copy=True):
    """Shifts a frequency domain waveform in time. The waveform is assumed to
    be sampled at equal frequency intervals.
    """
    if htilde.precision != 'single':
        raise NotImplementedError("CUDA version of apply_fseries_time_shift only supports single precision")

    if copy:
        out = htilde.copy()
    else:
        out = htilde

    kmin = numpy.int32(kmin)
    kmax = numpy.int32(len(htilde))
    nb = int(numpy.ceil(kmax / nt_float))
    if nb > 1024:
        raise ValueError("More than 1024 blocks not supported yet")

    phi = numpy.float32(-2 * numpy.pi * dt * htilde.delta_f)
    fseries_ts_fn.prepared_call((nb, 1), (nt, 1, 1), out.data.gpudata, phi, kmin, kmax)
    if copy:
        htilde = FrequencySeries(out, delta_f=htilde.delta_f, epoch=htilde.epoch,
                                 copy=False)
    return htilde
