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
"""This module contains the CPU-specific code for
   convenience utilities for manipulating waveforms
"""
from __future__ import absolute_import
from pycbc import WEAVE_FLAGS
from pycbc.types import FrequencySeries
from pycbc.opt import omp_libs, omp_flags
from pycbc.weave import inline
import numpy


_apply_shift_code = r"""
    #include <math.h>
    // cast the output to a float array for faster processing
    // this takes advantage of the fact that complex arrays store
    // their real and imaginary values next to each other in memory
    double* outptr = (double*) out;
    outptr += 2*kmin; // move to the start position
    double cphi = (double) phi;
    double re_h, im_h;
    int update_interval = 100;
    int jj = update_interval;
    double re_shift, im_shift, re_lastshift, im_lastshift;
    double re_inc = cos(cphi);
    double im_inc = sin(cphi);
    for (int kk=kmin; kk<kmax; kk++){
        if (jj == update_interval) {
            // recompute the added value to reduce numerical error
            re_shift = cos(cphi * (double) kk);
            im_shift = sin(cphi * (double) kk);
            jj = 0;
        }
        re_h = *outptr;
        im_h = *(outptr+1);
        *outptr = re_shift * re_h - im_shift * im_h; // the real part
        *(outptr+1) = re_shift * im_h + im_shift * re_h; // the imag part
        // increase the shift for the next element
        re_lastshift = re_shift;
        im_lastshift = im_shift;
        re_shift = re_lastshift * re_inc - im_lastshift * im_inc;
        im_shift = re_lastshift * im_inc + im_lastshift * re_inc;
        jj += 1;
        outptr += 2; // move to the next element
    }
    """
# for single precision
_apply_shift_code32 = _apply_shift_code.replace('double', 'float')

def apply_fseries_time_shift(htilde, dt, kmin=0, copy=True):
    """Shifts a frequency domain waveform in time. The waveform is assumed to
    be sampled at equal frequency intervals.
    """
    out = numpy.array(htilde.data, copy=copy)
    phi = -2 * numpy.pi * dt * htilde.delta_f # pylint:disable=unused-variable
    kmax = len(htilde) # pylint:disable=unused-variable
    if htilde.precision == 'single':
        code = _apply_shift_code32
    else:
        code = _apply_shift_code
    inline(code, ['out', 'phi', 'kmin', 'kmax'],
           extra_compile_args=[WEAVE_FLAGS]+omp_flags,
           libraries=omp_libs)
    if copy:
        htilde = FrequencySeries(out, delta_f=htilde.delta_f, epoch=htilde.epoch,
                                 copy=False)
    return htilde
