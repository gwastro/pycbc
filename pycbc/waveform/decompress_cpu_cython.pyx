# Copyright (C) 2019 Ian Harry
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

import cython
import numpy
cimport numpy
from pycbc.types import zeros, complex64, float32
from libc.stdint cimport int64_t, uint32_t

# This file just provides an interface from python to the C code. I choose to
# keep the bulk of the classes in pure python for ease of profiling, which will
# be important for this code.

cdef extern from "decompress_cpu_ccode.cpp":
    void _decomp_ccode_double(double complex * h,
                              double delta_f,
                              const int64_t hlen,
                              const int64_t start_index,
                              double * sample_frequencies,
                              double * amp,
                              double * phase,
                              const int64_t sflen,
                              const int64_t imin)
    void _decomp_ccode_float(float complex * h,
                             float delta_f,
                             const int64_t hlen,
                             const int64_t start_index,
                             float * sample_frequencies,
                             float * amp,
                             float * phase,
                             const int64_t sflen,
                             const int64_t imin)

# See simd_threshold_cython in events module for some guidance for how I
# constructed this in this way

@cython.boundscheck(False)
@cython.wraparound(False)
def decomp_ccode_double(numpy.ndarray[numpy.complex128_t, ndim=1, mode="c"] h not None,
                        double delta_f,
                        int hlen,
                        int start_index,
                        numpy.ndarray[double, ndim=1, mode="c"] sample_frequencies not None,
                        numpy.ndarray[double, ndim=1, mode="c"] amp not None,
                        numpy.ndarray[double, ndim=1, mode="c"] phase not None,
                        int sflen,
                        int imin):
    _decomp_ccode_double(&h[0], delta_f, hlen, start_index,
                         &sample_frequencies[0], &amp[0], &phase[0],
                         sflen, imin)

@cython.boundscheck(False)
@cython.wraparound(False)
def decomp_ccode_float(numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] h not None,
                       float delta_f,
                       int hlen,
                       int start_index,
                       numpy.ndarray[float, ndim=1, mode="c"] sample_frequencies not None,
                       numpy.ndarray[float, ndim=1, mode="c"] amp not None,
                       numpy.ndarray[float, ndim=1, mode="c"] phase not None,
                       int sflen,
                       int imin):
    _decomp_ccode_float(&h[0], delta_f, hlen, start_index,
                        &sample_frequencies[0], &amp[0], &phase[0],
                        sflen, imin)

