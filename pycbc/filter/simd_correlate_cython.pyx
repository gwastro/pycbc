# Copyright (C) 2014-2019 Josh Willis, Ian Harry
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

cdef extern from "simd_correlate_ccode.cpp":
    void _ccorrf_simd(float * inconj, float * innoconj,
                     float * out, const int64_t len)
    void _ccorrf_parallel(float complex * inconj,
                         float complex * innoconj,
                         float complex * out,
                         const int64_t arrlen, const int64_t segsize)

# See simd_threshold_cython in events module for some guidance for how I
# constructed this in this way

@cython.boundscheck(False)
@cython.wraparound(False)
def ccorrf_simd(numpy.ndarray[float, ndim=1, mode="c"] inconj not None,
                numpy.ndarray[float, ndim=1, mode="c"] innoconj not None,
                numpy.ndarray[float, ndim=1, mode="c"] out not None,
                int len):
    _ccorrf_simd(&inconj[0], &innoconj[0], &out[0], len)


@cython.boundscheck(False)
@cython.wraparound(False)
def ccorrf_parallel(numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] inconj not None,
                    numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] innoconj not None,
                    numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] out not None,
                    int arrlen, int segsize):
    _ccorrf_parallel(&inconj[0], &innoconj[0], &out[0], arrlen, segsize)
