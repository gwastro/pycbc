# Copyright (C) 2014 Josh Willis
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

cdef extern from "simd_threshold_ccode.cpp":
    void _parallel_threshold(int64_t N, float complex * arr,
                             float complex * outv,
                             uint32_t * outl,
                             uint32_t * count,
                             const float v)
    int _parallel_thresh_cluster(float complex * inarr,
                                 const uint32_t arrlen,
                                 float complex * values,
                                 uint32_t * locs,
                                 const float thresh, const uint32_t winsize,
                                 const uint32_t segsize)

# As a reference for sending numpy arrays onto C++
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
# http://docs.cython.org/en/latest/src/userguide/wrapping_CPlusPlus.html
# https://cython.readthedocs.io/en/latest/src/tutorial/numpy.html
# https://stackoverflow.com/questions/21242160/how-to-build-a-cython-wrapper-for-c-function-with-stl-list-parameter

@cython.boundscheck(False)
@cython.wraparound(False)
def parallel_threshold(int N,
                       numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] arr not None,
                       numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] outv not None,
                       numpy.ndarray[uint32_t, ndim=1, mode="c"] outl not None,
                       numpy.ndarray[uint32_t, ndim=1, mode="c"] count not None,
                       float threshold):
    _parallel_threshold(N, &arr[0], &outv[0], &outl[0], &count[0], threshold)


@cython.boundscheck(False)
@cython.wraparound(False)
def parallel_thresh_cluster(numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] series not None,
                            int slen,
                            numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] values not None,
                            numpy.ndarray[uint32_t, ndim=1, mode="c"] locs not None,
                            float thresh,
                            int window,
                            int segsize):
    return _parallel_thresh_cluster(&series[0], slen, &values[0], &locs[0],
                                    thresh, window, segsize)
