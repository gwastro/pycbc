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

from __future__ import absolute_import
import cython
import numpy
cimport numpy
from pycbc.types import zeros, complex64, float32
from libc.stdint cimport int64_t, uint32_t

# This file just provides an interface from python to the C code. I choose to
# keep the bulk of the classes in pure python for ease of profiling, which will
# be important for this code.

# FIXME: I had to remove the restrict keyword here
cdef extern from "simd_threshold_ccode.cpp":
    void _parallel_threshold(int64_t N, float * arr,
                             float complex * outv,
                             uint32_t * outl,
                             uint32_t * count,
                             const float v)
    void _max_simd(float * inarr, float * mval,
                   float * norm, int64_t * mloc,
                   int64_t nstart, int64_t howmany)
    void _windowed_max(float complex * inarr,
                       const int64_t arrlen,
                       float complex * cvals,
                       float * norms,
                       int64_t * locs, const int64_t winsize,
                       const int64_t startoffset)
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
                       numpy.ndarray[float, ndim=1, mode="c"] arr not None,
                       numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] outv not None,
                       numpy.ndarray[uint32_t, ndim=1, mode="c"] outl not None,
                       numpy.ndarray[uint32_t, ndim=1, mode="c"] count not None,
                       float threshold):
    _parallel_threshold(N, &arr[0], &outv[0], &outl[0], &count[0], threshold)

@cython.boundscheck(False)
@cython.wraparound(False)
def max_simd(numpy.ndarray[float, ndim=1, mode="c"] inarr not None,
             numpy.ndarray[float, ndim=1, mode="c"] mval not None,
             numpy.ndarray[float, ndim=1, mode="c"] norm not None,
             numpy.ndarray[int64_t, ndim=1, mode="c"] mloc not None,
             int nstart,
             int howmany):
    # FIXME: It's not clear to me if the last 3 are things that are actually
    #        returns.
    _max_simd(&inarr[0], &mval[0], &norm[0], &mloc[0], nstart, howmany)


@cython.boundscheck(False)
@cython.wraparound(False)
def windowed_max(numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] inarr not None,
                 int arrlen,
                 numpy.ndarray[numpy.complex64_t, ndim=1, mode="c"] cvals not None,
                 numpy.ndarray[float, ndim=1, mode="c"] norms not None,
                 numpy.ndarray[int64_t, ndim=1, mode="c"] locs not None,
                 int winsize,
                 int startoffset):
    # FIXME: As before, are any of these ints actually returns?
    _windowed_max(&inarr[0], arrlen, &cvals[0], &norms[0], &locs[0], winsize, startoffset)


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
