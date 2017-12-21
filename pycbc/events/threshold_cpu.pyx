# Copyright (C) 2012  Alex Nitz
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
from __future__ import absolute_import
import numpy
cimport numpy
from pycbc import WEAVE_FLAGS
from pycbc.weave import inline
from .simd_threshold import thresh_cluster_support, default_segsize
from .events import _BaseThresholdCluster
from pycbc.opt import omp_libs, omp_flags

def threshold_numpy(series, value):
    arr = series.data
    locs = numpy.where(arr.real**2 + arr.imag**2 > value**2)[0]
    vals = arr[locs]
    return locs, vals

threshold_only = threshold_numpy

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

def threshold_cython(numpy.ndarray[COMPLEXTYPE, ndim=1] series, float threshold):
    locs = []
    cdef unsigned int xmax = series.shape[0]
    cdef float magsq = 0
    cdef float thresholdsq = threshold * threshold
    for i in range(xmax):
        if series[i].real * series[i].real + series[i].imag * series[i].imag > thresholdsq:
            locs.append(i)

    return numpy.array(locs), series[locs]
    
def threshold_inline(series, value):
    return threshold_cython(series.data, value)

threshold=threshold_inline

class CPUThresholdCluster(_BaseThresholdCluster):
    def __init__(self, series):
        self.series = numpy.array(series.data, copy=False)

        self.slen = len(series)
        self.outv = numpy.zeros(self.slen, numpy.complex64)
        self.outl = numpy.zeros(self.slen, numpy.uint32)
        self.segsize = default_segsize
        self.code = """
             return_val = parallel_thresh_cluster(series, (uint32_t) slen, values, locs,
                                         (float) threshold, (uint32_t) window, (uint32_t) segsize);
              """
        self.support = thresh_cluster_support

    def threshold_and_cluster(self, threshold, window):
        series = self.series # pylint:disable=unused-variable
        slen = self.slen # pylint:disable=unused-variable
        values = self.outv
        locs = self.outl
        segsize = self.segsize # pylint:disable=unused-variable
        self.count = inline(self.code, ['series', 'slen', 'values', 'locs', 'threshold', 'window', 'segsize'],
                            extra_compile_args = [WEAVE_FLAGS] + omp_flags,
                            #extra_compile_args = ['-mno-avx -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -O2 -w'] + omp_flags,
                            #extra_compile_args = ['-msse3 -O3 -w'] + omp_flags,
                            support_code = self.support, libraries = omp_libs,
                            auto_downcast = 1)
        if self.count > 0:
            return values[0:self.count], locs[0:self.count]
        else:
            return numpy.array([], dtype = numpy.complex64), numpy.array([], dtype = numpy.uint32)

def _threshold_cluster_factory(series):
    return CPUThresholdCluster
