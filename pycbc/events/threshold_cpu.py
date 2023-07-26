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
import numpy
from .simd_threshold_cython import parallel_thresh_cluster, parallel_threshold
from .eventmgr import _BaseThresholdCluster
from .. import opt

if opt.HAVE_GETCONF:
    default_segsize = opt.LEVEL2_CACHE_SIZE / numpy.dtype('complex64').itemsize
else:
    # Seems to work for Sandy Bridge/Ivy Bridge/Haswell, for now?
    default_segsize = 32768

def threshold_numpy(series, value):
    arr = series.data
    locs = numpy.where(arr.real**2 + arr.imag**2 > value**2)[0]
    vals = arr[locs]
    return locs, vals


outl = None
outv = None
count = None
def threshold_inline(series, value):
    arr = numpy.array(series.data, copy=False, dtype=numpy.complex64)
    global outl, outv, count
    if outl is None or len(outl) < len(series):
        outl = numpy.zeros(len(series), dtype=numpy.uint32)
        outv = numpy.zeros(len(series), dtype=numpy.complex64)
        count = numpy.zeros(1, dtype=numpy.uint32)

    N = len(series)
    threshold = value**2.0
    parallel_threshold(N, arr, outv, outl, count, threshold)
    num = count[0]
    if num > 0:
        return outl[0:num], outv[0:num]
    else:
        return numpy.array([], numpy.uint32), numpy.array([], numpy.float32)

# threshold_numpy can also be used here, but for now we use the inline code
# in all instances. Not sure why we're defining threshold *and* threshold_only
# but we are, and I'm not going to change this at this point.
threshold = threshold_inline
threshold_only = threshold_inline

class CPUThresholdCluster(_BaseThresholdCluster):
    def __init__(self, series):
        self.series = numpy.array(series.data, copy=False,
                                  dtype=numpy.complex64)
        self.slen = numpy.uint32(len(series))
        self.outv = numpy.zeros(self.slen, numpy.complex64)
        self.outl = numpy.zeros(self.slen, numpy.uint32)
        self.segsize = numpy.uint32(default_segsize)

    def threshold_and_cluster(self, threshold, window):
        self.count = parallel_thresh_cluster(self.series, self.slen,
                                             self.outv, self.outl,
                                             numpy.float32(threshold),
                                             numpy.uint32(window),
                                             self.segsize)
        if self.count > 0:
            return self.outv[0:self.count], self.outl[0:self.count]
        else:
            return numpy.array([], dtype = numpy.complex64), numpy.array([], dtype = numpy.uint32)


def _threshold_cluster_factory(series):
    return CPUThresholdCluster
