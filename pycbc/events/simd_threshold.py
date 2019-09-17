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
from pycbc.types import zeros, complex64, float32
import numpy as _np
import pycbc.opt
from .simd_threshold_cython import (max_only_ccode,
                                    windowed_max_ccode,
                                    parallel_thresh_cluster_ccode)

"""
This module defines several C functions that compute
a combined thresholding and time clustering of a complex array using a
multithreaded, SIMD vectoried code.

There are four C functions defined:

max_simd: A single-threaded function that uses SIMD vectorization to compute
          the maximum of a complex time series over a given window.

max_simple: A single-threaded function that does NOT use explicit SIMD vectorization,
            to compute the maximum of a complex time series over a given
            window.

windowed_max: A single threaded (but not vectorized) function that finds the
              locations, norms, and complex values of the maxima in each of
              a set of predefined windows in a given complex array. It does
              this by calling max_simd or max_simple on each window with the
              appropriate parameters.

              At present, this function calls max_simple rather than max_simd because
              of bugs in certain versions of gcc on LDG clusters that can cause
              max_simd to be mis-compiled.

parallel_thresh_cluster: A multithreaded function that finds the maxima in
                         each of a set of contiguous, fixed-length windows
                         by calling windowed_max in parallel. It then sweeps
                         through the results of that (single-threaded) and
                         tests for above threshold, and time-clusters the
                         surviving triggers.

A user calls only the last function; the other three exist to conveniently
compartmentalize SIMD code from OpenMP code.
"""

### Now some actual code that just implements the different
### correlations in a parallelized fashion.

class MaxOnlyObject(object):
    def __init__(self, inarray, verbose=0):
        self.inarr = _np.array(inarray.data, copy=False, dtype=float32)
        self.howmany = _np.zeros(1, dtype=_np.int64)
        self.howmany[0] = len(self.inarr)
        self.nstart = _np.zeros(1, dtype=_np.int64)
        self.nstart[0] = 0
        self.cmplx_mval = zeros(1, dtype = complex64)
        self.mval = _np.array(self.cmplx_mval.data, copy=False, dtype=float32)
        self.norm = _np.zeros(1, dtype=float32)
        self.mloc = _np.zeros(1, dtype=_np.int64)
        self.code = max_only_code
        self.support = thresh_cluster_support
        self.verbose = verbose

    def execute(self):
        inarr = self.inarr
        mval = self.mval
        norm = self.norm
        mloc = self.mloc
        nstart = self.nstart
        howmany = self.howmany
        max_only(inarr, mval, norm, mloc, nstart[0], howmany[0])


class WindowedMaxObject(object):
    def __init__(self, inarray, winsize, verbose=0):
        self.inarr = _np.array(inarray.data, copy=False, dtype=complex64)
        self.arrlen = _np.zeros(1, dtype=_np.int64)
        self.arrlen[0] = len(self.inarr)
        self.len_win = winsize
        nwindows = int( len(self.inarr) / winsize)
        if (nwindows * winsize < len(self.inarr)):
            nwindows = nwindows + 1
        self.nwindows = nwindows
        self.cvals = _np.zeros(self.nwindows, dtype=complex64)
        self.norms = _np.zeros(self.nwindows, dtype=float32)
        self.locs = _np.zeros(self.nwindows, dtype=_np.int64)
        self.winsize = _np.zeros(1, dtype=_np.int64)
        self.winsize[0] = self.len_win
        self.startoffset = _np.zeros(1, dtype=_np.int64)
        self.code = windowed_max_code
        self.support = thresh_cluster_support
        self.verbose = verbose

    def execute(self):
        inarr = self.inarr
        arrlen = self.arrlen
        cvals = self.cvals
        norms = self.norms
        locs = self.locs
        winsize = self.winsize
        startoffset = self.startoffset
        windowed_max(inarr, arrlen[0], cvals, norms, locs, winsize[0],
                     startoffset[0])


if pycbc.opt.HAVE_GETCONF:
    default_segsize = pycbc.opt.LEVEL2_CACHE_SIZE / _np.dtype( _np.complex64).itemsize
else:
    # Seems to work for Sandy Bridge/Ivy Bridge/Haswell, for now?
    default_segsize = 32768

class ThreshClusterObject(object):
    """
    This class takes a complex SNR time series, a real threshold value, and a window size
    (expressed as a number of complex sample points), and optionally a segment size and
    verbosity indicator.

    The execute method returns two numpy arrays: one giving the location of clustered
    maxima above the threshold, and the other the corresponding (complex) values of the
    SNR time series at those clustered maxima.  The expectation is that the memory of
    the SNR time series will be filled new inputs many times, and execute() called
    repeatedly.
    """
    def __init__(self, series, window, segsize = default_segsize, verbose=0):
        self.series = _np.array(series.data, copy=False, dtype=complex64)
        self.slen = len(self.series)
        nwindows = int( self.slen / window)
        if (nwindows * window < self.slen):
            nwindows = nwindows + 1
        self.nwindows = nwindows
        self.values = _np.zeros(self.nwindows, dtype=complex64)
        self.locs = _np.zeros(self.nwindows, dtype=_np.uint32)
        self.window = window
        self.segsize = segsize
        self.code = thresh_cluster_code
        self.support = thresh_cluster_support
        self.verbose = verbose

    def execute(self, thresh):
        series = self.series
        slen = self.slen
        values = self.values
        locs = self.locs
        window = self.window
        segsize = self.segsize
        nthr = parallel_thresh_cluster(series, slen, values, locs,
                                       (float) thresh, window, segsize)
        if nthr > 0:
            return self.values[0:nthr], self.locs[0:nthr]
        else:
            return _np.array([], dtype=complex64),\
                 _np.array([], dtype=_np.uint32)

