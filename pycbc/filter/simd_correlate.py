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
from pycbc.types import float32, complex64
import numpy as _np
from .. import opt
from .simd_correlate_cython import ccorrf_simd, ccorrf_parallel

"""
This module interfaces to C functions for multiplying
the complex conjugate of one vector by a second vector, writing the output
to a third vector. They do this multi-threaded and with SIMD vectorization.

The code defined here, and the other calling that function,
are imported and used in the CPUCorrelator class defined in
matchedfilter_cpu.py.

Two functions are defined in the 'support' C/Cython module:

ccorrf_simd: Runs on a single core, but vectorized
ccorrf_parallel: Runs multicore, but not explicitly vectorized.
                 Parallelized using OpenMP, and calls ccorrf_simd
"""


def correlate_simd(ht, st, qt):
    htilde = _np.array(ht.data, copy=False, dtype=float32)
    stilde = _np.array(st.data, copy=False, dtype=float32)
    qtilde = _np.array(qt.data, copy=False, dtype=float32)
    arrlen = len(htilde)
    ccorrf_simd(htilde, stilde, qtilde, arrlen)


# We need a segment size (number of complex elements) such that *three* segments
# of that size will fit in the L2 cache. We also want it to be a power of two.
# We are dealing with single-precision complex numbers, which each require 8 bytes.
#
# Our kernel is written to assume a complex correlation of single-precision vectors,
# so that's all we support here.  Note that we are assuming that the correct target
# is that the vectors should fit in L2 cache.  Figuring out cache topology dynamically
# is a harder problem than we attempt to solve here.

if opt.HAVE_GETCONF:
    # Since we need 3 vectors fitting in L2 cache, divide by 3
    # We find the nearest power-of-two that fits, and the length
    # of the single-precision complex array that fits into that size.
    pow2 = int(_np.log(opt.LEVEL2_CACHE_SIZE/3.0)/_np.log(2.0))
    default_segsize = pow(2, pow2)/_np.dtype(_np.complex64).itemsize
else:
    # Seems to work for Sandy Bridge/Ivy Bridge/Haswell, for now?
    default_segsize = 8192

def correlate_parallel(ht, st, qt):
    htilde = _np.array(ht.data, copy=False, dtype=complex64)
    stilde = _np.array(st.data, copy=False, dtype=complex64)
    qtilde = _np.array(qt.data, copy=False, dtype=complex64)
    arrlen = len(htilde)
    segsize = default_segsize
    ccorrf_parallel(htilde, stilde, qtilde, arrlen, segsize)
