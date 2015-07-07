# Copyright (C) 2015 Joshua Willis
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

"""
This module defines optimization flags and determines hardware features that some
other modules and packages may use.
"""
import os, sys
import pycbc

# Work around different Python versions to get runtime
# info on hardware cache sizes
_USE_SUBPROCESS = False
HAVE_GETCONF = False
if sys.version_info >= (2, 7):
    import subprocess
    _USE_SUBPROCESS = True
    HAVE_GETCONF = True
else:
    try:
        import commands # Only available on Unix?
        HAVE_GETCONF = True
    except:
        pass

if pycbc.HAVE_OMP:
    omp_support = """
#include <omp.h>
"""
    omp_libs = ['gomp']
    omp_flags = ['-fopenmp']
else:
    omp_support = ""
    omp_libs = []
    omp_flags = []

# The following are intended to also be used by weave
# functions. It would be good to figure out a way to
# test for the presence of <x86intrin.h>, and make
# this conditional, or else define something to
# flag to other codes that it cannot be used.
#
# The constant ALGN is how much to align in bytes, and
# ALGN_FLT and ALGN_DBL are that value as a number of
# floats and doubles, respectively.

simd_intel_intrin_support = """
#include <x86intrin.h>

#ifdef __AVX2__
#define _HAVE_AVX2 1
#else
#define _HAVE_AVX2 0
#endif

#ifdef __AVX__
#define _HAVE_AVX 1
#else
#define _HAVE_AVX 0
#endif

#ifdef __SSE4_1__
#define _HAVE_SSE4_1 1
#else
#define _HAVE_SSE4_1 0
#endif

#ifdef __SSE3__
#define _HAVE_SSE3 1
#else
#define _HAVE_SSE3 0
#endif

#if _HAVE_AVX
#define ALGN 32
#define ALGN_FLT 8
#define ALGN_DBL 4
#else
#define ALGN 16
#define ALGN_FLT 4
#define ALGN_DBL 2
#endif

"""

if HAVE_GETCONF:
    if _USE_SUBPROCESS:
        def getconf(confvar):
            return int(subprocess.check_output(['getconf', confvar]))
    else:
        def getconf(confvar):
            retlist = commands.getstatusoutput('getconf ' + confvar)
            return int(retlist[1])

    LEVEL1_DCACHE_SIZE = getconf('LEVEL1_DCACHE_SIZE')
    LEVEL1_DCACHE_ASSOC = getconf('LEVEL1_DCACHE_ASSOC')
    LEVEL1_DCACHE_LINESIZE = getconf('LEVEL1_DCACHE_LINESIZE')
    LEVEL2_CACHE_SIZE = getconf('LEVEL2_CACHE_SIZE')
    LEVEL2_CACHE_ASSOC = getconf('LEVEL2_CACHE_ASSOC')
    LEVEL2_CACHE_LINESIZE = getconf('LEVEL2_CACHE_LINESIZE')
    LEVEL3_CACHE_SIZE = getconf('LEVEL3_CACHE_SIZE')
    LEVEL3_CACHE_ASSOC = getconf('LEVEL3_CACHE_ASSOC')
    LEVEL3_CACHE_LINESIZE = getconf('LEVEL3_CACHE_LINESIZE')
