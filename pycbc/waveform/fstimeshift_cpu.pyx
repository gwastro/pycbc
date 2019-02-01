# Copyright (C) 2019 Collin Capano
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
# cython: embedsignature=True
from __future__ import absolute_import
import numpy
cimport numpy, cython

from libc.math cimport (sin, cos)

# automatic switch for complex types
ctypedef fused COMPLEXTYPE:
    float complex
    double complex

# automatic switch for real types
ctypedef fused REALTYPE:
    float
    double

def fstimeshift(numpy.ndarray [COMPLEXTYPE, ndim=1] freqseries,
                REALTYPE phi,  # this should have same precision as freqseries
                int kmin,
                int kmax):
    cdef REALTYPE re_h
    cdef REALTYPE im_h
    cdef unsigned int update_interval = 100
    cdef unsigned int jj = update_interval

    cdef REALTYPE re_shift
    cdef REALTYPE im_shift
    cdef REALTYPE re_lastshift
    cdef REALTYPE im_lastshift

    cdef REALTYPE re_inc = cos(phi)
    cdef REALTYPE im_inc = sin(phi)

    for kk in range(kmin, kmax):
        if jj == update_interval:
            # recompute the added value to reduce numerical error
            re_shift = cos(phi * <REALTYPE>kk)
            im_shift = sin(phi * <REALTYPE>kk)
            jj = 0
        re_h = freqseries[kk].real
        im_h = freqseries[kk].imag
        freqseries[kk].real = re_shift * re_h - im_shift * im_h
        freqseries[kk].imag = re_shift * im_h + im_shift * re_h
        # increase the shift for the next element
        re_lastshift = re_shift
        im_lastshift = im_shift
        re_shift = re_lastshift * re_inc - im_lastshift * im_inc
        im_shift = re_lastshift * im_inc + im_lastshift * re_inc
        jj += 1
