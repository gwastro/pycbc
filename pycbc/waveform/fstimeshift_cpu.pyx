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

# for double precision
def fstimeshift(numpy.ndarray [double complex, ndim=1] freqseries,
                double phi,
                int kmin,
                int kmax):
    cdef double re_h
    cdef double im_h
    cdef unsigned int update_interval = 100
    cdef unsigned int jj = update_interval

    cdef double re_shift
    cdef double im_shift
    cdef double re_lastshift
    cdef double im_lastshift

    cdef double re_inc = cos(phi)
    cdef double im_inc = sin(phi)

    for kk in range(kmin, kmax):
        if jj == update_interval:
            # recompute the added value to reduce numerical error
            re_shift = cos(phi * <double>kk)
            im_shift = sin(phi * <double>kk)
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


# for single precision
def fstimeshift32(numpy.ndarray [float complex, ndim=1] freqseries,
                float phi,
                int kmin,
                int kmax):
    cdef float re_h
    cdef float im_h
    cdef unsigned int update_interval = 100
    cdef unsigned int jj = update_interval

    cdef float re_shift
    cdef float im_shift
    cdef float re_lastshift
    cdef float im_lastshift

    cdef float re_inc = cos(phi)
    cdef float im_inc = sin(phi)

    for kk in range(kmin, kmax):
        if jj == update_interval:
            # recompute the added value to reduce numerical error
            re_shift = cos(phi * <float>kk)
            im_shift = sin(phi * <float>kk)
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
