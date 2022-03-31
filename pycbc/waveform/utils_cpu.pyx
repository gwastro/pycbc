# Copyright (C) 2018 Collin Capano, Josh Willis
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""This module contains the CPU-specific code for
   convenience utilities for manipulating waveforms
"""
from pycbc.types import FrequencySeries
import numpy

# cython: embedsignature=True
cimport numpy, cython

from libc.math cimport (sin, cos)


def apply_fseries_time_shift(htilde, dt, kmin=0, copy=True):
    """Shifts a frequency domain waveform in time. The waveform is assumed to
    be sampled at equal frequency intervals.
    """
    out = numpy.array(htilde.data, copy=copy)
    phi = -2 * numpy.pi * dt * htilde.delta_f
    kmax = len(htilde)
    # make phi have the same precision as htilde
    if htilde.precision == 'single':
        phi = numpy.float32(phi)
        fstimeshift32(out, phi, kmin, kmax)
    else:
        fstimeshift(out, phi, kmin, kmax)
    if copy:
        htilde = FrequencySeries(out, delta_f=htilde.delta_f,
                                 epoch=htilde.epoch, copy=False)
    return htilde


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
