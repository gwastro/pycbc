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
"""This module contains the CuPy-specific code for
   convenience utilities for manipulating waveforms
"""
from pycbc.types import FrequencySeries
import cupy as xp



def apply_fseries_time_shift(htilde, dt, kmin=0, copy=True):
    """Shifts a frequency domain waveform in time. The waveform is assumed to
    be sampled at equal frequency intervals.
    """
    out = xp.array(htilde.data, copy=copy)
    phi = -2j * xp.pi * dt * htilde.delta_f
    kmax = len(htilde)
    fstimeshift(out, phi, kmin, kmax)
    if copy:
        htilde = FrequencySeries(out, delta_f=htilde.delta_f,
                                 epoch=htilde.epoch, copy=False)
    return htilde


def fstimeshift(freqseries, phi, kmin, kmax):
    # Please note that the CPU version has a number of optimizations with respect
    # to this.
    # FIXME: Convert to ElementwiseKernel and use kmin and max in that.
    idx = xp.arange(len(freqseries))
    phase_shift = xp.exp(phi * idx)
    freqseries[:] = freqseries[:] * phase_shift 

