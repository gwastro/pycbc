# Copyright (C) 2016  Alex Nitz, Collin Capano
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
""" Utilities for handling frequency compressed an unequally spaced frequency
domain waveforms.
"""
from __future__ import absolute_import
import numpy
from .decompress_cpu_cython import decomp_ccode_double, decomp_ccode_float

def inline_linear_interp(amp, phase, sample_frequencies, output,
                         df, f_lower, imin, start_index):
    # FIXME: This function needs to enforce that *all* variables are the same
    #        precision and fail if they're not. Otherwise we'll get nasty C
    #        errors (or worse, unpredictable behaviour if it tries to read
    #        random memory)

    sample_frequencies = numpy.array(sample_frequencies)
    amp = numpy.array(amp)
    phase = numpy.array(phase)
    sflen = len(sample_frequencies)
    h = numpy.array(output.data, copy=False)
    hlen = len(output)
    delta_f = float(df)
    if output.precision == 'single':
        decomp_ccode_float(h, hlen, sflen, delta_f, sample_frequencies, amp,
                           phase, start_index, imin)
    else:
        decomp_ccode_double(h, hlen, sflen, delta_f, sample_frequencies, amp,
                            phase, start_index, imin)

    return output
