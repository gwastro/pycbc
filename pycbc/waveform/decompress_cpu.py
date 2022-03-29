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
import numpy
from ..types import real_same_precision_as
from ..types import complex_same_precision_as
from .decompress_cpu_cython import decomp_ccode_double, decomp_ccode_float

def inline_linear_interp(amp, phase, sample_frequencies, output,
                         df, f_lower, imin, start_index):

    rprec = real_same_precision_as(output)
    cprec = complex_same_precision_as(output)
    sample_frequencies = numpy.array(sample_frequencies, copy=False,
                                     dtype=rprec)
    amp = numpy.array(amp, copy=False, dtype=rprec)
    phase = numpy.array(phase, copy=False, dtype=rprec)
    sflen = len(sample_frequencies)
    h = numpy.array(output.data, copy=False, dtype=cprec)
    hlen = len(output)
    delta_f = float(df)
    if output.precision == 'single':
        decomp_ccode_float(h, delta_f, hlen, start_index, sample_frequencies,
                           amp, phase, sflen, imin)
    else:
        decomp_ccode_double(h, delta_f, hlen, start_index, sample_frequencies,
                            amp, phase, sflen, imin)

    return output
