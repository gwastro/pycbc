# Copyright (C) 2024 Y Ddraig Goch
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

import cupy as cp
from .matchedfilter import _BaseCorrelator

# Here X,Y,Z are "type placeholder"s, so this covers 32 and 64 bit inputs.
# It would work for real -> real as well, except maybe conj would fail.
# I've also made this work for mixed types (32 bit -> 64 bit), but this means
# we need to always supply the output, which we do.
correlate_kernel = cp.ElementwiseKernel(
    "X x, Y y",
    "Z z",
    "z = conj(x) * y",
    "correlate_kernel"
)

def correlate(a, b, out):
    correlate_kernel(a.data, b.data, out.data)

class CUPYCorrelator(_BaseCorrelator):
    def __init__(self, x, y, z):
        self.x = x.data
        self.y = y.data
        self.z = z.data

    def correlate(self):
        correlate_kernel(self.x, self.y, self.z)

def _correlate_factory(x, y, z):
    return CUPYCorrelator


