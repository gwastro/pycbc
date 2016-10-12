# Copyright (C) 2016  Collin Capano
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
"""
This modules provides utilities for manipulating parameter boundaries. Namely,
classes are offered that will map values to a specified domain using either
cyclic boundaries or reflected boundaries.
"""

import numpy

class CyclicBounds(object):
    """Applies cyclic bounds to the given values."""

    def __init__(self, min_bound, max_bound):
        self._min = min_bound
        self._max = max_bound
        self._mod = max_bound - min_bound

    @property
    def min_bound(self):
        return self._min

    @property
    def max_bound(self):
        return self._max

    def __call__(self, value):
        return (value - self._min) %(self._mod) + self._min


def _reflect_left(value, boundary):
    """Given a value and a boundary, reflects the value to the left if
    the value is > the boundary."""
    if value > boundary:
        value = 2*boundary - value
    return value


def _reflect_right(value, boundary):
    """Given a value and a boundary, reflects the value to the right if the
    value is < the boundary."""
    if value < boundary:
        value = 2*boundary - value
    return value


def _reflect_well(value, min_bound, max_bound):
    """Given two boundaries, reflects the value until it falls within both
    boundaries. This is done iteratively, reflecting left off of the
    `max_bound`, then right off of the `min_bound`, etc."""
    while (value < min_bound or value > max_bound):
        value = _reflect_left(value, max_bound)
        value = _reflect_right(value, min_bound)
    return value

reflect_left = numpy.vectorize(_reflect_left)
reflect_right = numpy.vectorize(_reflect_right)
reflect_well = numpy.vectorize(_reflect_well)


class ReflectedBounds(object):
    """Applies reflected bounds to the given values."""

    def __init__(self, min_bound=None, max_bound=None):
        # check that at least one bound was provided
        if min_bound is None and max_bound is None:
            raise ValueError("must provide at least one bound")
        if min_bound is not None and max_bound is not None and \
                min_bound >= max_bound:
            raise ValueError("min_bound must be < max_bound")
        self._min = min_bound
        self._max = max_bound

    @property
    def min_bound(self):
        return self._min

    @property
    def max_bound(self):
        return self._max

    def __call__(self, value):
        left_bounded = self._min is not None
        right_bounded = self._max is not None
        if left_bounded and right_bounded:
            return reflect_well(value, self._min, self._max)
        elif left_bounded:
            return reflect_right(value, self._min)
        else:
            return reflect_left(value, self._max)
