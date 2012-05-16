# Copyright (C) 2012  Tito Dal Canton
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
"""
This module provides a class representing a frequency series.
"""

from pycbc.array import Array
import swiglal as _swiglal

class FrequencySeries(Array):
    def __init__(self, initial_array, delta_f, epoch=None, dtype=None, copy=True):
        """initial_array: array containing sampled data.
        delta_f: frequency between consecutive samples in Hertz.
        epoch: start time of the associated time domain data, in seconds.
               Must be a swiglal.LIGOTimeGPS object.
        """
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if not delta_f > 0:
            raise ValueError('delta_f must be a positive number')
        if epoch is not None and not isinstance(epoch, _swiglal.LIGOTimeGPS):
            raise TypeError('epoch must be either None or a swiglal.LIGOTimeGPS')
        Array.__init__(self, initial_array, dtype=dtype, copy=copy)
        self._delta_f = delta_f
        self._epoch = epoch
    
    def _return(self, ary):
        return FrequencySeries(ary, self._delta_f, epoch=self._epoch, copy=False)

    def _typecheck(self, other):
        if not isinstance(other, Array):
            return NotImplemented
        if isinstance(other, FrequencySeries):
            if other._delta_f != self._delta_f:
                raise ValueError('different delta_f')
            # consistency of _epoch is not required because we may want
            # to combine frequency series estimated at different times
            # (e.g. PSD estimation)
    
    def get_delta_f(self):
        "Return frequency between consecutive samples in Hertz."
        return self._delta_f
    delta_f = property(get_delta_f)

    def get_sample_frequencies(self):
        "Return an Array containing the sample frequencies."
        return Array(range(len(self))) * self._delta_f
    sample_frequencies = property(get_sample_frequencies)

