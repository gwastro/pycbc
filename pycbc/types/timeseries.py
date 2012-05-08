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
This module provides a class representing a time series.
"""

from pycbc.types.array import Array,_convert
import swiglal as _swiglal
import numpy as _numpy

class TimeSeries(Array):
    def __init__(self, initial_array, delta_t, epoch=None, dtype=None, copy=True):
        """initial_array: array containing sampled data.        
        delta_t: time between consecutive samples in seconds.
        epoch: time series start time in seconds. Must be a swiglal.LIGOTimeGPS object.
        """
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if not delta_t > 0:
            raise ValueError('delta_t must be a positive number')
        if epoch is not None and not isinstance(epoch, _swiglal.LIGOTimeGPS):
            raise TypeError('epoch must be either None or a swiglal.LIGOTimeGPS')
        Array.__init__(self, initial_array, dtype=dtype, copy=copy)
        self._delta_t = delta_t
        self._epoch = epoch
    
    def _return(self, ary):
        return TimeSeries(ary, self._delta_t, epoch=self._epoch, copy=False)
    
    def _typecheck(self, other):
        if not isinstance(other, Array):
            return NotImplemented
        if isinstance(other, TimeSeries):
            if other._delta_t != self._delta_t:
                raise ValueError('different delta_t')
            if self._epoch != other._epoch:
                raise ValueError('different epoch')

    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.start > index.stop:
                raise ValueError('start index must be smaller than stop index')
            if index.step is None:
                new_delta_t = self._delta_t
            elif index.step < 0:
                raise ValueError('negative step is not supported')
            else:
                new_delta_t = index.step * self._delta_t
            if self._epoch is None:
                new_epoch = None
            else:
                new_epoch = self._epoch + index.start * self._delta_t
            return TimeSeries(Array.__getitem__(self, index), new_delta_t, new_epoch)
        else:
            return Array.__getitem__(self, index)

    def get_delta_t(self):
        "Return time between consecutive samples in seconds."
        return self._delta_t
    delta_t = property(get_delta_t)

    def get_duration(self):
        "Return duration of time series in seconds."
        return len(self) * self._delta_t
    duration = property(get_duration)

    def get_start_time(self):
        "Return time series start time as a LIGOTimeGPS."
        if self._epoch:
            return self._epoch
        else:
            return None
    start_time = property(get_start_time)

    def get_end_time(self):
        "Return time series end time as a LIGOTimeGPS."
        if self._epoch:
            return self._epoch + self.get_duration()
        else:
            return None
    end_time = property(get_end_time)

    def get_sample_times(self):
        "Return an Array containing the sample times."
        if self._epoch is None:
            return Array(range(len(self))) * self._delta_t
        else:
            return Array(range(len(self))) * self._delta_t + float(self._epoch)
    sample_times = property(get_sample_times)

    def resample(self, new_delta_t):
        "Return a resampled TimeSeries with the specified delta_t."
        return NotImplemented
            
    @property
    @_convert
    def  lal(self):
        """ Returns a LAL Object that contains this data """
        lal_data = None
        epoch = _swiglal.LIGOTimeGPS()
        if self._epoch is not None:
            epoch = self._epoch            
        
        if type(self._data) is not _numpy.ndarray:
            raise TypeError("Cannot return lal type from the GPU")
        elif self._data.dtype == _numpy.float32:
            lal_data = _swiglal.XLALCreateREAL4TimeSeries("",epoch,
                                    0,self._delta_t,_swiglal.LALUnit(),len(self))
        elif self._data.dtype == _numpy.float64:
            lal_data = _swiglal.XLALCreateREAL8TimeSeries("",epoch,
                                    0,self._delta_t,_swiglal.LALUnit(),len(self))
        elif self._data.dtype == _numpy.complex64:
            lal_data = _swiglal.XLALCreateCOMPLEX8TimeSeries("",epoch,
                                    0,self._delta_t,_swiglal.LALUnit(),len(self))
        elif self._data.dtype == _numpy.complex128:
            lal_data = _swiglal.XLALCreateCOMPLEX16TimeSeries("",epoch,
                                    0,self._delta_t,_swiglal.LALUnit(),len(self))

        lal_data.data.data = self._data
        self.data = lal_data.data.data

        return lal_data
