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

"""
Provides a class representing a time series.
"""

from pycbc.types.array import Array,_convert
import lal as _lal
import numpy as _numpy

class TimeSeries(Array):
    """Models a time series consisting of uniformly sampled scalar values.

    Parameters
    ----------
    initial_array : array-like
        Array containing sampled data.
    delta_t : float
        Time between consecutive samples in seconds.
    epoch : {None, lal.LIGOTimeGPS}, optional
        Time of the first sample in seconds.
    dtype : {None, data-type}, optional
        Sample data type.
    copy : boolean, optional
        If True, samples are copied to a new array.

    Attributes
    ----------
    delta_t
    duration
    start_time
    end_time
    sample_times
    """

    def __init__(self, initial_array, delta_t, epoch=None, dtype=None, copy=True):
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if not delta_t > 0:
            raise ValueError('delta_t must be a positive number')
        if epoch is not None and not isinstance(epoch, _lal.LIGOTimeGPS):
            raise TypeError('epoch must be either None or a lal.LIGOTimeGPS')
        if epoch is None:
            epoch = _lal.LIGOTimeGPS(0,0)
        else:
            epoch = _lal.LIGOTimeGPS(epoch)
        Array.__init__(self, initial_array, dtype=dtype, copy=copy)
        self._delta_t = delta_t
        self._epoch = epoch

    def _return(self, ary):
        return TimeSeries(ary, self._delta_t, epoch=self._epoch, copy=False)

    def _typecheck(self, other):
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
            return TimeSeries(Array.__getitem__(self, index), new_delta_t, new_epoch, copy=False)
        else:
            return Array.__getitem__(self, index)

    def get_delta_t(self):
        """Return time between consecutive samples in seconds.
        """
        return self._delta_t
    delta_t = property(get_delta_t)

    def get_duration(self):
        """Return duration of time series in seconds.
        """
        return len(self) * self._delta_t
    duration = property(get_duration)

    def get_start_time(self):
        """Return time series start time as a LIGOTimeGPS.
        """
        return self._epoch
    start_time = property(get_start_time)

    def get_end_time(self):
        """Return time series end time as a LIGOTimeGPS.
        """
        return self._epoch + self.get_duration()
    end_time = property(get_end_time)

    def get_sample_times(self):
        """Return an Array containing the sample times.
        """
        if self._epoch is None:
            return Array(range(len(self))) * self._delta_t
        else:
            return Array(range(len(self))) * self._delta_t + float(self._epoch)
    sample_times = property(get_sample_times)

    def __eq__(self,other):
        """
        This is the Python special method invoked whenever the '=='
        comparison is used.  It will return true if the data of two
        time series are identical, and all of the numeric meta-data
        are identical, irrespective of whether or not the two
        instances live in the same memory (for that comparison, the
        Python statement 'a is b' should be used instead).

        Thus, this method returns 'True' if the types of both 'self'
        and 'other' are identical, as well as their lengths, dtypes,
        epochs, delta_ts and the data in the arrays, element by element.
        It will always do the comparison on the CPU, but will *not* move
        either object to the CPU if it is not already there, nor change
        the scheme of either object. It is possible to compare a CPU
        object to a GPU object, and the comparison should be true if the
        data and meta-data of the two objects are the same.

        Note in particular that this function returns a single boolean,
        and not an array of booleans as Numpy does.  If the numpy
        behavior is instead desired it can be obtained using the numpy()
        method of the PyCBC type to get a numpy instance from each
        object, and invoking '==' on those two instances.

        Parameters
        ----------
        other: another Python object, that should be tested for equality
            with 'self'.

        Returns
        -------
        boolean: 'True' if the types, dtypes, lengths, epochs, delta_ts
            and data of the two objects are each identical.
        """
        if super(TimeSeries,self).__eq__(other):
            return (self._epoch == other._epoch and self._delta_t == other._delta_t)
        else:
            return False

    def almost_equal_elem(self,other,tol,relative=True):
        """
        Compare whether two time series are almost equal, element
        by element.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
             abs(self[i]-other[i]) <= tol*abs(self[i])
        for all elements of the series.

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
             abs(self[i]-other[i]) <= tol
        for all elements of the series.

        Other meta-data (type, dtype, length, epoch, and delta_t) must
        be exactly equal.  If either object's memory lives on the GPU it
        will be copied to the CPU for the comparison, which may be slow.
        But the original object itself will not have its memory relocated
        nor scheme changed.

        Parameters
        ----------
        other: another Python object, that should be tested for
            almost-equality with 'self', element-by-element.
        tol: a non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative: A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).

        Returns
        -------
        boolean: 'True' if the data agree within the tolerance, as
            interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, epochs, and delta_ts are exactly the same.
        """
        if super(TimeSeries,self).almost_equal_elem(other):
            return (self._epoch == other._epoch and self._delta_t == other._delta_t)
        else:
            return False

    def almost_equal_norm(self,other,tol,relative=True):
        """
        Compare whether two time series are almost equal, normwise.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
             abs(norm(self-other)) <= tol*abs(norm(self)).

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
             abs(norm(self-other)) <= tol

        Other meta-data (type, dtype, length, epoch, and delta_t) must
        be exactly equal.  If either object's memory lives on the GPU it
        will be copied to the CPU for the comparison, which may be slow.
        But the original object itself will not have its memory relocated
        nor scheme changed.

        Parameters
        ----------
        other: another Python object, that should be tested for
            almost-equality with 'self', based on their norms.
        tol: a non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative: A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).

        Returns
        -------
        boolean: 'True' if the data agree within the tolerance, as
            interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, epochs, and delta_ts are exactly the same.
        """
        if super(TimeSeries,self).almost_equal_norm(other):
            return (self._epoch == other._epoch and self._delta_t == other._delta_t)
        else:
            return False

    @_convert
    def lal(self):
        """Produces a LAL time series object equivalent to self.

        Returns
        -------
        lal_data : {lal.*TimeSeries}
            LAL time series object containing the same data as self.
            The actual type depends on the sample's dtype.

        Raises
        ------
        TypeError
            If time series is stored in GPU memory.
        """
        lal_data = None
        if type(self._data) is not _numpy.ndarray:
            raise TypeError("Cannot return lal type from the GPU")
        elif self._data.dtype == _numpy.float32:
            lal_data = _lal.CreateREAL4TimeSeries("",self._epoch,0,self.delta_t,_lal.lalSecondUnit,len(self))
        elif self._data.dtype == _numpy.float64:
            lal_data = _lal.CreateREAL8TimeSeries("",self._epoch,0,self.delta_t,_lal.lalSecondUnit,len(self))
        elif self._data.dtype == _numpy.complex64:
            lal_data = _lal.CreateCOMPLEX8TimeSeries("",self._epoch,0,self.delta_t,_lal.lalSecondUnit,len(self))
        elif self._data.dtype == _numpy.complex128:
            lal_data = _lal.CreateCOMPLEX16TimeSeries("",self._epoch,0,self.delta_t,_lal.lalSecondUnit,len(self))

        lal_data.data.data[:] = self._data

        return lal_data

