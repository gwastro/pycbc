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
Provides a class representing a frequency series.
"""

from pycbc.types.array import Array,_convert
import lal as _lal
import numpy as _numpy

class FrequencySeries(Array):
    """Models a frequency series consisting of uniformly sampled scalar values.

    Parameters
    ----------
    initial_array : array-like
        Array containing sampled data.
    delta_f : float
        Frequency between consecutive samples in Hertz.
    epoch : {None, lal.LIGOTimeGPS}, optional
        Start time of the associated time domain data in seconds.
    dtype : {None, data-type}, optional
        Sample data type.
    copy : boolean, optional
        If True, samples are copied to a new array.

    Attributes
    ----------
    delta_f
    epoch
    sample_frequencies
    """

    def __init__(self, initial_array, delta_f, epoch=None, dtype=None, copy=True):
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if not delta_f > 0:
            raise ValueError('delta_f must be a positive number')
        if epoch is not None and not isinstance(epoch, _lal.LIGOTimeGPS):
            raise TypeError('epoch must be either None or a lal.LIGOTimeGPS')
        if epoch is None:
            epoch = _lal.LIGOTimeGPS(0,0)
        else:
            epoch = _lal.LIGOTimeGPS(epoch)
        Array.__init__(self, initial_array, dtype=dtype, copy=copy)
        self._delta_f = delta_f
        self._epoch = epoch

    def _return(self, ary):
        return FrequencySeries(ary, self._delta_f, epoch=self._epoch, copy=False)

    def _typecheck(self, other):
        if isinstance(other, FrequencySeries):
            if other._delta_f != self._delta_f:
                raise ValueError('different delta_f')
            # consistency of _epoch is not required because we may want
            # to combine frequency series estimated at different times
            # (e.g. PSD estimation)

    def get_delta_f(self):
        """Return frequency between consecutive samples in Hertz.
        """
        return self._delta_f
    delta_f = property(get_delta_f)

    def get_epoch(self):
        """Return frequency series epoch as a LIGOTimeGPS.
        """
        return self._epoch
    epoch = property(get_epoch)

    def get_sample_frequencies(self):
        """Return an Array containing the sample frequencies.
        """
        return Array(range(len(self))) * self._delta_f
    sample_frequencies = property(get_sample_frequencies)

    @_convert
    def  lal(self):
        """Produces a LAL frequency series object equivalent to self.

        Returns
        -------
        lal_data : {lal.*FrequencySeries}
            LAL frequency series object containing the same data as self.
            The actual type depends on the sample's dtype.

        Raises
        ------
        TypeError
            If frequency series is stored in GPU memory.
        """

        lal_data = None
        if type(self._data) is not _numpy.ndarray:
            raise TypeError("Cannot return lal type from the GPU")
        elif self._data.dtype == _numpy.float32:
            lal_data = _lal.CreateREAL4FrequencySeries("",self._epoch,0,self.delta_f,_lal.lalSecondUnit,len(self))
        elif self._data.dtype == _numpy.float64:
            lal_data = _lal.CreateREAL8FrequencySeries("",self._epoch,0,self.delta_f,_lal.lalSecondUnit,len(self))
        elif self._data.dtype == _numpy.complex64:
            lal_data = _lal.CreateCOMPLEX8FrequencySeries("",self._epoch,0,self.delta_f,_lal.lalSecondUnit,len(self))
        elif self._data.dtype == _numpy.complex128:
            lal_data = _lal.CreateCOMPLEX16FrequencySeries("",self._epoch,0,self.delta_f,_lal.lalSecondUnit,len(self))

        lal_data.data.data[:] = self._data

        return lal_data

