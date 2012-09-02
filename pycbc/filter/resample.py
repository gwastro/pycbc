# Copyright (C) 2012  Alex Nitz
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
import lal
import numpy
from pycbc.types import TimeSeries

_resample_func = {numpy.dtype('float32'): lal.ResampleREAL4TimeSeries,
                 numpy.dtype('float64'): lal.ResampleREAL8TimeSeries}

def resample_to_delta_t(timeseries, delta_t):
    """Return a new time series that is resampled to the given delta_t. Only powers
       of two are currently supported.
    """
    if not isinstance(timeseries,TimeSeries):
        raise TypeError("Can only resample time series")

    if timeseries.kind is not 'real':
        raise TypeError("Time series must be real")

    if timeseries.delta_t == delta_t:
        return timeseries * 1

    lal_data = timeseries.lal()
    _resample_func[timeseries.dtype](lal_data, delta_t)

    return TimeSeries(lal_data.data.data, delta_t = lal_data.deltaT,
                      dtype=timeseries.dtype, epoch=timeseries._epoch)

    

__all__ = ['resample_to_delta_t']

