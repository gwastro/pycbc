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

def resample_to_delta_t(timeseries, delta_t, method='butterworth'):
    """Resmple the time_series to delta_t

    Resamples the TimeSeries instance time_series to the given time step, 
    delta_t. Only powers of two and real valued time series are supported 
    at this time. 

    Parameters
    ----------
    time_series: TimeSeries
        The time series to be resampled
    delta_t: float
        The desired time step 

    Returns
    -------
    Time Series: TimeSeries
        A TimeSeries that has been resampled to delta_t.

    Raises
    ------
    TypeError: 
        time_series is not an instance of TimeSeries.
    TypeError: 
        time_series is not real valued

    Examples
    --------

    >>> h_plus_sampled = resample_to_delta_t(h_plus, 1.0/2048)
    """

    if not isinstance(timeseries,TimeSeries):
        raise TypeError("Can only resample time series")

    if timeseries.kind is not 'real':
        raise TypeError("Time series must be real")

    if timeseries.delta_t == delta_t:
        return timeseries * 1

    if method == 'butterworth':
        lal_data = timeseries.lal()
        _resample_func[timeseries.dtype](lal_data, delta_t)
        data = lal_data.data.data 
    elif method == 'ldas':
        pass
    elif method == 'broken_butterworth':    
        # Do a low pass filter to remove power at frequencies above the
        # new nyquist frequency
        
        #Yep this doesn't make sense, but is what the code does here
        # http://www.lsc-group.phys.uwm.edu/cgit/lalsuite/tree/lal/packages/tools/src/ResampleTimeSeries.c#n408
        lowpass_params = lal.PassBandParamStruc()
        lowpass_params.nMax = 20
        lowpass_params.f1 = 0.5 / delta_t
        lowpass_params.a1 = 0.1
        lowpass_params.f2 = lal.LAL_REAL4_MAX
        lowpass_params.a1 = 0
        
        lal_data = timeseries.lal()
        lal.DButterworthREAL4TimeSeries(lal_data, lowpass_params)       
        
        # Decimate the time series
        factor = int(delta_t / timeseries.delta_t)
        data = lal_data.data.data[::factor] * 1
        
    return TimeSeries(data, delta_t = delta_t,
                      dtype=timeseries.dtype, 
                      epoch=timeseries._epoch)
       

_highpass_func = {numpy.dtype('float32'): lal.HighPassREAL4TimeSeries,
                 numpy.dtype('float64'): lal.HighPassREAL8TimeSeries}

def highpass(timeseries, frequency, filter_order=8, attenuation=0.1):
    """Return a new timeseries that is highpassed.

    Return a new time series that is highpassed above the `frequency`. 

    Parameters
    ----------
    Time Series: TimeSeries
        The time series to be high-passed.
    frequency: float
        The frequency below which is suppressed. 
    filter_order: {8, int}, optional
        The order of the filter to use when high-passing the time series.
    attenuation: {0.1, float}, optional
        The attenuation of the filter. 

    Returns
    -------
    Time Series: TimeSeries
        A  new TimeSeries that has been high-passed. 

    Raises
    ------
    TypeError: 
        time_series is not an instance of TimeSeries.
    TypeError: 
        time_series is not real valued

    """

    if not isinstance(timeseries, TimeSeries):
        raise TypeError("Can only resample time series")

    if timeseries.kind is not 'real':
        raise TypeError("Time series must be real")

    lal_data = timeseries.lal()
    _highpass_func[timeseries.dtype](lal_data, frequency, 1-attenuation, filter_order)

    return TimeSeries(lal_data.data.data, delta_t = lal_data.deltaT,
                      dtype=timeseries.dtype, epoch=timeseries._epoch)

    

__all__ = ['resample_to_delta_t', 'highpass']

