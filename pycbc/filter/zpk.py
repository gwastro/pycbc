# Copyright (C) 2014  Christopher M. Biwer
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

import numpy
import scipy.signal

from pycbc.types import TimeSeries

def filter_zpk(timeseries, z, p, k):
    """Return a new timeseries that is filter with zpk (zeros, poles, gain)
       parameters.

    Parameters
    ----------
    timeseries: TimeSeries
        The time series to be filtered.
    z: array
        Array of zeros to include in zpk filter design, eg. 3 zeroes at 1Hz
        would be array([1., 1., 1.])
    p: array
        Array of poles to include in zpk filter design, eg. 3 poles at 100Hz
        would be array([-100., -100., -100.])
    k: float
        Gain to include in zpk filter design. This gain is a contast
        multiplied to the transfer function.

    Returns
    -------
    Time Series: TimeSeries
        A  new TimeSeries that has been filtered. 
    """

    if not isinstance(timeseries,TimeSeries):
        raise TypeError("Can only filter time series.")

    if len(z) > len(p):
        raise TypeError("Improper tansfer function, more numerator terms than denominator.")

    # get digital filter coefficients
    _, digital_coefficients = coefficients_zpk(z, p, k, timeseries.sample_rate)

    # apply the filter
    series = scipy.signal.lfilter(digital_coefficients[0], digital_coefficients[1], \
                                      timeseries.numpy())

    return TimeSeries(series, delta_t = timeseries.delta_t,
                      dtype=timeseries.dtype,
                      start_time=timeseries.start_time)

def filter_zpk_factored(timeseries, z, p, k, increment=3):
    """Return a new timeseries that is filter with zpk (zeros, poles, gain)
       parameters. This function factors the zpk filter into smaller filters,
       and applies each factor seperately. Only use this function instead of
       filter_zpk with higher order filters; this is because in the frequency
       response of the digital filter there may be large error at low frequency.

    Parameters
    ----------
    timeseries: TimeSeries
        The time series to be filtered.
    z: array
        Array of zeros to include in zpk filter design, eg. 3 zeroes at 1Hz
        would be array([1., 1., 1.])
    p: array
        Array of poles to include in zpk filter design, eg. 3 poles at 100Hz
        would be array([-100., -100., -100.])
    k: float
        Gain to include in zpk filter design. This gain is a contast
        multiplied to the transfer function.

    Returns
    -------
    Time Series: TimeSeries
        A  new TimeSeries that has been filtered. 
    """

    if not isinstance(timeseries,TimeSeries):
        raise TypeError("Can only filter time series.")

    if len(z) > len(p):
        raise TypeError("Improper tansfer function, more numerator terms than denominator.")

    # get max of fraction denominator and numerator
    length = max(len(z), len(p))

    # split up input into smaller filters
    # and filter the data with each small filter
    i = 0
    j = increment
    series = timeseries.numpy()
    for r in range(length/increment+1):
        # get digital filter coefficients
        _, digital_coefficients = coefficients_zpk(z[i:j], p[i:j], k, \
                                                 timeseries.sample_rate)

        # apply the filter
        series = scipy.signal.lfilter(digital_coefficients[0], \
                                        digital_coefficients[1], series)

        # increment
        i += increment
        j += increment

        # for subsequent filters do not apply gain,
        # since it was applied in first filter
        k = 1

    return TimeSeries(series, delta_t = timeseries.delta_t,
                      dtype=timeseries.dtype,
                      start_time=timeseries.start_time)
 
def coefficients_zpk(z, p, k, sample_rate):
    """Return anolog and digital coefficients for a zero-pole-gain filter.

    Parameters
    ----------
    z: list
        List of zeros to include in zpk filter design, eg. 3 zeroes at 1Hz
        would be array([1., 1., 1.])
    p: list
        list of poles to include in zpk filter design, eg. 3 poles at 100Hz
        would be array([-100., -100., -100.])
    k: float
        Gain to include in zpk filter design. This gain is a contast
        multiplied to the transfer function.

    Returns
    -------
    analog_coefficients: tuple
        A  tuple of (numerator_coefficients, denominator_coefficients).
    digital_coefficients: tuple
        A  tuple of (numerator_coefficients, denominator_coefficients). 
    """

    # create transfer function coefficients using analog filter
    analog_coefficients = scipy.signal.zpk2tf(z, p, k)

    # convert to digital filter
    digital_coefficients = scipy.signal.bilinear(analog_coefficients[0], \
                                       analog_coefficients[1], sample_rate)

    return analog_coefficients, digital_coefficients
