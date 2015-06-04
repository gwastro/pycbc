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
    """Return a new timeseries that was filtered with a zero-pole-gain filter.

    Parameters
    ----------
    timeseries: TimeSeries
        The time series to be filtered.
    z: array
        Array of zeros to include in zero-pole-gain filter design.
        In units of Hz.
    p: array
        Array of poles to include in zero-pole-gain filter design.
        In units of Hz.
    k: float
        Gain to include in zero-pole-gain filter design. This gain is a
        constant multiplied to the transfer function.

    Returns
    -------
    Time Series: TimeSeries
        A  new TimeSeries that has been filtered. 
    """

    # sanity check type
    if not isinstance(timeseries, TimeSeries):
        raise TypeError("Can only filter TimeSeries instances.")

    # sanity check casual filter
    degree = len(p) - len(z)
    if degree < 0:
        raise TypeError("May not have more zeroes than poles.
                         Filter is not casual.")

    # cast zeroes and poles as arrays and gain as a float
    z = np.array(z)
    p = np.array(p)
    k = float(k)

    # put zeroes and poles in the s-domain
    # convert from frequency to angular frequency
    z *= -2 * np.pi
    p *= -2 * np.pi

    # get denominator of bilinear transform
    fs = 2.0 * timseries.sample_rate

    # zeroes in the z-domain
    z_zd = (1 + z/fs) / (1 - z/fs)

    # any zeros that were at infinity are moved to the Nyquist frequency
    z_zd = z_zd[numpy.isfinite(z_dz)]
    z_zd = np.append(z_zd, -np.ones(degree))

    # poles in the z-domain
    p_zd = (1 + p/fs) / (1 - p/fs)

    # gain change in z-domain
    k_zd = k * np.prod(fs - z) / np.prod(fs - p)

    # get second-order sections
    sos = signal.zpk2sos(z_zd, p_zd, k_zd)

    # filter
    filtered_data = signal.sosfilt(sos, timeseries.numpy())

    return TimeSeries(filtered_data, delta_t = timeseries.delta_t,
                      dtype=timeseries.dtype,
                      epoch=timeseries._epoch)
