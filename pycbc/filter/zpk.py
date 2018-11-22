# Copyright (C) 2014  Christopher M. Biwer
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

import numpy as np

from scipy.signal import zpk2sos, sosfilt
from pycbc.types import TimeSeries

def filter_zpk(timeseries, z, p, k):
    """Return a new timeseries that was filtered with a zero-pole-gain filter.
    The transfer function in the s-domain looks like:
    .. math::
    \\frac{H(s) = (s - s_1) * (s - s_3) * ... * (s - s_n)}{(s - s_2) * (s - s_4) * ... * (s - s_m)}, m >= n

    The zeroes, and poles entered in Hz are converted to angular frequency,
    along the imaginary axis in the s-domain s=i*omega.  Then the zeroes, and
    poles are bilinearly transformed via:
    .. math::
    z(s) = \\frac{(1 + s*T/2)}{(1 - s*T/2)}

    Where z is the z-domain value, s is the s-domain value, and T is the
    sampling period.  After the poles and zeroes have been bilinearly
    transformed, then the second-order sections are found and filter the data
    using scipy.

    Parameters
    ----------
    timeseries: TimeSeries
        The TimeSeries instance to be filtered.
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

    Examples
    --------
    To apply a 5 zeroes at 100Hz, 5 poles at 1Hz, and a gain of 1e-10 filter
    to a TimeSeries instance, do:
    >>> filtered_data = zpk_filter(timeseries, [100]*5, [1]*5, 1e-10)
    """

    # sanity check type
    if not isinstance(timeseries, TimeSeries):
        raise TypeError("Can only filter TimeSeries instances.")

    # sanity check casual filter
    degree = len(p) - len(z)
    if degree < 0:
        raise TypeError("May not have more zeroes than poles. \
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
    fs = 2.0 * timeseries.sample_rate

    # zeroes in the z-domain
    z_zd = (1 + z/fs) / (1 - z/fs)

    # any zeros that were at infinity are moved to the Nyquist frequency
    z_zd = z_zd[np.isfinite(z_zd)]
    z_zd = np.append(z_zd, -np.ones(degree))

    # poles in the z-domain
    p_zd = (1 + p/fs) / (1 - p/fs)

    # gain change in z-domain
    k_zd = k * np.prod(fs - z) / np.prod(fs - p)

    # get second-order sections
    sos = zpk2sos(z_zd, p_zd, k_zd)

    # filter
    filtered_data = sosfilt(sos, timeseries.numpy())

    return TimeSeries(filtered_data, delta_t = timeseries.delta_t,
                      dtype=timeseries.dtype,
                      epoch=timeseries._epoch)
