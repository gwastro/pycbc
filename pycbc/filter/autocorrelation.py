# Copyright (C) 2016  Christopher M. Biwer
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
This modules provides functions for calculating the autocorrelation function
and length of a data series.
"""

import numpy
from math import isnan
from pycbc.filter.matchedfilter import correlate
from pycbc.types import FrequencySeries, TimeSeries, zeros

def calculate_acf(data, delta_t=1.0, unbiased=False):
    r"""Calculates the one-sided autocorrelation function.

    Calculates the autocorrelation function (ACF) and returns the one-sided
    ACF. The ACF is defined as the autocovariance divided by the variance. The
    ACF can be estimated using

    .. math::

        \hat{R}(k) = \frac{1}{n \sigma^{2}} \sum_{t=1}^{n-k} \left( X_{t} - \mu \right) \left( X_{t+k} - \mu \right) 

    Where :math:`\hat{R}(k)` is the ACF, :math:`X_{t}` is the data series at
    time t, :math:`\mu` is the mean of :math:`X_{t}`, and :math:`\sigma^{2}` is
    the variance of :math:`X_{t}`.

    Parameters
    -----------
    data : TimeSeries or numpy.array
        A TimeSeries or numpy.array of data.
    delta_t : float
        The time step of the data series if it is not a TimeSeries instance.
    unbiased : bool
        If True the normalization of the autocovariance function is n-k
        instead of n. This is called the unbiased estimation of the
        autocovariance. Note that this does not mean the ACF is unbiased.

    Returns
    -------
    acf : numpy.array
        If data is a TimeSeries then acf will be a TimeSeries of the
        one-sided ACF. Else acf is a numpy.array.
    """

    # if given a TimeSeries instance then get numpy.array
    if isinstance(data, TimeSeries):
        y = data.numpy()
        delta_t = data.delta_t
    else:
        y = data

    # FFT data minus the mean
    fdata = TimeSeries(y-y.mean(), delta_t=delta_t).to_frequencyseries()

    # correlate
    # do not need to give the congjugate since correlate function does it
    cdata = FrequencySeries(zeros(len(fdata), dtype=numpy.complex64),
                           delta_f=fdata.delta_f, copy=False)
    correlate(fdata, fdata, cdata)

    # IFFT correlated data to get unnormalized autocovariance time series
    acf = cdata.to_timeseries()

    # normalize the autocovariance
    # note that dividing by acf[0] is the same as ( y.var() * len(acf) )
    if unbiased:
        acf /= ( y.var() * numpy.arange(len(acf), 0, -1) )
    else:
        acf /= acf[0]

    # return input datatype
    if isinstance(data, TimeSeries):
        return TimeSeries(acf, delta_t=delta_t)
    else:
        return acf

def calculate_acl(data, m=5, k=2, dtype=int):
    """ Calculates the autocorrelation length (ACL).

    ACL is estimated using

        r = 1 + 2 \sum_{i=1}^{m*s} \hat{R}(i) < s

    Where r is the ACL and \hat{R}(i) is the ACF that has been normalized so
    that \hat{R}(0) is 1.0. And s is equal to i/m.

    The parameter k sets the maximum samples to use in calculation of ACL. The
    maximum number of samples will be the length of the ACL divided by k.

    The parameter m controls the length of the window that is summed to
    compute the ACL.

    Parameters
    -----------
    data : {TimeSeries, numpy.array}
        A TimeSeries or numpy.array of data.
    dtype : {int, float}
        The datatype of the output. If the dtype was set to int, then the
        ceiling is returned.

    Returns
    -------
    acl : {int, float}
        The length s which is longer than the ACL. If ACL can not be estimated
        then returns numpy.inf.
    """

    # sanity check output data type
    if dtype not in [int, float]:
        raise ValueError("The dtype must be either int or float.")

    # calculate ACF that is normalized by the zero-lag value
    acf = calculate_acf(data)

    # multiply all values beyond the zero-lag value by 2.0
    acf[1:] *= 2.0

    # sanity check ACF
    if isnan(acf[0]):
        return numpy.inf
    assert acf[0] == 1.0

    # the maximum index to calculate ACF
    imax = int(len(acf)/k)

    # calculate cumlative ACL until s is less than the cumulative ACL
    cum_acl = 0.0
    for i,val in enumerate(acf[:imax]):
        s = float(i)/m
        if cum_acl+val < s:
            if dtype == int:
                return numpy.ceil(s)
            elif dtype == float:
                return s
        cum_acl += val

    return numpy.inf


