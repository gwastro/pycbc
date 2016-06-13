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
from pycbc.types import TimeSeries

def calculate_acf(data):
    """ Calculates the autocorrelation function (ACF) and returns the one-sided
    ACF.

    ACF is estimated using

        \hat{R}(k) = \frac{1}{\left( n-k \right) \sigma^{2}} \sum_{t=1}^{n-k} \left( X_{t} - \mu \right) \left( X_{t+k} - \mu \right) 

    Where \hat{R}(k) is the ACF, X_{t} is the data series at time t, \mu is the
    mean of X_{t}, and \sigma^{2} is the variance of X_{t}.

    Parameters
    -----------
    data : {TimeSeries, numpy.array}
        A TimeSeries or numpy.array of data.

    Returns
    -------
    acf : numpy.array
        If data is a TimeSeries then acf will be a TimeSeries of the
        one-sided ACF. Else acf is a numpy.array.
    """

    # if given a TimeSeries instance then get numpy.array
    if isinstance(data, TimeSeries):
        y = data.numpy()
    else:
        y = data

    # subtract mean
    z =  y - y.mean()

    # autocorrelate
    acf = numpy.correlate(z, z, mode="full")

    # take only the second half of the autocorrelation function
    acf = acf[acf.size/2:]

    # normalize
    # note that ACF is function of k and we have a factor of n-k
    # at each k so the array here is a vectorized version of computing it
    acf /= ( y.var() * numpy.arange(y.size, 0, -1) )

    # return a TimeSeries if input was a TimeSeries
    # otherwise return the numpy.array
    if isinstance(data, TimeSeries):
        return TimeSeries(acf, delta_t=data.delta_t)
    else:
        return acf

def calculate_acl(data, m=5, k=2, dtype=int):
    """ Calculates the autocorrelation length (ACL).

    ACL is estimated using

        r = 1 + 2 \sum_{k=1}^{n} \hat{R}(k)

    Where r is the ACL and \hat{R}(k) is the ACF.

    The parameter k sets the maximum samples to use in calculation of ACL.

    The parameter m controls the length of the window that is summed to
    compute the ACL.

    Parameters
    -----------
    data : {TimeSeries, numpy.array}
        A TimeSeries or numpy.array of data.
    dtype : {int, float}
        The datatype of the output.

    Returns
    -------
    acl : {int, float}
        The ACL. If ACL can not be estimated then returns numpy.inf. If data
        was a TimeSeries then the ACL will be the number of samples times the
        delta_t attribute.
    """

    # calculate ACF
    acf = calculate_acf(data)

    # get maximum index in ACF for calculation
    # will terminate at this index
    n = len(data)
    imax = n / k

    # loop over ACF indexes
    acl = 1.0
    for i in range(n):

        # see if we have reached terminating condition
        s = float(i+1) / m
        if not acl >= s:
            break

        # see if ACL is indeterminate
        if i > imax:
            return numpy.inf

        # add term for ACL
        acl += 2 * acf[i]

    # if input was a TimeSeries then multiply ACL by the delta_t attribute
    if isinstance(data, TimeSeries):
        acl *= data.delta_t

    # typecast and return
    if dtype == int:
        return int(numpy.ceil(acl))
    elif dtype == float:
        return acl
    else:
        raise ValueError("Invalid value for dtype")
