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
and length of time series.
"""

import numpy

def calculate_acf(data):
    """ Calculates the autocorrelation function (ACF) and returns the one-sided
    ACF.
    """

    # subtract mean
    z =  data - data.mean()

    # autocorrelate
    acf = numpy.correlate(z, z, mode="full")

    # take only the second half of the autocorrelation time series
    acf = acf[acf.size/2:]

    # normalize
    acf /= data.var() * numpy.arange(data.size, 0, -1)

    return acf

def calculate_acl(data, m=5, k=2):
    """ Calculates the autocorrelation length (ACL).
    """

    # calculate ACF
    # get number of points in ACF
    acf = calculate_acf(data)
    n = len(acf)

    print acf

    # get maximum index in ACF for calculation
    # will terminate at this index
    imax = n / k

    # loop over ACF indexes
    acl = 1.0
    for i in range(n):

        # see if we have reached terminating condition
        s = float(i+1) / m
        if not acl >= s:
            print acl, s
            break

        # see if ACL is indeterminate
        if s > imax:
            return numpy.inf

        # add term for ACL
        acl += 2 * acf[i]

    return numpy.ceil(acl)
