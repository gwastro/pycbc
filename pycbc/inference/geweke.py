# Copyright (C) 2017 Christopher M. Biwer
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
""" Functions for computing the Geweke convergence statistic.
"""

import numpy


def geweke(x, seg_length, seg_stride, end_idx, ref_start,
           ref_end=None, seg_start=0):
    """ Calculates Geweke conervergence statistic for a chain of data.
    This function will advance along the chain and calculate the
    statistic for each step.

    Parameters
    ----------
    x : numpy.array
        A one-dimensional array of data.
    seg_length : int
        Number of samples to use for each Geweke calculation.
    seg_stride : int
        Number of samples to advance before next Geweke calculation.
    end_idx : int
        Index of last start.
    ref_start : int
        Index of beginning of end reference segment.
    ref_end : int
        Index of end of end reference segment. Default is None which
        will go to the end of the data array.
    seg_start : int
        What index to start computing the statistic. Default is 0 which
        will go to the beginning of the data array.

    Returns
    -------
    starts : numpy.array
        The start index of the first segment in the chain.
    ends : numpy.array
        The end index of the first segment in the chain.
    stats : numpy.array
        The Geweke convergence diagnostic statistic for the segment.
    """

    # lists to hold statistic and end index
    stats = []
    ends = []

    # get the beginning of all segments
    starts = numpy.arange(seg_start, end_idx, seg_stride)

    # get second segment of data at the end to compare
    x_end = x[ref_start:ref_end]

    # loop over all segments
    for start in starts:

        # find the end of the first segment
        x_start_end = int(start + seg_length)

        # get first segment
        x_start = x[start:x_start_end]

        # compute statistic
        stats.append((x_start.mean() - x_end.mean()) / numpy.sqrt(
            x_start.var() + x_end.var()))

        # store end of first segment
        ends.append(x_start_end)

    return numpy.array(starts), numpy.array(ends), numpy.array(stats)
