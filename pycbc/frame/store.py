# Copyright (C) 2019 Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Generals
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
This modules contains functions for reading in data from hdf stores
"""
from __future__ import division
import h5py
import numpy
from pycbc.types import TimeSeries


def read_store(fname, channel, start_time, end_time):
    """ Read time series data from hdf store

    Parameters
    ----------
    fname: str
        Name of hdf store file
    channel: str
        Channel name to read
    start_time: int
        GPS time to start reading from
    end_time: int
        GPS time to end time series

    Returns
    -------
    ts: pycbc.types.TimeSeries
        Time series containing the requested data

    """
    fhandle = h5py.File(fname, 'r')
    if channel not in fhandle:
        raise ValueError('Could not find channel name {}'.format(channel))

    # Determine which segment data lies in (can only read contiguous data now)
    starts = fhandle[channel]['segments']['start'][:]
    ends = fhandle[channel]['segments']['end'][:]

    diff = start_time - starts
    loc = numpy.where(diff >= 0)[0]
    sidx = loc[diff[loc].argmin()]

    stime = starts[sidx]
    etime = ends[sidx]

    if etime < end_time:
        raise ValueError("Cannot read data segment past {}".format(etime))

    data = fhandle[channel][str(sidx)]
    sample_rate = len(data) / (etime - stime)

    start = int((start_time - stime) * sample_rate)
    end = int((end_time - stime) * sample_rate)
    return TimeSeries(data[start:end], delta_t=1.0/sample_rate,
                      epoch=start_time)
