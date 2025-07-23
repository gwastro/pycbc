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
import logging
import numpy

from pycbc.types import TimeSeries
from pycbc.io.hdf import HFile

logger = logging.getLogger('pycbc.frame.store')


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
    fhandle = HFile(fname, 'r')

    starts = fhandle['meta']['GPSstart'][()]
    duration = fhandle['meta']['Duration'][()]
    ends = starts + duration

    if starts > start_time:
        raise ValueError("Cannot read data segment before {}".format(starts))

    if ends < end_time:
        raise ValueError("Cannot read data segment past {}".format(ends))

    data = fhandle['strain']['Strain'][:]
    sample_rate = len(data) / (ends - starts)

    start_inx = int((start_time - starts) * sample_rate)
    end_inx = int((end_time - starts) * sample_rate)
    return TimeSeries(data[start_inx:end_inx], delta_t=1.0/sample_rate,
                      epoch=start_time)
