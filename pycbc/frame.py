# Copyright (C) 2012 Andrew Miller, Alex Nitz, Tito Dal Canton
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
This modules contains functions for reading in data from frame files or caches
"""
import lalframe
import lal
import numpy
import os.path
from pycbc.types import TimeSeries
import copy


# map frame vector types to corresponding functions and Numpy types
_fr_type_map = {
    lalframe.LAL_FRAMEU_FR_VECT_2S: [
        lalframe.FrReadINT2TimeSeries, numpy.int16,
        lal.CreateINT2TimeSeries,
        lalframe.FrGetINT2TimeSeriesMetadata
    ],
    lalframe.LAL_FRAMEU_FR_VECT_4S: [
        lalframe.FrReadINT4TimeSeries, numpy.int32,
        lal.CreateINT4TimeSeries,
        lalframe.FrGetINT4TimeSeriesMetadata
    ],
    lalframe.LAL_FRAMEU_FR_VECT_8S: [
        lalframe.FrReadINT8TimeSeries, numpy.int64,
        lal.CreateINT8TimeSeries,
        lalframe.FrGetINT8TimeSeriesMetadata
    ],
    lalframe.LAL_FRAMEU_FR_VECT_4U: [
        lalframe.FrReadREAL4TimeSeries, numpy.float32,
        lal.CreateREAL4TimeSeries,
        lalframe.FrGetREAL4TimeSeriesMetadata
    ],
    lalframe.LAL_FRAMEU_FR_VECT_8U: [
        lalframe.FrReadREAL8TimeSeries, numpy.float64,
        lal.CreateREAL8TimeSeries,
        lalframe.FrGetREAL8TimeSeriesMetadata
    ],
    lalframe.LAL_FRAMEU_FR_VECT_8C: [
        lalframe.FrReadCOMPLEX8TimeSeries, numpy.complex64,
        lal.CreateCOMPLEX8TimeSeries,
        lalframe.FrGetCOMPLEX8TimeSeriesMetadata
    ],
    lalframe.LAL_FRAMEU_FR_VECT_16C: [
        lalframe.FrReadCOMPLEX16TimeSeries, numpy.complex128,
        lal.CreateCOMPLEX16TimeSeries,
        lalframe.FrGetCOMPLEX16TimeSeriesMetadata
    ],
}

def _read_channel(channel, stream, start, duration):
    channel_type = lalframe.FrGetTimeSeriesType(channel, stream)
    read_func = _fr_type_map[channel_type][0]
    d_type = _fr_type_map[channel_type][1]
    data = read_func(stream, channel, start, duration, 0)
    return TimeSeries(data.data.data, delta_t=data.deltaT, epoch=start, dtype=d_type)

def read_frame(location, channels, start_time=None, end_time=None, duration=None):
    """Read time series from frame data.

    Using a the `location`, which can either be a frame file ".gwf" or a 
    frame cache ".gwf", read in the data for the given channel(s) and output
    as a TimeSeries or list of TimeSeries. 

    Parameters
    ----------
    location : string
        Either a frame filename (can include pattern) or the cache filename.  
    channels : string or list of strings
        Either a string that contains the channel name or a list of channel name
        strings.
    start_time : {None, LIGOTimeGPS}, optional
        The gps start time of the time series. Defaults to reading from the 
        beginning of the available frame(s). 
    end_time : {None, LIGOTimeGPS}, optional
        The gps end time of the time series. Defaults to the end of the frame(s).
        Note, this argument is incompatible with `duration`.
    duration : {None, float}, optional
        The amount of data to read in seconds. Note, this argument is incompatible
        with `end`.

    Returns
    -------
    Frame Data: TimeSeries or list of TimeSeries
        A TimeSeries or a list of TimeSeries, corresponding to the data from the
        frame file/cache for a given channel or channels. 
    """

    if end_time and duration:
        raise ValueError("end time and duration are mutually exclusive")

    dir_name, file_name = os.path.split(location)
    base_name, file_extension = os.path.splitext(file_name)

    if file_extension == ".lcf":
        cache = lalframe.FrImportCache(location)
        stream = lalframe.FrCacheOpen(cache)
    elif file_extension == ".gwf": 
        stream = lalframe.FrOpen(dir_name, file_name)
    else:
        raise TypeError("Invalid location name")
        
    stream.mode = lalframe.LAL_FR_VERBOSE_MODE
    lalframe.FrSetMode(stream, stream.mode | lalframe.LAL_FR_CHECKSUM_MODE)

    # determine duration of data
    if type(channels) is list:
        first_channel = channels[0]
    else:
        first_channel = channels
    data_length = lalframe.FrGetVectorLength(first_channel, stream)
    channel_type = lalframe.FrGetTimeSeriesType(first_channel, stream)
    create_series_func = _fr_type_map[channel_type][2]
    get_series_metadata_func = _fr_type_map[channel_type][3]
    series = create_series_func(first_channel, stream.epoch, 0, 0,
                                lal.lalADCCountUnit, 0)
    get_series_metadata_func(series, stream)
    data_duration = data_length * series.deltaT

    if start_time is None:
        start_time = stream.epoch*1
    if end_time is None:
        end_time = start_time + data_duration

    if start_time is not lal.LIGOTimeGPS:
        start_time = lal.LIGOTimeGPS(start_time)
    if end_time is not lal.LIGOTimeGPS:
        end_time = lal.LIGOTimeGPS(end_time)

    if duration is None:
        duration = float(end_time - start_time)
    else:
        duration = float(duration)

    # lalframe behaves dangerously with invalid duration so catch it here
    if duration <= 0:
        raise ValueError("Negative or null duration")
    if duration > data_duration:
        raise ValueError("Requested duration longer than available data")

    if type(channels) is list:
        all_data = []
        for channel in channels:
            channel_data = _read_channel(channel, stream, start_time, duration)
            lalframe.FrSeek(stream, start_time)
            all_data.append(channel_data)
        return all_data
    else:
        return _read_channel(channels, stream, start_time, duration)
