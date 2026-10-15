# Copyright (C) 2014 Andrew Miller, Alex Nitz, Tito Dal Canton, Christopher M. Biwer
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
"""
Read and write GWF frame files using lalframe.

This is the only module in pycbc.frame that uses lalframe. It provides the
same functions as `pycbc.frame.gwosc_hdf` (is_format, span, channel_info and
read), each acting on a single file.
"""

import os.path
from pathlib import Path

import lal
import lalframe
import numpy as np

# map LAL series types to corresponding functions and Numpy types
_fr_type_map = {
    lal.S_TYPE_CODE: [
        lalframe.FrStreamReadREAL4TimeSeries,
        np.float32,
        lal.CreateREAL4TimeSeries,
        lalframe.FrStreamGetREAL4TimeSeriesMetadata,
        lal.CreateREAL4Sequence,
        lalframe.FrameAddREAL4TimeSeriesProcData,
    ],
    lal.D_TYPE_CODE: [
        lalframe.FrStreamReadREAL8TimeSeries,
        np.float64,
        lal.CreateREAL8TimeSeries,
        lalframe.FrStreamGetREAL8TimeSeriesMetadata,
        lal.CreateREAL8Sequence,
        lalframe.FrameAddREAL8TimeSeriesProcData,
    ],
    lal.C_TYPE_CODE: [
        lalframe.FrStreamReadCOMPLEX8TimeSeries,
        np.complex64,
        lal.CreateCOMPLEX8TimeSeries,
        lalframe.FrStreamGetCOMPLEX8TimeSeriesMetadata,
        lal.CreateCOMPLEX8Sequence,
        lalframe.FrameAddCOMPLEX8TimeSeriesProcData,
    ],
    lal.Z_TYPE_CODE: [
        lalframe.FrStreamReadCOMPLEX16TimeSeries,
        np.complex128,
        lal.CreateCOMPLEX16TimeSeries,
        lalframe.FrStreamGetCOMPLEX16TimeSeriesMetadata,
        lal.CreateCOMPLEX16Sequence,
        lalframe.FrameAddCOMPLEX16TimeSeriesProcData,
    ],
    lal.U4_TYPE_CODE: [
        lalframe.FrStreamReadUINT4TimeSeries,
        np.uint32,
        lal.CreateUINT4TimeSeries,
        lalframe.FrStreamGetUINT4TimeSeriesMetadata,
        lal.CreateUINT4Sequence,
        lalframe.FrameAddUINT4TimeSeriesProcData,
    ],
    lal.I4_TYPE_CODE: [
        lalframe.FrStreamReadINT4TimeSeries,
        np.int32,
        lal.CreateINT4TimeSeries,
        lalframe.FrStreamGetINT4TimeSeriesMetadata,
        lal.CreateINT4Sequence,
        lalframe.FrameAddINT4TimeSeriesProcData,
    ],
}


def is_format(path):
    """
    Test if a file is a frame file by checking if its contents begins with
    the magic string 'IGWD'.
    """
    try:
        with Path(path).open("rb") as f:
            return f.read(4) == b"IGWD"
    except OSError:
        return False


def _open(path, check_integrity=False):
    """Open a lalframe stream reading the single file `path`."""
    dir_name, file_name = os.path.split(str(Path(path).resolve()))
    stream = lalframe.FrStreamOpen(dir_name, file_name)
    mode = lalframe.FR_STREAM_VERBOSE_MODE
    if check_integrity:
        mode |= lalframe.FR_STREAM_CHECKSUM_MODE
    lalframe.FrStreamSetMode(stream, mode)
    return stream


def span(path):
    """Return the GPS start time and duration of the data in a file."""
    dir_name, file_name = os.path.split(str(Path(path).resolve()))
    entry = lalframe.FrOpen(dir_name, file_name).cache.list
    return float(entry.t0), float(entry.dt)


def channel_info(path, channel):
    """Return the numpy dtype and sample rate of a channel in a file."""
    stream = _open(path)
    lalframe.FrStreamGetVectorLength(channel, stream)
    channel_type = lalframe.FrStreamGetTimeSeriesType(channel, stream)
    create_series_func = _fr_type_map[channel_type][2]
    get_series_metadata_func = _fr_type_map[channel_type][3]
    series = create_series_func(channel, stream.epoch, 0, 0, lal.ADCCountUnit, 0)
    get_series_metadata_func(series, stream)
    return np.dtype(_fr_type_map[channel_type][1]), 1.0 / series.deltaT


def read(path, channel, start, end, check_integrity=False):
    """
    Read the samples of a channel between GPS times `start` and `end`,
    which must lie within the file.

    Returns
    -------
    data : np.ndarray
    epoch : float
        GPS time of the first sample.
    delta_t : float
        Sample spacing in seconds.

    """
    stream = _open(path, check_integrity=check_integrity)
    channel_type = lalframe.FrStreamGetTimeSeriesType(channel, stream)
    read_func = _fr_type_map[channel_type][0]
    series = read_func(stream, channel, lal.LIGOTimeGPS(start), end - start, 0)
    return series.data.data, float(series.epoch), series.deltaT


def write(location, channels, timeseries, start_time, duration):
    """
    Write a list of time series to a single frame file. See
    `pycbc.frame.write_frame`.
    """
    # create frame
    frame = lalframe.FrameNew(
        epoch=start_time,
        duration=duration,
        project="",
        run=1,
        frnum=1,
        detectorFlags=lal.LALDETECTORTYPE_ABSENT,
    )

    for i, tseries in enumerate(timeseries):
        # get data type
        for seriestype in _fr_type_map:
            if _fr_type_map[seriestype][1] == tseries.dtype:
                create_series_func = _fr_type_map[seriestype][2]
                create_sequence_func = _fr_type_map[seriestype][4]
                add_series_func = _fr_type_map[seriestype][5]
                break

        # add time series to frame
        series = create_series_func(
            channels[i],
            tseries.start_time,
            0,
            tseries.delta_t,
            lal.ADCCountUnit,
            len(tseries.numpy()),
        )
        series.data = create_sequence_func(len(tseries.numpy()))
        series.data.data = tseries.numpy()
        add_series_func(frame, series)

    # write frame
    lalframe.FrameWrite(frame, location)
