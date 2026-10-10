# Copyright (C) 2014 Andrew Miller, Alex Nitz, Tito Dal Canton, Christopher M. Biwer
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

import glob
import logging
import math
import os.path
import re
import warnings
from urllib.parse import urlparse

import lal
import lalframe
import numpy
from gwdatafind import find_urls as find_frame_urls

from pycbc.types import TimeSeries

logger = logging.getLogger("pycbc.frame.frame")

# map LAL series types to corresponding functions and Numpy types
_fr_type_map = {
    lal.S_TYPE_CODE: [
        lalframe.FrStreamReadREAL4TimeSeries,
        numpy.float32,
        lal.CreateREAL4TimeSeries,
        lalframe.FrStreamGetREAL4TimeSeriesMetadata,
        lal.CreateREAL4Sequence,
        lalframe.FrameAddREAL4TimeSeriesProcData,
    ],
    lal.D_TYPE_CODE: [
        lalframe.FrStreamReadREAL8TimeSeries,
        numpy.float64,
        lal.CreateREAL8TimeSeries,
        lalframe.FrStreamGetREAL8TimeSeriesMetadata,
        lal.CreateREAL8Sequence,
        lalframe.FrameAddREAL8TimeSeriesProcData,
    ],
    lal.C_TYPE_CODE: [
        lalframe.FrStreamReadCOMPLEX8TimeSeries,
        numpy.complex64,
        lal.CreateCOMPLEX8TimeSeries,
        lalframe.FrStreamGetCOMPLEX8TimeSeriesMetadata,
        lal.CreateCOMPLEX8Sequence,
        lalframe.FrameAddCOMPLEX8TimeSeriesProcData,
    ],
    lal.Z_TYPE_CODE: [
        lalframe.FrStreamReadCOMPLEX16TimeSeries,
        numpy.complex128,
        lal.CreateCOMPLEX16TimeSeries,
        lalframe.FrStreamGetCOMPLEX16TimeSeriesMetadata,
        lal.CreateCOMPLEX16Sequence,
        lalframe.FrameAddCOMPLEX16TimeSeriesProcData,
    ],
    lal.U4_TYPE_CODE: [
        lalframe.FrStreamReadUINT4TimeSeries,
        numpy.uint32,
        lal.CreateUINT4TimeSeries,
        lalframe.FrStreamGetUINT4TimeSeriesMetadata,
        lal.CreateUINT4Sequence,
        lalframe.FrameAddUINT4TimeSeriesProcData,
    ],
    lal.I4_TYPE_CODE: [
        lalframe.FrStreamReadINT4TimeSeries,
        numpy.int32,
        lal.CreateINT4TimeSeries,
        lalframe.FrStreamGetINT4TimeSeriesMetadata,
        lal.CreateINT4Sequence,
        lalframe.FrameAddINT4TimeSeriesProcData,
    ],
}


def _read_channel(channel, stream, start, duration):
    """Get channel using lalframe"""
    channel_type = lalframe.FrStreamGetTimeSeriesType(channel, stream)
    read_func = _fr_type_map[channel_type][0]
    d_type = _fr_type_map[channel_type][1]
    data = read_func(stream, channel, start, duration, 0)
    return TimeSeries(data.data.data, delta_t=data.deltaT, epoch=start, dtype=d_type)


def _is_gwf(file_path):
    """
    Test if a file is a frame file by checking if its contents begins with
    the magic string 'IGWD'.
    """
    try:
        with open(file_path, "rb") as f:
            if f.read(4) == b"IGWD":
                return True
    except OSError:
        pass
    return False


def locations_to_cache(locations, latest=False):
    """
    Return a cumulative cache file build from the list of locations

    Parameters
    ----------
    locations : list
        A list of strings containing files, globs, or cache files used to
        build a combined lal cache file object.
    latest : Optional, {False, Boolean}
        Only return a cache with the most recent frame in the locations.
        If false, all results are returned.

    Returns
    -------
    cache : lal.Cache
        A cumulative lal cache object containing the files derived from the
        list of locations.

    """
    cum_cache = lal.Cache()
    for source in locations:
        flist = glob.glob(source)
        if latest:

            def relaxed_getctime(fn):
                # when building a cache from a directory of temporary
                # low-latency frames, files might disappear between
                # the glob() and getctime() calls
                try:
                    return os.path.getctime(fn)
                except OSError:
                    return 0

            if not flist:
                raise ValueError("no frame or cache files found in " + source)
            flist = [max(flist, key=relaxed_getctime)]

        for file_path in flist:
            dir_name, file_name = os.path.split(file_path)
            _, file_extension = os.path.splitext(file_name)

            if file_extension in [".lcf", ".cache"]:
                cache = lal.CacheImport(file_path)
            elif file_extension == ".gwf" or _is_gwf(file_path):
                cache = lalframe.FrOpen(str(dir_name), str(file_name)).cache
            else:
                raise TypeError("Invalid location name")

            cum_cache = lal.CacheMerge(cum_cache, cache)
    return cum_cache


def read_frame(
    location,
    channels,
    start_time=None,
    end_time=None,
    duration=None,
    check_integrity=False,
    sieve=None,
):
    """
    Read time series from frame data.

    Using the `location`, which can either be a frame file ".gwf" or a
    frame cache ".gwf", read in the data for the given channel(s) and output
    as a TimeSeries or list of TimeSeries.

    Parameters
    ----------
    location : string
        A source of gravitational wave frames. Either a frame filename
        (can include pattern), a list of frame files, or frame cache file.
    channels : string or list of strings
        Either a string that contains the channel name or a list of channel
        name strings.
    start_time : {None, LIGOTimeGPS}, optional
        The gps start time of the time series. Defaults to reading from the
        beginning of the available frame(s).
    end_time : {None, LIGOTimeGPS}, optional
        The gps end time of the time series. Defaults to the end of the frame.
        Note, this argument is incompatible with `duration`.
    duration : {None, float}, optional
        The amount of data to read in seconds. Note, this argument is
        incompatible with `end`.
    check_integrity : {True, bool}, optional
        Test the frame files for internal integrity.
    sieve : string, optional
        Selects only frames where the frame URL matches the regular
        expression sieve

    Returns
    -------
    Frame Data: TimeSeries or list of TimeSeries
        A TimeSeries or a list of TimeSeries, corresponding to the data from
        the frame file/cache for a given channel or channels.

    """
    if end_time and duration:
        raise ValueError("end time and duration are mutually exclusive")

    if type(location) is list:
        locations = location
    else:
        locations = [location]

    cum_cache = locations_to_cache(locations)
    if sieve:
        logger.info("Using frames that match regexp: %s", sieve)
        lal.CacheSieve(cum_cache, 0, 0, None, None, sieve)
    if start_time is not None and end_time is not None:
        # Before sieving, check if this is sane. Otherwise it will fail later.
        if (int(math.ceil(end_time)) - int(start_time)) <= 0:
            raise ValueError("Negative or null duration")
        lal.CacheSieve(
            cum_cache, int(start_time), int(math.ceil(end_time)), None, None, None
        )

    stream = lalframe.FrStreamCacheOpen(cum_cache)
    stream.mode = lalframe.FR_STREAM_VERBOSE_MODE

    if check_integrity:
        stream.mode = stream.mode | lalframe.FR_STREAM_CHECKSUM_MODE

    lalframe.FrStreamSetMode(stream, stream.mode)

    # determine duration of data
    if type(channels) is list:
        first_channel = channels[0]
    else:
        first_channel = channels

    data_length = lalframe.FrStreamGetVectorLength(first_channel, stream)
    channel_type = lalframe.FrStreamGetTimeSeriesType(first_channel, stream)
    create_series_func = _fr_type_map[channel_type][2]
    get_series_metadata_func = _fr_type_map[channel_type][3]
    series = create_series_func(first_channel, stream.epoch, 0, 0, lal.ADCCountUnit, 0)
    get_series_metadata_func(series, stream)
    data_duration = (data_length + 0.5) * series.deltaT

    if start_time is None:
        start_time = stream.epoch * 1
    if end_time is None:
        end_time = start_time + data_duration

    if type(start_time) is not lal.LIGOTimeGPS:
        start_time = lal.LIGOTimeGPS(start_time)
    if type(end_time) is not lal.LIGOTimeGPS:
        end_time = lal.LIGOTimeGPS(end_time)

    if duration is None:
        duration = float(end_time - start_time)
    else:
        duration = float(duration)

    # lalframe behaves dangerously with invalid duration so catch it here
    if duration <= 0:
        raise ValueError("Negative or null duration")
    # if duration > data_duration:
    #    raise ValueError("Requested duration longer than available data")

    if type(channels) is list:
        all_data = []
        for channel in channels:
            channel_data = _read_channel(channel, stream, start_time, duration)
            lalframe.FrStreamSeek(stream, start_time)
            all_data.append(channel_data)
        return all_data
    return _read_channel(channels, stream, start_time, duration)


def frame_paths(
    frame_type, start_time, end_time, server=None, url_type="file", site=None
):
    """
    Return the paths to a span of frame files.

    Parameters
    ----------
    frame_type : string
        The string representation of the frame type (ex. 'H1_ER_C00_L1').
    start_time : int
        The start time that we need the frames to span.
    end_time : int
        The end time that we need the frames to span.
    server : {None, SERVER:PORT string}, optional
        Optional string to specify the datafind server to use. By default an
        attempt is made to use a local datafind server.
    url_type : string
        Returns only frame URLs with a particular scheme or head such
        as "file" or "https". Default is "file", which queries locally
        stored frames. Option can be disabled if set to None.
    site : string, optional
        One-letter string specifying which site you want data from (H, L, V,
        etc).  If not given, the site is assumed to be the first letter of
        `frame_type`, which is usually (but not always) a safe assumption.

    Returns
    -------
    paths : list of paths
        The list of paths to the frame files.

    Examples
    --------
    >>> paths = frame_paths('H1_LDAS_C02_L2', 968995968, 968995968+2048)

    """
    if site is None:
        # this case is tolerated for backward compatibility
        site = frame_type[0]
        warnings.warn(
            f"Guessing site {site} from frame type {frame_type}", DeprecationWarning
        )
    cache = find_frame_urls(
        site, frame_type, start_time, end_time, urltype=url_type, host=server
    )
    return [urlparse(entry).path for entry in cache]


def get_site_from_type_or_channel(frame_type, channels):
    """
    Determine the site for querying gwdatafind (H, L, V, etc) based on
    substrings of the frame type and channel(s).

    The type should begin with S: or SN:, in which case S is taken as the
    site.  Otherwise, the same is done with the channel (with the first
    channel if more than one are given). If that also fails, the site is
    taken to be the first letter of the frame type, which is usually
    (but not always) a correct assumption.

    Parameters
    ----------
    frame_type : string
        The frame type, ideally prefixed by the site indicator.
    channels : string or list of strings
        The channel name or names.

    Returns
    -------
    site : string
        The site letter.
    frame_type : string
        The frame type with the site prefix (if any) removed.

    """
    site_re = "^([^:])[^:]?:"
    m = re.match(site_re, frame_type)
    if m:
        return m.groups(1)[0], frame_type[m.end() :]
    chan = channels
    if isinstance(chan, list):
        chan = channels[0]
    m = re.match(site_re, chan)
    if m:
        return m.groups(1)[0], frame_type
    warnings.warn(
        f"Guessing site {frame_type[0]} from frame type {frame_type}",
        DeprecationWarning,
    )
    return frame_type[0], frame_type


def query_and_read_frame(
    frame_type, channels, start_time, end_time, sieve=None, check_integrity=False
):
    """
    Read time series from frame data.

    Query for the location of physical frames matching the frame type. Return
    a time series containing the channel between the given start and end times.

    Parameters
    ----------
    frame_type : string
        The type of frame file that we are looking for. The string should begin
        with S: or SN:, in which case S is taken as the site to query. If this
        is not the case, the site will be guessed from the channel name or from
        the type in a different way, which may not work.
    channels : string or list of strings
        Either a string that contains the channel name or a list of channel
        name strings.
    start_time : LIGOTimeGPS or int
        The gps start time of the time series. Defaults to reading from the
        beginning of the available frame(s).
    end_time : LIGOTimeGPS or int
        The gps end time of the time series. Defaults to the end of the frame.
    sieve : string, optional
        Selects only frames where the frame URL matches the regular
        expression sieve
    check_integrity : boolean
        Do an expensive checksum of the file before returning.

    Returns
    -------
    Frame Data: TimeSeries or list of TimeSeries
        A TimeSeries or a list of TimeSeries, corresponding to the data from
        the frame file/cache for a given channel or channels.

    Examples
    --------
    >>> ts = query_and_read_frame('H1_LDAS_C02_L2', 'H1:LDAS-STRAIN',
    >>>                               968995968, 968995968+2048)

    """
    site, frame_type = get_site_from_type_or_channel(frame_type, channels)

    # Allows compatibility with our standard tools
    # We may want to place this into a higher level frame getting tool
    if frame_type in ["LOSC_STRAIN", "GWOSC_STRAIN"]:
        from pycbc.frame.gwosc import read_strain_gwosc

        if not isinstance(channels, list):
            channels = [channels]
        data = [read_strain_gwosc(c[:2], start_time, end_time) for c in channels]
        return data if len(data) > 1 else data[0]
    if frame_type in ["LOSC", "GWOSC"]:
        from pycbc.frame.gwosc import read_frame_gwosc

        return read_frame_gwosc(channels, start_time, end_time)

    logger.info("Querying datafind server")
    paths = frame_paths(
        frame_type, int(start_time), int(numpy.ceil(end_time)), site=site
    )
    logger.info("Found frame file paths: %s", " ".join(paths))
    return read_frame(
        paths,
        channels,
        start_time=start_time,
        end_time=end_time,
        sieve=sieve,
        check_integrity=check_integrity,
    )


def write_frame(location, channels, timeseries):
    """
    Write a list of time series to a single frame file.

    Parameters
    ----------
    location : string
        A frame filename.
    channels : string or list of strings
        Either a string that contains the channel name or a list of channel
        name strings.
    timeseries: TimeSeries
        A TimeSeries or list of TimeSeries, corresponding to the data to be
        written to the frame file for a given channel.

    """
    # check if a single channel or a list of channels
    if type(channels) is list and type(timeseries) is list:
        channels = channels
        timeseries = timeseries
    else:
        channels = [channels]
        timeseries = [timeseries]

    # check that timeseries have the same start and end time
    gps_start_times = {series.start_time for series in timeseries}
    gps_end_times = {series.end_time for series in timeseries}
    if len(gps_start_times) != 1 or len(gps_end_times) != 1:
        raise ValueError("Start and end times of TimeSeries must be identical.")

    # check that start, end time, and duration are integers
    gps_start_time = gps_start_times.pop()
    gps_end_time = gps_end_times.pop()
    duration = int(gps_end_time - gps_start_time)
    if gps_start_time % 1 or gps_end_time % 1:
        raise ValueError("Start and end times of TimeSeries must be integer seconds.")

    # create frame
    frame = lalframe.FrameNew(
        epoch=gps_start_time,
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


__all__ = [
    "frame_paths",
    "locations_to_cache",
    "query_and_read_frame",
    "read_frame",
    "write_frame",
]
