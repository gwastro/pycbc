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

import logging
import re
import warnings
from pathlib import Path
from urllib.parse import urlparse

import numpy
from gwdatafind import find_urls as find_frame_urls

from pycbc.frame import gwf, gwosc_hdf
from pycbc.frame.cache import file_format, locations_to_cache
from pycbc.types import TimeSeries

logger = logging.getLogger(__name__)


def _read_channel(cache, channel, start, end, check_integrity=False):
    """
    Read one channel from the files in `cache`, which must hold contiguous
    data from `start` to `end`. If `end` is None, read up to the end of the
    contiguous data.
    """
    start = float(start)
    end = None if end is None else float(end)
    pieces = []
    for entry in cache.sieve(start=start, end=end):
        data, epoch, delta_t = file_format(entry.path).read(
            entry.path,
            channel,
            max(start, entry.start),
            entry.end if end is None else min(end, entry.end),
            check_integrity=check_integrity,
        )
        if pieces:
            prev_data, prev_epoch, prev_delta_t = pieces[-1]
            expected = prev_epoch + len(prev_data) * prev_delta_t
            if abs(epoch - expected) > delta_t / 2 or delta_t != prev_delta_t:
                if end is None:
                    break
                raise ValueError(  # noqa: TRY003
                    f"{channel} has a gap, overlap or change of "
                    f"sample rate at GPS time {expected}"
                )
        pieces.append((data, epoch, delta_t))

    if not pieces:
        raise ValueError(f"No data found for {channel} at GPS time {start}")  # noqa: TRY003
    data, epoch, delta_t = pieces[0]
    last_data, last_epoch, _ = pieces[-1]
    if epoch - start >= delta_t:
        raise ValueError(  # noqa: TRY003
            f"{channel} data starts at {epoch}, after the requested start time {start}"
        )
    if end is not None and end - (last_epoch + len(last_data) * delta_t) >= delta_t:
        raise ValueError(f"{channel} data ends before the requested end time {end}")  # noqa: TRY003
    return TimeSeries(
        numpy.concatenate([p[0] for p in pieces]), delta_t=delta_t, epoch=epoch
    )


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
        The gps end time of the time series. Defaults to the end of the
        contiguous data following `start_time`.
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

    cache = locations_to_cache(locations)
    if sieve:
        logger.info("Using frames that match regexp: %s", sieve)
        cache = cache.sieve(regex=sieve)
    if not cache:
        raise ValueError(f"No frame files found in {location}")  # noqa: TRY003

    if start_time is None:
        start_time = cache[0].start
    if duration is not None:
        end_time = float(start_time) + float(duration)
    if end_time is not None and end_time <= start_time:
        raise ValueError("Negative or null duration")

    if type(channels) is list:
        return [
            _read_channel(
                cache, channel, start_time, end_time, check_integrity=check_integrity
            )
            for channel in channels
        ]
    return _read_channel(
        cache, channels, start_time, end_time, check_integrity=check_integrity
    )


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

    if location.endswith((".hdf", ".hdf5", ".h5")) or (
        Path(location).exists() and gwosc_hdf.is_format(location)
    ):
        gwosc_hdf.write(location, channels, timeseries, gps_start_time, duration)
    else:
        gwf.write(location, channels, timeseries, gps_start_time, duration)


__all__ = [
    "frame_paths",
    "locations_to_cache",
    "query_and_read_frame",
    "read_frame",
    "write_frame",
]
