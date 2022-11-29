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

import lalframe, logging
import lal
import numpy
import math
import os.path, glob, time
import gwdatafind
import pycbc
from urllib.parse import urlparse
from pycbc.types import TimeSeries, zeros


# map LAL series types to corresponding functions and Numpy types
_fr_type_map = {
    lal.S_TYPE_CODE: [
        lalframe.FrStreamReadREAL4TimeSeries, numpy.float32,
        lal.CreateREAL4TimeSeries,
        lalframe.FrStreamGetREAL4TimeSeriesMetadata,
        lal.CreateREAL4Sequence,
        lalframe.FrameAddREAL4TimeSeriesProcData
    ],
    lal.D_TYPE_CODE: [
        lalframe.FrStreamReadREAL8TimeSeries, numpy.float64,
        lal.CreateREAL8TimeSeries,
        lalframe.FrStreamGetREAL8TimeSeriesMetadata,
        lal.CreateREAL8Sequence,
        lalframe.FrameAddREAL8TimeSeriesProcData
    ],
    lal.C_TYPE_CODE: [
        lalframe.FrStreamReadCOMPLEX8TimeSeries, numpy.complex64,
        lal.CreateCOMPLEX8TimeSeries,
        lalframe.FrStreamGetCOMPLEX8TimeSeriesMetadata,
        lal.CreateCOMPLEX8Sequence,
        lalframe.FrameAddCOMPLEX8TimeSeriesProcData
    ],
    lal.Z_TYPE_CODE: [
        lalframe.FrStreamReadCOMPLEX16TimeSeries, numpy.complex128,
        lal.CreateCOMPLEX16TimeSeries,
        lalframe.FrStreamGetCOMPLEX16TimeSeriesMetadata,
        lal.CreateCOMPLEX16Sequence,
        lalframe.FrameAddCOMPLEX16TimeSeriesProcData
    ],
    lal.U4_TYPE_CODE: [
        lalframe.FrStreamReadUINT4TimeSeries, numpy.uint32,
        lal.CreateUINT4TimeSeries,
        lalframe.FrStreamGetUINT4TimeSeriesMetadata,
        lal.CreateUINT4Sequence,
        lalframe.FrameAddUINT4TimeSeriesProcData
    ],
    lal.I4_TYPE_CODE: [
        lalframe.FrStreamReadINT4TimeSeries, numpy.int32,
        lal.CreateINT4TimeSeries,
        lalframe.FrStreamGetINT4TimeSeriesMetadata,
        lal.CreateINT4Sequence,
        lalframe.FrameAddINT4TimeSeriesProcData
    ],
}

def _read_channel(channel, stream, start, duration):
    """ Get channel using lalframe """
    channel_type = lalframe.FrStreamGetTimeSeriesType(channel, stream)
    read_func = _fr_type_map[channel_type][0]
    d_type = _fr_type_map[channel_type][1]
    data = read_func(stream, channel, start, duration, 0)
    return TimeSeries(data.data.data, delta_t=data.deltaT, epoch=start,
                      dtype=d_type)


def _is_gwf(file_path):
    """Test if a file is a frame file by checking if its contents begins with
    the magic string 'IGWD'."""
    try:
        with open(file_path, 'rb') as f:
            if f.read(4) == b'IGWD':
                return True
    except IOError:
        pass
    return False


def locations_to_cache(locations, latest=False):
    """ Return a cumulative cache file build from the list of locations

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
                raise ValueError('no frame or cache files found in ' + source)
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

def read_frame(location, channels, start_time=None,
               end_time=None, duration=None, check_integrity=False,
               sieve=None):
    """Read time series from frame data.

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
        logging.info("Using frames that match regexp: %s", sieve)
        lal.CacheSieve(cum_cache, 0, 0, None, None, sieve)
    if start_time is not None and end_time is not None:
        # Before sieving, check if this is sane. Otherwise it will fail later.
        if (int(math.ceil(end_time)) - int(start_time)) <= 0:
            raise ValueError("Negative or null duration")
        lal.CacheSieve(cum_cache, int(start_time), int(math.ceil(end_time)),
                       None, None, None)

    stream = lalframe.FrStreamCacheOpen(cum_cache)
    stream.mode = lalframe.FR_STREAM_VERBOSE_MODE

    if check_integrity:
        stream.mode = (stream.mode | lalframe.FR_STREAM_CHECKSUM_MODE)

    lalframe.FrSetMode(stream.mode, stream)

    # determine duration of data
    if type(channels) is list:
        first_channel = channels[0]
    else:
        first_channel = channels

    data_length = lalframe.FrStreamGetVectorLength(first_channel, stream)
    channel_type = lalframe.FrStreamGetTimeSeriesType(first_channel, stream)
    create_series_func = _fr_type_map[channel_type][2]
    get_series_metadata_func = _fr_type_map[channel_type][3]
    series = create_series_func(first_channel, stream.epoch, 0, 0,
                                lal.ADCCountUnit, 0)
    get_series_metadata_func(series, stream)
    data_duration = (data_length + 0.5) * series.deltaT

    if start_time is None:
        start_time = stream.epoch*1
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
    #if duration > data_duration:
    #    raise ValueError("Requested duration longer than available data")

    if type(channels) is list:
        all_data = []
        for channel in channels:
            channel_data = _read_channel(channel, stream, start_time, duration)
            lalframe.FrStreamSeek(stream, start_time)
            all_data.append(channel_data)
        return all_data
    else:
        return _read_channel(channels, stream, start_time, duration)

def datafind_connection(server=None):
    """ Return a connection to the datafind server

    Parameters
    -----------
    server : {SERVER:PORT, string}, optional
       A string representation of the server and port.
       The port may be ommitted.

    Returns
    --------
    connection
        The open connection to the datafind server.
    """
    return gwdatafind.connect(host=server)

def frame_paths(frame_type, start_time, end_time, server=None, url_type='file'):
    """Return the paths to a span of frame files

    Parameters
    ----------
    frame_type : string
        The string representation of the frame type (ex. 'H1_ER_C00_L1')
    start_time : int
        The start time that we need the frames to span.
    end_time : int
        The end time that we need the frames to span.
    server : {None, SERVER:PORT string}, optional
        Optional string to specify the datafind server to use. By default an
        attempt is made to use a local datafind server.
    url_type : string
        Returns only frame URLs with a particular scheme or head such
        as "file" or "gsiftp". Default is "file", which queries locally
        stored frames. Option can be disabled if set to None.
    Returns
    -------
    paths : list of paths
        The list of paths to the frame files.

    Examples
    --------
    >>> paths = frame_paths('H1_LDAS_C02_L2', 968995968, 968995968+2048)
    """
    site = frame_type[0]
    connection = datafind_connection(server)
    connection.find_times(site, frame_type,
                          gpsstart=start_time, gpsend=end_time)
    cache = connection.find_frame_urls(site, frame_type, start_time, end_time,urltype=url_type)
    return [urlparse(entry).path for entry in cache]

def query_and_read_frame(frame_type, channels, start_time, end_time,
                         sieve=None, check_integrity=False):
    """Read time series from frame data.

    Query for the locatin of physical frames matching the frame type. Return
    a time series containing the channel between the given start and end times.

    Parameters
    ----------
    frame_type : string
        The type of frame file that we are looking for.
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
    # Allows compatibility with our standard tools
    # We may want to place this into a higher level frame getting tool
    if frame_type == 'LOSC_STRAIN':
        from pycbc.frame.losc import read_strain_losc
        if not isinstance(channels, list):
            channels = [channels]
        data = [read_strain_losc(c[:2], start_time, end_time)
                for c in channels]
        return data if len(data) > 1 else data[0]
    if frame_type == 'LOSC':
        from pycbc.frame.losc import read_frame_losc
        return read_frame_losc(channels, start_time, end_time)

    logging.info('querying datafind server')
    paths = frame_paths(frame_type, int(start_time), int(numpy.ceil(end_time)))
    logging.info('found files: %s' % (' '.join(paths)))
    return read_frame(paths, channels,
                      start_time=start_time,
                      end_time=end_time,
                      sieve=sieve,
                      check_integrity=check_integrity)

__all__ = ['read_frame', 'frame_paths',
           'datafind_connection',
           'query_and_read_frame']

def write_frame(location, channels, timeseries):
    """Write a list of time series to a single frame file.

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
    frame = lalframe.FrameNew(epoch=gps_start_time, duration=duration,
                              project='', run=1, frnum=1,
                              detectorFlags=lal.LALDETECTORTYPE_ABSENT)

    for i,tseries in enumerate(timeseries):
        # get data type
        for seriestype in _fr_type_map.keys():
            if _fr_type_map[seriestype][1] == tseries.dtype:
                create_series_func = _fr_type_map[seriestype][2]
                create_sequence_func = _fr_type_map[seriestype][4]
                add_series_func = _fr_type_map[seriestype][5]
                break

        # add time series to frame
        series = create_series_func(channels[i], tseries.start_time,
                       0, tseries.delta_t, lal.ADCCountUnit,
                       len(tseries.numpy()))
        series.data = create_sequence_func(len(tseries.numpy()))
        series.data.data = tseries.numpy()
        add_series_func(frame, series)

    # write frame
    lalframe.FrameWrite(frame, location)

class DataBuffer(object):

    """A linear buffer that acts as a FILO for reading in frame data
    """

    def __init__(self, frame_src,
                 channel_name,
                 start_time,
                 max_buffer=2048,
                 force_update_cache=True,
                 increment_update_cache=None,
                 dtype=numpy.float64):
        """ Create a rolling buffer of frame data

        Parameters
        ---------
        frame_src: str of list of strings
            Strings that indicate where to read from files from. This can be a
        list of frame files, a glob, etc.
        channel_name: str
            Name of the channel to read from the frame files
        start_time:
            Time to start reading from.
        max_buffer: {int, 2048}, Optional
            Length of the buffer in seconds
        dtype: {dtype, numpy.float32}, Optional
            Data type to use for the interal buffer
        """
        self.frame_src = frame_src
        self.channel_name = channel_name
        self.read_pos = start_time
        self.force_update_cache = force_update_cache
        self.increment_update_cache = increment_update_cache
        self.detector = channel_name.split(':')[0]

        self.update_cache()
        self.channel_type, self.raw_sample_rate = self._retrieve_metadata(self.stream, self.channel_name)

        raw_size = self.raw_sample_rate * max_buffer
        self.raw_buffer = TimeSeries(zeros(raw_size, dtype=dtype),
                                     copy=False,
                                     epoch=start_time - max_buffer,
                                     delta_t=1.0/self.raw_sample_rate)

    def update_cache(self):
        """Reset the lal cache. This can be used to update the cache if the
        result may change due to more files being added to the filesystem,
        for example.
        """
        cache = locations_to_cache(self.frame_src, latest=True)
        stream = lalframe.FrStreamCacheOpen(cache)
        self.stream = stream

    @staticmethod
    def _retrieve_metadata(stream, channel_name):
        """Retrieve basic metadata by reading the first file in the cache

        Parameters
        ----------
        stream: lal stream object
            Stream containing a channel we want to learn about
        channel_name: str
            The name of the channel we want to know the dtype and sample rate of

        Returns
        -------
        channel_type: lal type enum
            Enum value which indicates the dtype of the channel
        sample_rate: int
            The sample rate of the data within this channel
        """
        lalframe.FrStreamGetVectorLength(channel_name, stream)
        channel_type = lalframe.FrStreamGetTimeSeriesType(channel_name, stream)
        create_series_func = _fr_type_map[channel_type][2]
        get_series_metadata_func = _fr_type_map[channel_type][3]
        series = create_series_func(channel_name, stream.epoch, 0, 0,
                            lal.ADCCountUnit, 0)
        get_series_metadata_func(series, stream)
        return channel_type, int(1.0/series.deltaT)

    def _read_frame(self, blocksize):
        """Try to read the block of data blocksize seconds long

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel

        Returns
        -------
        data: TimeSeries
            TimeSeries containg 'blocksize' seconds of frame data

        Raises
        ------
        RuntimeError:
            If data cannot be read for any reason
        """
        try:
            read_func = _fr_type_map[self.channel_type][0]
            dtype = _fr_type_map[self.channel_type][1]
            data = read_func(self.stream, self.channel_name,
                             self.read_pos, int(blocksize), 0)
            return TimeSeries(data.data.data, delta_t=data.deltaT,
                              epoch=self.read_pos,
                              dtype=dtype)
        except Exception:
            raise RuntimeError('Cannot read {0} frame data'.format(self.channel_name))

    def null_advance(self, blocksize):
        """Advance and insert zeros

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel
        """
        self.raw_buffer.roll(-int(blocksize * self.raw_sample_rate))
        self.read_pos += blocksize
        self.raw_buffer.start_time += blocksize

    def advance(self, blocksize):
        """Add blocksize seconds more to the buffer, push blocksize seconds
        from the beginning.

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel
        """
        ts = self._read_frame(blocksize)

        self.raw_buffer.roll(-len(ts))
        self.raw_buffer[-len(ts):] = ts[:]
        self.read_pos += blocksize
        self.raw_buffer.start_time += blocksize
        return ts

    def update_cache_by_increment(self, blocksize):
        """Update the internal cache by starting from the first frame
        and incrementing.

        Guess the next frame file name by incrementing from the first found
        one. This allows a pattern to be used for the GPS folder of the file,
        which is indicated by `GPSX` where x is the number of digits to use.

        Parameters
        ----------
        blocksize: int
            Number of seconds to increment the next frame file.
        """
        start = float(self.raw_buffer.end_time)
        end = float(start + blocksize)

        if not hasattr(self, 'dur'):
            fname = glob.glob(self.frame_src[0])[0]
            fname = os.path.splitext(os.path.basename(fname))[0].split('-')

            self.beg = '-'.join([fname[0], fname[1]])
            self.ref = int(fname[2])
            self.dur = int(fname[3])

        fstart = int(self.ref + numpy.floor((start - self.ref) / float(self.dur)) * self.dur)
        starts = numpy.arange(fstart, end, self.dur).astype(int)

        keys = []
        for s in starts:
            pattern = self.increment_update_cache
            if 'GPS' in pattern:
                n = int(pattern[int(pattern.index('GPS') + 3)])
                pattern = pattern.replace('GPS%s' % n, str(s)[0:n])

            name = f'{pattern}/{self.beg}-{s}-{self.dur}.gwf'
            # check that file actually exists, else abort now
            if not os.path.exists(name):
                raise RuntimeError

            keys.append(name)
        cache = locations_to_cache(keys)
        stream = lalframe.FrStreamCacheOpen(cache)
        self.stream = stream
        self.channel_type, self.raw_sample_rate = \
            self._retrieve_metadata(self.stream, self.channel_name)

    def attempt_advance(self, blocksize, timeout=10):
        """ Attempt to advance the frame buffer. Retry upon failure, except
        if the frame file is beyond the timeout limit.

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel
        timeout: {int, 10}, Optional
            Number of seconds before giving up on reading a frame

        Returns
        -------
        data: TimeSeries
            TimeSeries containg 'blocksize' seconds of frame data
        """
        if self.force_update_cache:
            self.update_cache()
        while True:
            try:
                if self.increment_update_cache:
                    self.update_cache_by_increment(blocksize)
                return DataBuffer.advance(self, blocksize)
            except RuntimeError:
                if pycbc.gps_now() > timeout + self.raw_buffer.end_time:
                    # The frame is not there and it should be by now,
                    # so we give up and treat it as zeros
                    DataBuffer.null_advance(self, blocksize)
                    return None
                # I am too early to give up on this frame,
                # so we should try again
                time.sleep(0.1)

class StatusBuffer(DataBuffer):

    """ Read state vector or DQ information from a frame file """

    def __init__(self, frame_src,
                       channel_name,
                       start_time,
                       max_buffer=2048,
                       valid_mask=3,
                       force_update_cache=False,
                       increment_update_cache=None,
                       valid_on_zero=False):
        """ Create a rolling buffer of status data from a frame

        Parameters
        ---------
        frame_src: str of list of strings
            Strings that indicate where to read from files from. This can be a
            list of frame files, a glob, etc.
        channel_name: str
            Name of the channel to read from the frame files
        start_time:
            Time to start reading from.
        max_buffer: {int, 2048}, Optional
            Length of the buffer in seconds
        valid_mask: {int, HOFT_OK | SCIENCE_INTENT}, Optional
            Set of flags that must be on to indicate valid frame data.
        valid_on_zero: bool
            If True, `valid_mask` is ignored and the status is considered
            "good" simply when the channel is zero.
        """
        DataBuffer.__init__(self, frame_src, channel_name, start_time,
                            max_buffer=max_buffer,
                            force_update_cache=force_update_cache,
                            increment_update_cache=increment_update_cache,
                            dtype=numpy.int32)

        self.valid_mask = valid_mask
        self.valid_on_zero = valid_on_zero

    def check_valid(self, values, flag=None):
        """Check if the data contains any non-valid status information

        Parameters
        ----------
        values: pycbc.types.Array
            Array of status information
        flag: str, optional
            Override the default valid mask with a user defined mask.

        Returns
        -------
        status: boolean
            Returns True if all of the status information if valid,
             False if any is not.
        """
        if self.valid_on_zero:
            valid = values.numpy() == 0
        else:
            if flag is None:
                flag = self.valid_mask
            valid = numpy.bitwise_and(values.numpy(), flag) == flag
        return bool(numpy.all(valid))

    def is_extent_valid(self, start_time, duration, flag=None):
        """Check if the duration contains any non-valid frames

        Parameters
        ----------
        start_time: int
            Beginning of the duration to check in gps seconds
        duration: int
            Number of seconds after the start_time to check
        flag: str, optional
            Override the default valid mask with a user defined mask.

        Returns
        -------
        status: boolean
            Returns True if all of the status information if valid,
            False if any is not.
        """
        sr = self.raw_buffer.sample_rate
        s = int((start_time - self.raw_buffer.start_time) * sr)
        e = s + int(duration * sr) + 1
        data = self.raw_buffer[s:e]
        return self.check_valid(data, flag=flag)

    def indices_of_flag(self, start_time, duration, times, padding=0):
        """ Return the indices of the times lying in the flagged region

        Parameters
        ----------
        start_time: int
            Beginning time to request for
        duration: int
            Number of seconds to check.
        padding: float
            Number of seconds to add around flag inactive times to be considered
        inactive as well.

        Returns
        -------
        indices: numpy.ndarray
            Array of indices marking the location of triggers within valid
        time.
        """
        from pycbc.events.veto import indices_outside_times
        sr = self.raw_buffer.sample_rate
        s = int((start_time - self.raw_buffer.start_time - padding) * sr) - 1
        e = s + int((duration + padding) * sr) + 1
        data = self.raw_buffer[s:e]
        stamps = data.sample_times.numpy()

        if self.valid_on_zero:
            invalid = data.numpy() != 0
        else:
            invalid = numpy.bitwise_and(data.numpy(), self.valid_mask) \
                    != self.valid_mask

        starts = stamps[invalid] - padding
        ends = starts + 1.0 / sr + padding * 2.0
        idx = indices_outside_times(times, starts, ends)
        return idx

    def advance(self, blocksize):
        """ Add blocksize seconds more to the buffer, push blocksize seconds
        from the beginning.

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel

        Returns
        -------
        status: boolean
            Returns True if all of the status information if valid,
            False if any is not.
        """
        try:
            if self.increment_update_cache:
                self.update_cache_by_increment(blocksize)
            ts = DataBuffer.advance(self, blocksize)
            return self.check_valid(ts)
        except RuntimeError:
            self.null_advance(blocksize)
            return False

class iDQBuffer(object):

    """ Read iDQ timeseries from a frame file """

    def __init__(self, frame_src,
                 idq_channel_name,
                 idq_status_channel_name,
                 idq_threshold,
                 start_time,
                 max_buffer=512,
                 force_update_cache=False,
                 increment_update_cache=None):
        """
        Parameters
        ----------
        frame_src: str of list of strings
            Strings that indicate where to read from files from. This can be a
            list of frame files, a glob, etc.
        idq_channel_name: str
            Name of the channel to read the iDQ statistic from
        idq_status_channel_name: str
            Name of the channel to read the iDQ status from
        idq_threshold: float
            Threshold which triggers a veto if iDQ channel falls below this threshold
        start_time:
            Time to start reading from.
        max_buffer: {int, 512}, Optional
            Length of the buffer in seconds
        force_update_cache: {boolean, True}, Optional
            Re-check the filesystem for frame files on every attempt to
            read more data.
        increment_update_cache: {str, None}, Optional
            Pattern to look for frame files in a GPS dependent directory. This
            is an alternate to the forced updated of the frame cache, and
            apptempts to predict the next frame file name without probing the
            filesystem.
        """
        self.threshold = idq_threshold
        self.idq = DataBuffer(frame_src, idq_channel_name, start_time,
                                max_buffer=max_buffer,
                                force_update_cache=force_update_cache,
                                increment_update_cache=increment_update_cache)
        self.idq_state = DataBuffer(frame_src, idq_status_channel_name, start_time,
                                        max_buffer=max_buffer,
                                        force_update_cache=force_update_cache,
                                        increment_update_cache=increment_update_cache)

    def indices_of_flag(self, start_time, duration, times, padding=0):
        """ Return the indices of the times lying in the flagged region

        Parameters
        ----------
        start_time: int
            Beginning time to request for
        duration: int
            Number of seconds to check.
        padding: float
            Number of seconds to add around flag inactive times to be considered
        inactive as well.

        Returns
        -------
        indices: numpy.ndarray
            Array of indices marking the location of triggers within valid
        time.
        """
        from pycbc.events.veto import indices_outside_times
        sr = self.idq.raw_buffer.sample_rate
        s = int((start_time - self.idq.raw_buffer.start_time - padding) * sr) - 1
        e = s + int((duration + padding) * sr) + 1
        idq_fap = self.idq.raw_buffer[s:e]
        stamps = idq_fap.sample_times.numpy()
        low_fap = idq_fap.numpy() <= self.threshold
        idq_valid = self.idq_state.raw_buffer[s:e]
        idq_valid = idq_valid.numpy().astype(bool)
        valid_low_fap = numpy.logical_and(idq_valid, low_fap)
        glitch_idx = numpy.flatnonzero(valid_low_fap)
        glitch_times = stamps[glitch_idx]
        starts = glitch_times - padding
        ends = starts + 1.0 / sr + padding * 2.0
        idx = indices_outside_times(times, starts, ends)
        return idx

    def advance(self, blocksize):
        """ Add blocksize seconds more to the buffer, push blocksize seconds
        from the beginning.

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read

        Returns
        -------
        status: boolean
            Returns True if advance is succesful,
            False if not.
        """
        idq_ts = self.idq.attempt_advance(blocksize)
        idq_state_ts = self.idq_state.attempt_advance(blocksize)
        return (idq_ts is not None) and (idq_state_ts is not None)

    def null_advance(self, blocksize):
        """Advance and insert zeros

        Parameters
        ----------
        blocksize: int
            The number of seconds to advance the buffers
        """
        self.idq.null_advance(blocksize)
        self.idq_state.null_advance(blocksize)
