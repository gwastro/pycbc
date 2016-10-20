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
import os.path, glob, time
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
    channel_type = lalframe.FrStreamGetTimeSeriesType(channel, stream)
    read_func = _fr_type_map[channel_type][0]
    d_type = _fr_type_map[channel_type][1]
    data = read_func(stream, channel, start, duration, 0)
    return TimeSeries(data.data.data, delta_t=data.deltaT, epoch=start, 
                      dtype=d_type)

def locations_to_cache(locations):
    """ Return a cumulative cache file build from the list of locations
    
    Parameters
    ----------
    locations : list
        A list of strings containing files, globs, or cache files used to build
    a combined lal cache file object.

    Returns
    -------
    cache : lal.Cache
        A cumulative lal cache object containing the files derived from the
    list of locations
    """
    cum_cache = lal.Cache()
    for source in locations:
        for file_path in glob.glob(source):
            dir_name, file_name = os.path.split(file_path)
            base_name, file_extension = os.path.splitext(file_name)

            if file_extension == ".lcf" or file_extension == ".cache":
                cache = lal.CacheImport(file_path)
            elif file_extension == ".gwf": 
                cache = lalframe.FrOpen(dir_name, file_name).cache
            else:
                raise TypeError("Invalid location name")
                
            cum_cache = lal.CacheMerge(cum_cache, cache)
    return cum_cache

def read_frame(location, channels, start_time=None, 
               end_time=None, duration=None, check_integrity=True):
    """Read time series from frame data.

    Using a the `location`, which can either be a frame file ".gwf" or a 
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
    data_duration = data_length * series.deltaT

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
    # import inside function to avoid adding M2Crypto
    # as a general PyCBC requirement
    import glue.datafind

    if server:
        datafind_server = server
    else:
        # Get the server name from the environment
        if os.environ.has_key("LIGO_DATAFIND_SERVER"):
            datafind_server = os.environ["LIGO_DATAFIND_SERVER"]
        else:
            err = "Trying to obtain the ligo datafind server url from "
            err += "the environment, ${LIGO_DATAFIND_SERVER}, but that "
            err += "variable is not populated."
            raise ValueError(err)

    # verify authentication options
    if not datafind_server.endswith("80"):
        cert_file, key_file = glue.datafind.find_credential()
    else:
        cert_file, key_file = None, None

    # Is a port specified in the server URL
    server, port = datafind_server.split(':',1)
    if port == "":
        port = None
    else:
        port = int(port)

    # Open connection to the datafind server
    if cert_file and key_file:
        connection = glue.datafind.GWDataFindHTTPSConnection(
                host=server, port=port, cert_file=cert_file, key_file=key_file)
    else:
        connection = glue.datafind.GWDataFindHTTPConnection(
                host=server, port=port)
    return connection
    
def frame_paths(frame_type, start_time, end_time, server=None):
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
    times = connection.find_times(site, frame_type, 
                                  gpsstart=start_time, 
                                  gpsend=end_time)
    cache = connection.find_frame_urls(site, frame_type, start_time, end_time)
    paths = [entry.path for entry in cache]
    return paths    
    
def query_and_read_frame(frame_type, channels, start_time, end_time):
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
    logging.info('querying datafind server')
    paths = frame_paths(frame_type, int(start_time), int(end_time))
    logging.info('found files: %s' % (' '.join(paths)))
    return read_frame(paths, channels, 
                      start_time=start_time, 
                      end_time=end_time)
    
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
    gps_start_times = set([series.start_time for series in timeseries])
    gps_end_times = set([series.end_time for series in timeseries])
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
        cache = locations_to_cache(self.frame_src)
        stream = lalframe.FrStreamCacheOpen(cache)
        self.stream = stream

    def _retrieve_metadata(self, stream, channel_name):
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
        data_length = lalframe.FrStreamGetVectorLength(channel_name, stream)
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
        except Exception as e:
            raise RuntimeError('Cannot read requested frame data') 

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
        starts = numpy.arange(fstart, end, self.dur).astype(numpy.int)
        
        keys = []
        for s in starts:
            pattern = self.increment_update_cache
            if 'GPS' in pattern:
                n = int(pattern[int(pattern.index('GPS') + 3)])
                pattern = pattern.replace('GPS%s' % n, str(s)[0:n])
                
            name = '%s/%s-%s-%s.gwf' % (pattern, self.beg, s, self.dur)
            # check that file actually exists, else abort now
            if not os.path.exists(name):
                logging.info("%s does not seem to exist yet" % name)
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
        
        try:
            if self.increment_update_cache:
                self.update_cache_by_increment(blocksize)

            return DataBuffer.advance(self, blocksize)

        except RuntimeError:
            if lal.GPSTimeNow() > timeout + self.raw_buffer.end_time:
                # The frame is not there and it should be by now, so we give up
                # and treat it as zeros
                logging.info('Frame missing, giving up...')
                DataBuffer.null_advance(self, blocksize)
                return None
            else:
                # I am too early to give up on this frame, so we should try again
                logging.info('Frame missing, waiting a bit more...')
                time.sleep(1)
                return self.attempt_advance(blocksize, timeout=timeout)

# Status flags for the calibration state vector 
HOFT_OK = 1
SCIENCE_INTENT = 2
SCIENCE_QUALITY = 4
HOFT_PROD = 8
FILTERS_OK = 16
NO_STOCH_HW_INJ = 32
NO_CBC_HW_INJ = 64
NO_BURST_HW_INJ = 128
NO_DETCHAR_HW_INJ = 256
KAPPA_A_OK = 512
KAPPA_PU_OK = 1024
KAPPA_TST_OK = 2048
KAPPA_C_OK = 4096
FCC_OK = 8192
NO_GAP = 16384    
NO_HWINJ = NO_STOCH_HW_INJ | NO_CBC_HW_INJ | \
           NO_BURST_HW_INJ | NO_DETCHAR_HW_INJ

class StatusBuffer(DataBuffer):

    """ Read state vector information from a frame file """

    def __init__(self, frame_src, 
                       channel_name,
                       start_time,
                       max_buffer=2048,
                       valid_mask=HOFT_OK | SCIENCE_INTENT,
                       force_update_cache=False,
                       increment_update_cache=None):
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
        """
        DataBuffer.__init__(self, frame_src, channel_name, start_time,
                            max_buffer=max_buffer,
                            force_update_cache=force_update_cache,
                            increment_update_cache=increment_update_cache,
                            dtype=numpy.int32)
        self.valid_mask = valid_mask

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
        if flag is None:
            flag = self.valid_mask

        if numpy.any(numpy.bitwise_and(values.numpy(), flag) != flag):
            return False
        else:
            return True

    def is_extent_valid(self, start_time, duration, flag=None):
        """Check if the duration contains any non-valid frames

        Parameters
        ----------
        start_time: int
            Begging of the duration to check in gps seconds
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
        e = s + int(duration * sr)
        data = self.raw_buffer[s:e]
        return self.check_valid(data, flag=flag)

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
        if self.increment_update_cache:
            self.update_cache_by_increment(blocksize)

        try:
            ts = DataBuffer.advance(self, blocksize)
            return self.check_valid(ts)
        except RuntimeError:
            return False
