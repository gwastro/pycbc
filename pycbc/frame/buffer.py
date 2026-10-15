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
Rolling buffers of frame data, used by low-latency analyses.
"""

import glob
import time
from pathlib import Path

import numpy as np

import pycbc
from pycbc.events.veto import indices_outside_times, indices_within_times
from pycbc.frame.cache import file_format, locations_to_cache
from pycbc.frame.frame import _read_channel
from pycbc.types import TimeSeries, zeros


class DataBuffer:
    """A linear buffer that acts as a FILO for reading in frame data"""

    def __init__(
        self,
        frame_src,
        channel_name,
        start_time,
        max_buffer=2048,
        force_update_cache=True,
        increment_update_cache=None,
        dtype=np.float64,
    ):
        """
        Create a rolling buffer of frame data

        Parameters
        ----------
        frame_src: str of list of strings
            Strings that indicate where to read from files from. This can be a
        list of frame files, a glob, etc.
        channel_name: str
            Name of the channel to read from the frame files
        start_time:
            Time to start reading from.
        max_buffer: {int, 2048}, Optional
            Length of the buffer in seconds
        dtype: {dtype, np.float32}, Optional
            Data type to use for the interal buffer

        """
        self.frame_src = frame_src
        self.channel_name = channel_name
        self.read_pos = start_time
        self.force_update_cache = force_update_cache
        self.increment_update_cache = increment_update_cache
        self.detector = channel_name.split(":")[0]

        self.update_cache()
        self.channel_type, self.raw_sample_rate = self._retrieve_metadata(
            self.cache, self.channel_name
        )

        raw_size = self.raw_sample_rate * max_buffer
        self.raw_buffer = TimeSeries(
            zeros(raw_size, dtype=dtype),
            copy=False,
            epoch=start_time - max_buffer,
            delta_t=1.0 / self.raw_sample_rate,
        )

    def update_cache(self):
        """
        Reset the cache. This can be used to update the cache if the
        result may change due to more files being added to the filesystem,
        for example.
        """
        self.cache = locations_to_cache(self.frame_src, latest=True)

    @staticmethod
    def _retrieve_metadata(cache, channel_name):
        """
        Retrieve basic metadata by reading the first file in the cache

        Parameters
        ----------
        cache: FrameCache
            Cache containing a channel we want to learn about
        channel_name: str
            The name of the channel we want to know the dtype and sample rate of

        Returns
        -------
        channel_type: np.dtype
            The dtype of the channel
        sample_rate: int
            The sample rate of the data within this channel

        """
        path = cache[0].path
        dtype, sample_rate = file_format(path).channel_info(path, channel_name)
        return dtype, int(sample_rate)

    def _read_frame(self, blocksize):
        """
        Try to read the block of data blocksize seconds long

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
            return _read_channel(
                self.cache,
                self.channel_name,
                self.read_pos,
                self.read_pos + int(blocksize),
            )
        except Exception as exc:
            raise RuntimeError(f"Cannot read {self.channel_name} frame data") from exc  # noqa: TRY003

    def null_advance(self, blocksize):
        """
        Advance and insert zeros

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel

        """
        self.raw_buffer.roll(-int(blocksize * self.raw_sample_rate))
        self.read_pos += blocksize
        self.raw_buffer.start_time += blocksize

    def advance(self, blocksize):
        """
        Add blocksize seconds more to the buffer, push blocksize seconds
        from the beginning.

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel

        """
        ts = self._read_frame(blocksize)

        self.raw_buffer.roll(-len(ts))
        self.raw_buffer[-len(ts) :] = ts[:]
        self.read_pos += blocksize
        self.raw_buffer.start_time += blocksize
        return ts

    def update_cache_by_increment(self, blocksize):
        """
        Update the internal cache by starting from the first frame
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

        if not hasattr(self, "dur"):
            fname = glob.glob(self.frame_src[0])[0]  # noqa: PTH207
            fname = Path(fname).stem.split("-")

            self.beg = "-".join([fname[0], fname[1]])
            self.ref = int(fname[2])
            self.dur = int(fname[3])

        fstart = int(
            self.ref + np.floor((start - self.ref) / float(self.dur)) * self.dur
        )
        starts = np.arange(fstart, end, self.dur).astype(int)

        keys = []
        for s in starts:
            pattern = self.increment_update_cache
            if "GPS" in pattern:
                n = int(pattern[int(pattern.index("GPS") + 3)])
                pattern = pattern.replace(f"GPS{n}", str(s)[0:n])

            name = f"{pattern}/{self.beg}-{s}-{self.dur}.gwf"
            # check that file actually exists, else abort now
            if not Path(name).exists():
                raise RuntimeError

            keys.append(name)
        self.cache = locations_to_cache(keys)
        self.channel_type, self.raw_sample_rate = self._retrieve_metadata(
            self.cache, self.channel_name
        )

    def attempt_advance(self, blocksize, timeout=10):
        """
        Attempt to advance the frame buffer. Retry upon failure, except
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
            except RuntimeError:  # noqa: PERF203
                if pycbc.gps_now() > timeout + self.raw_buffer.end_time:
                    # The frame is not there and it should be by now,
                    # so we give up and treat it as zeros
                    DataBuffer.null_advance(self, blocksize)
                    return None
                # I am too early to give up on this frame,
                # so we should try again
                time.sleep(0.1)


class StatusBuffer(DataBuffer):
    """Read state vector or DQ information from a frame file"""

    def __init__(
        self,
        frame_src,
        channel_name,
        start_time,
        max_buffer=2048,
        valid_mask=3,
        force_update_cache=False,
        increment_update_cache=None,
        valid_on_zero=False,
    ):
        """
        Create a rolling buffer of status data from a frame

        Parameters
        ----------
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
        DataBuffer.__init__(
            self,
            frame_src,
            channel_name,
            start_time,
            max_buffer=max_buffer,
            force_update_cache=force_update_cache,
            increment_update_cache=increment_update_cache,
            dtype=np.int32,
        )

        self.valid_mask = valid_mask
        self.valid_on_zero = valid_on_zero

    def check_valid(self, values, flag=None):
        """
        Check if the data contains any non-valid status information

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
            valid = np.bitwise_and(values.numpy(), flag) == flag
        return bool(np.all(valid))

    def is_extent_valid(self, start_time, duration, flag=None):
        """
        Check if the duration contains any non-valid frames

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
        """
        Return the indices of the times lying in the flagged region

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
        indices: np.ndarray
            Array of indices marking the location of triggers within valid
        time.

        """
        sr = self.raw_buffer.sample_rate
        s = int((start_time - self.raw_buffer.start_time - padding) * sr) - 1
        e = s + int((duration + padding) * sr) + 1
        data = self.raw_buffer[s:e]
        stamps = data.sample_times.numpy()

        if self.valid_on_zero:
            invalid = data.numpy() != 0
        else:
            invalid = np.bitwise_and(data.numpy(), self.valid_mask) != self.valid_mask

        starts = stamps[invalid] - padding
        ends = starts + 1.0 / sr + padding * 2.0
        return indices_outside_times(times, starts, ends)

    def advance(self, blocksize):
        """
        Add blocksize seconds more to the buffer, push blocksize seconds
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


class iDQBuffer:  # noqa: N801
    """Read iDQ timeseries from a frame file"""

    def __init__(
        self,
        frame_src,
        idq_channel_name,
        idq_status_channel_name,
        idq_threshold,
        start_time,
        max_buffer=512,
        force_update_cache=False,
        increment_update_cache=None,
    ):
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
        self.idq = DataBuffer(
            frame_src,
            idq_channel_name,
            start_time,
            max_buffer=max_buffer,
            force_update_cache=force_update_cache,
            increment_update_cache=increment_update_cache,
        )
        self.idq_state = DataBuffer(
            frame_src,
            idq_status_channel_name,
            start_time,
            max_buffer=max_buffer,
            force_update_cache=force_update_cache,
            increment_update_cache=increment_update_cache,
        )

    def flag_at_times(self, start_time, duration, times, padding=0):
        """
        Check whether the idq flag was on at given times

        Parameters
        ----------
        start_time: int
            Beginning time to request for
        duration: int
            Number of seconds to check.
        times: array of floats
            Times to check for an active flag
        padding: float
            Amount of time in seconds to flag around samples
            below the iDQ FAP threshold

        Returns
        -------
        flag_state: np.ndarray
            Boolean array of whether flag was on at given times

        """
        # convert start and end times to buffer indices
        sr = self.idq.raw_buffer.sample_rate
        s = int((start_time - self.idq.raw_buffer.start_time - padding) * sr) - 1
        e = s + int((duration + padding) * sr) + 1

        # find samples when iDQ FAP is below threshold and state is valid
        idq_fap = self.idq.raw_buffer[s:e]
        low_fap = idq_fap.numpy() <= self.threshold
        idq_valid = self.idq_state.raw_buffer[s:e]
        idq_valid = idq_valid.numpy().astype(bool)
        valid_low_fap = np.logical_and(idq_valid, low_fap)

        # find times corresponding to the valid low FAP samples
        glitch_idx = np.flatnonzero(valid_low_fap)
        stamps = idq_fap.sample_times.numpy()
        glitch_times = stamps[glitch_idx]

        # construct start and end times of flag segments
        starts = glitch_times - padding
        ends = starts + 1.0 / sr + padding * 2.0

        # check if times were flagged
        idx = indices_within_times(times, starts, ends)
        flagged_bool = np.zeros(len(times), dtype=bool)
        flagged_bool[idx] = True
        return flagged_bool

    def advance(self, blocksize):
        """
        Add blocksize seconds more to the buffer, push blocksize seconds
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
        """
        Advance and insert zeros

        Parameters
        ----------
        blocksize: int
            The number of seconds to advance the buffers

        """
        self.idq.null_advance(blocksize)
        self.idq_state.null_advance(blocksize)
