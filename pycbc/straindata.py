# Copyright (C) 2011 Karsten Wiesner
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
Base class of strain data
"""

import numpy as np

class TimeSeriesGeneric(object):
    def __init__(self, data, length, fs):
        self.data = None
        self.length = None
        self.fs = None
        
    def display(self):
        print self.data, self.data.dtype, self.length
    
    def to_float(self):
        self.data = self.data.astype(np.float32)

class StrainDataGeneric(object):
    def __init__(self):
        print 'instanciated StrainData'
        self.time_series = TimeSeriesGeneric(None, 1, 1)
        self.__frequency_series = None
        self.__global_time_intervall = None

    def __iter__(self):
        """
        provide functionality to iterate over the array of frequency elements
        """
        pass

    def __rmul__(self):
        """
        overload multiply for data objects, e.g. to allow multiplcation by a psd
        """
        pass

    @property
    def time_series(self):
        return self.__time_series
        
    @time_series.setter
    def time_series(self, value):
        self.__time_series = value
    
    def read_frames(self, channel_name, gps_start_time, gps_end_time, cache_url):
        """
        @type  channel_name: string
        @param channel_name: input gravitational wave strain channel name 
        @type gps_start_time: int
        @param gps_start_time: gps start_time of data to be read in
        @type gps_end_time:  int
        @param gps_end_time: gps end_time of data to be read in
        @type  cache_url: string
        @param cache_url: URL of a lal frame cache file
        
        This method fills self.__time_series_data with the data read in from the
        frame. It is responsible for allocating memory in the C layer 
        for the input data in a real_vector_t.
        """
        
        self.__global_time_intervall = gps_end_time - gps_start_time
        
        print 'reading frame data of {0} from {1} to {2} in total {3} s'.format(
               channel_name, gps_start_time, gps_end_time, \
               self.__global_time_intervall)
        
        self.time_series.fs = 8192
        self.time_series.length = self.__global_time_intervall * self.time_series.fs
        self.time_series.data = np.random.rand(self.time_series.length)
        #self.time_series.data = np.ones(self.time_series.length)
                        

    def fft_segments(self, segment_length, segment_overlap):
        """
        split the time_series data into segments and transform each segment into
        the frequency domain using fftw on the cpu.
        """
        pass

    def frequency_series_array_to_cuda_global_memory(self):
        """
        move the frequency series data from cpu memory to the gpu
        """
        pass

