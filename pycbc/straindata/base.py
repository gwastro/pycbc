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

from abc import ABCMeta, abstractmethod, abstractproperty

from math import *
 
import logging

# Swigged C-layer functions
from pycbc.straindata.clayer.cpu import frame_cpp_read_frames

class StrainDataBase:
    """
    Holds straindata from interferometers in different representations:
    - Initial data
    - Segmented and fourier transformed data
    Transfers data to the appropriate processing architecture 
    after preprocessing (Render data)
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self, t_start, t_end, n_segments, sample_freq,
                 interferometer, 
                 initial_time_series_t= None,
                 time_series_t=         None, 
                 frequency_series_t=    None,
                 fft_segments_impl_t=   None):
        """
        Constructor of the straindata base class
        @type t_start: int
        @param t_start: GPS start time of the strain data 
        @type t_end: int
        @param t_end: GPS start time of the strain data
        @type n_segments: int
        @param n_segments: Number of the overlapping segments to transform
        @type sample_freq: double
        @param sample_freq: Sample frequency in Hz
        @type interferometer: string
        @param interferometer: Interferometer name (H0, H1, L1 etc.)
        @type initial_time_series_t: class
        @param initial_time_series_t: real_vector_double_<arch>_t initial straindata
        @type time_series_t: class
        @param time_series_t: real_vector_single_<arch>_t straindata after being 
                              converted to single precision
        @type frequency_series_t: class
        @param frequency_series_t: complex_vector_single_<arch>_t segmented straindata
                                   transformed into frequency domain
        @type fft_segments_impl_t: class
        @param fft_segments_impl_t: implementation class of segmenting (and 
                                    fft transforming)
        @rtype: none
        @return: none
        """
   
        self.__logger= logging.getLogger('pycbc.StrainDataBase')
        
        # init members
        self.__sample_freq= sample_freq
        self.__length= (t_end - t_start) * self.__sample_freq
        self.__segments= n_segments
        self.__segments_index = 0
        self.__overlap_fact= 0.5
        self.__interferometer = interferometer

        # calculate and check segment length -----------------------------------
        self.__segments_length = \
            self.__length / (1-self.__overlap_fact) / (self.__segments + 1)
            
        seg_len_base = log(self.__segments_length,2)
        assert seg_len_base == int(seg_len_base), "calculated segment length "+\
        "from parameters: {0} not radix 2 ".format(self.__segments_length)
        
        self.__segments_length = int(self.__segments_length)

        # setup datavectors
        # save the type of the time_series for convert_to_single_preci()
        self.__time_series_t = time_series_t

        # setup initial data time series            
        self.__time_series = initial_time_series_t(self.__length, 
                                                   1.0/self.__sample_freq)
        
        # TODO This is currently only prototyping a call to the swigged 
        # C-layer function to read in real straindata from a 3rd party 
        # library like frame-cpp
        self.read_frames("channel_name", t_start, t_end, "cache_filename")
        
        #setup segmented frequency series stilde(f)
        self.__frequency_series= []
        for i in range(self.__segments):
            tmp_series = frequency_series_t(self.__segments_length, 
                                            self.__sample_freq)
            self.__frequency_series.append(tmp_series)

        # instanciate the (fft) segmenting-implementation object
        self.__fft_segments_impl = fft_segments_impl_t(self.__segments_length, 
                    self.__overlap_fact, time_series_t, frequency_series_t)
        assert isinstance(self.__fft_segments_impl, 
                          FftSegmentsImplementationBase), "fft_segments "+\
                          "implementation class is not derivate \n\
                          FftSegmentsImplementationBase"
        
        self.__logger.debug("Instanciated with:")
        self.__logger.debug("self.__sample_freq: {0}".
                            format(self.__sample_freq))
        self.__logger.debug("self.__length: {0}".
                            format(self.__length))
        self.__logger.debug("self.__segments: {0}".
                            format(self.__segments))
        self.__logger.debug("self.__segments_index: {0}".
                             format(self.__segments_index))
        self.__logger.debug("self.__segments_length: {0}".
                             format(self.__segments_length))
        self.__logger.debug("self.__overlap_fact: {0}".
                            format(self.__overlap_fact))
        self.__logger.debug("self.__interferometer: {0}".
                            format(self.__interferometer))
        self.__logger.debug("self.__time_series_t: {0}".
                            format(self.__time_series_t))

    # Define the iterater of StrainData. Other access patterns to the data 
    # should be implemented by generators (i.g. reverse())
    def __iter__(self):
        """
        define straindata itself to iterate over it's inherent list of segments
        """
        return self

    def next(self):
        """
        define straindata's next methode to iterate over it's inherent list 
        of segments
        """

        if self.__segments_index == self.__segments:
            self.__segments_index = 0
            raise StopIteration
        self.__segments_index = self.__segments_index + 1
        return self.__frequency_series[self.__segments_index-1]

    # Multiplication with itself will be used to apply the overwhitening filter. 
    # (complex *= real operation)
    def __rmul__(self):
        """
        Overload mul for data objects, to allow multiplcation by a psd series
        """
        # TODO
        pass

    #-properties----------------------------------------------------------------

    @property
    def time_series(self):
        return self.__time_series

    #@time_series.setter
    #def time_series(self, value):
    #    self.__time_series = value

    @property
    def frequency_series(self):
        return self.__frequency_series

    #@frequency_series.setter
    #def frequency_series(self, value):
    #    self.__frequency_series = value

    @property
    def segments_length(self):
        return self.__segments_length

    @property
    def segments_delta_x(self):
        return self.__frequency_series[0].get_delta_x()

    #---------------------------------------------------------------------------
    
    @abstractmethod
    def render(self):
        """
        Force every derivate of StrainDataBase to implement a render() methode
        """
        pass  

    def perform_fft_segments(self):
        """
        Call the fft_segments() methode has been specified by the constructor
        of StrainDataBase
        
        """
        return self.__fft_segments_impl.fft_segments(self.__time_series, 
               self.__frequency_series)
    
    def convert_to_single_preci(self):
        """
        Convert the initial double precision timeseries to single precision
        
        """
        tmp_series= self.__time_series_t(self.__length, 1.0/self.__sample_freq)
        for i in range(self.__length):
            tmp_series[i] = self.__time_series[i]
        self.__time_series = tmp_series
        
                                                  
    def read_frames(self, channel_name, gps_start_time, gps_end_time,cache_url):
        """
        @type  channel_name: string
        @param channel_name: input gw strain channel name 
        @type gps_start_time: unsigned long
        @param gps_start_time: gps start_time of data to be read in
        @type gps_end_time:  unsigned long
        @param gps_end_time: gps end_time of data to be read in
        @type  cache_url: string
        @param cache_url: URL of a lal frame cache file
        """

        # Redirection to swigged C-layer function:
 
        frame_cpp_read_frames(self.time_series, channel_name, 
                              gps_start_time, gps_end_time, 
                              cache_url)

 
class FftSegmentsImplementationBase:
    """
    Base class of fft segmenting 
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        pass
        
    @abstractmethod
    def fft_segments(self, input_buf, output_buf):
        """
        Force every derivate of FftSegmentsImplementationBase to implement 
        a fft_segments() methode
        """
        pass
                       
