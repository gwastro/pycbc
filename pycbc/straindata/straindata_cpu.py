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
Cpu version of strain data class
"""

from straindata_base import StrainDataBase
from straindata_base import FftSegmentsImplementationBase

from pycbc.datavector.datavectorcpu import real_vector_double_t as InitialTimeSeriesDoublePreci
from pycbc.datavector.datavectorcpu import real_vector_single_t as TimeSeriesSinglePreci
from pycbc.datavector.datavectorcpu import complex_vector_single_t as FrequencySeries

from straindatacpu import fftw_generate_plan
from straindatacpu import fftw_transform_segments

import logging

class StrainDataCpu(StrainDataBase):

    def __init__(self, t_start, t_end, n_segments, sample_freq, interferometer):
        
        super(StrainDataCpu, self).__init__(t_start, t_end, n_segments, sample_freq,
                                            interferometer,
                                            InitialTimeSeriesDoublePreci,
                                            TimeSeriesSinglePreci,
                                            FrequencySeries,
                                            FftSegmentsImplementationFftw)

    def render(self):
        pass         # nothing to render in cpu version (data remaines in place)
                                            
                                                                                                                                                                          
class  FftSegmentsImplementationFftw(FftSegmentsImplementationBase):

    def __init__(self, length, overlap_fact, input_buf_t, output_buffer_t):
        print "instanciated FftSegmentsImplementationCpu" 

        assert repr(input_buf_t).find("datavectorcpu") >= 0, "try to instanciate FftSegmentsImplementationFftw with wrong type of datavector for input_buf"
        assert repr(output_buffer_t).find("datavectorcpu") >= 0, "try to instanciate FftSegmentsImplementationFftw with wrong type of datavector for output_buffers_t"

        super(FftSegmentsImplementationFftw, self).__init__()
    
        self.__logger= logging.getLogger('pycbc.FftSegmentsImplementationFftw')
    
        self.__length = length
        self.__overlap_fact = overlap_fact
        self.__input_buf_t = input_buf_t
        self.__output_buffer_t =  output_buffer_t

        # create plan
        in_tmp  = self.__input_buf_t(self.__length)
        out_tmp = self.__output_buffer_t(self.__length)
        self.__fft_forward_plan= fftw_generate_plan(self.__length, in_tmp, out_tmp, "FFTW_FORWARD", "FFTW_ESTIMATE")
    
    def fft_segments(self, input_buf, output_buf):
        """
        Process ffts of strain data segments
        """
        
        self.__logger.debug("self.__fft_forward_plan: {0}".format(self.__fft_forward_plan))

        input_buf_offset = 0
        for output_buffer_segment in output_buf:
            self.__logger.debug("input_buf_offset: {0}".format(input_buf_offset))
            fftw_transform_segments(self.__fft_forward_plan, input_buf, input_buf_offset, output_buffer_segment)
            input_buf_offset = int( input_buf_offset + self.__length * (1 - self.__overlap_fact) )
