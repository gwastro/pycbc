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
MatchedFilter Cpu implementation class for the pycbc package
"""

from pycbc.matchedfilter.base import MatchedFilterBase
from pycbc.matchedfilter.base import GenSnrImplementationBase
from pycbc.matchedfilter.base import MaxImplementationBase

# import processing functions from the clayer 
from pycbc.matchedfilter.clayer.cpu import gen_snr_cpu

# import member datavectors
from pycbc.datavector.clayer.cpu import complex_vector_single_cpu_t
from pycbc.datavector.clayer.cpu import real_vector_single_cpu_t

# for testing data_in
#from pycbc.datavector.datavectoropencl import real_vector_single_opencl_t as alien_datavector_t

from pycbc.cpu import CpuProcessingObj

# import FFT members

# Josh to fix the fftwf fftw issue when the setup.py system allows to link to fftw
#from pycbc.fft.fftcpu import fft_real_single_plan
#from pycbc.fft.fftcpu import execute_real_single_reverse_fft_cpu

import logging

class MatchedFilterCpu(MatchedFilterBase, CpuProcessingObj):

    def __init__(self, context, length=0, delta_x=1):

        self.__logger= logging.getLogger('pycbc.MatchedFilterCpu')
        self.__logger.debug("instanciated MatchedFilterCpu")

        self.__length= length
        self.__delta_x= delta_x
        
        super(MatchedFilterCpu, self).__init__(context, self.__length, 
                               self.__delta_x,
                               GenSnrImplementationCpu, MaxImplementationCpu, 
                               snr_vector_t=    real_vector_single_cpu_t, 
                               qtilde_vector_t= complex_vector_single_cpu_t, 
                               q_vector_t =     complex_vector_single_cpu_t, 
                               derived_mfilt =  self)


class  GenSnrImplementationCpu(GenSnrImplementationBase):

    def __init__(self, owner_mfilt):

        self.__logger= logging.getLogger('pycbc.GenSnrImplementationCpu')
        self.__logger.debug("instanciated GenSnrImplementationCpu")

        super(GenSnrImplementationCpu, self).__init__(owner_mfilt)
    
    def generate_snr(self, stilde, htilde):
        """
        Process matched filtering by generating snr timeseries \rho(t)
        """
                
        # check and possibly transfer input datavectors
        stilde= self._owner_mfilt.data_in(stilde)
        htilde= self._owner_mfilt.data_in(htilde)
            
        gen_snr_cpu(self._owner_mfilt._devicecontext, self._owner_mfilt._snr, 
                    stilde, htilde, 
                    self._owner_mfilt._q, 
                    self._owner_mfilt._qtilde, 1.0, 1.0) # just for prototyping
                                                         # f_min and sigma_sq
                                                         # are magic numbers

class  MaxImplementationCpu(MaxImplementationBase):

    def __init__(self, owner_mfilt):

        self.__logger= logging.getLogger('pycbc.MaxImplementationCpu')
        self.__logger.debug("instanciated MaxImplementationCpu") 

        super(MaxImplementationCpu, self).__init__(owner_mfilt)
    
    def max(self, snr):
        """
        Find the maximum in the generated snr timeseries \rho(t)
        """
        assert repr(snr).find("pycbc.datavector.clayer.cpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for snr"
        
        self.__logger.debug("finding maximum of snr series - to be implemented in clayer")            

        return 5.5


