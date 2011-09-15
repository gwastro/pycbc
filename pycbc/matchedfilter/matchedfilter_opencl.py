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
MatchedFilter OpenCl implementation class for the pycbc package
"""

# import framework related parent
from pycbc.pycbc import OpenClProcessingObj

# import algorithm related parents
from matchedfilter_base import MatchedFilterBase
from matchedfilter_base import GenSnrImplementationBase
from matchedfilter_base import MaxImplementationBase

# import processing functions from the clayer 
from matchedfilteropencl import gen_snr_opencl
from matchedfilteropencl import matched_filter_opencl_t


# import member datavectors
from pycbc.datavector.datavectoropencl import complex_vector_single_opencl_t
from pycbc.datavector.datavectoropencl import real_vector_single_opencl_t


import logging

class MatchedFilterOpenCl(MatchedFilterBase, OpenClProcessingObj):

    def __init__(self, context, length=0, delta_x=1):
        self.__logger= logging.getLogger('pycbc.MatchedFilterOpenCl')
        self.__logger.debug("instanciated MatchedFilterOpenCl") 

        self.__length= length
        self.__delta_x= delta_x
        
        self.__mf_opencl_clayer_struct= matched_filter_opencl_t()
        print matched_filter_opencl_t
        
        super(MatchedFilterOpenCl, self).__init__(context, self.__length, 
                self.__delta_x, 
                GenSnrImplementationOpenCl, MaxImplementationOpenCl,
                snr_vector_t=    real_vector_single_opencl_t, 
                qtilde_vector_t= complex_vector_single_opencl_t, 
                q_vector_t =     complex_vector_single_opencl_t, 
                derived_mfilt =  self)


class  GenSnrImplementationOpenCl(GenSnrImplementationBase):

    def __init__(self, owner_mfilt):
    
        self.__logger= logging.getLogger('pycbc.GenSnrImplementationOpenCl')
        self.__logger.debug("instanciated GenSnrImplementationOpenCl")
    
        super(GenSnrImplementationOpenCl, self).__init__(owner_mfilt)
    
    def generate_snr(self, stilde, htilde):
        """
        Process matched filtering by generating snr timeseries \rho(t)
        """

        # check and possibly transfer input datavectors
        stilde= self._owner_mfilt.data_in(stilde)
        htilde= self._owner_mfilt.data_in(htilde)
        
        # this finally calls the clayer function:
        gen_snr_opencl(self._owner_mfilt._devicecontext, self._owner_mfilt._snr,
                       stilde, htilde)
        
class  MaxImplementationOpenCl(MaxImplementationBase):

    def __init__(self, owner_mfilt):
    
        self.__logger= logging.getLogger('pycbc.MaxImplementationOpenCl')
        self.__logger.debug("instanciated MaxImplementationOpenCl") 
    
        super(MaxImplementationOpenCl, self).__init__(owner_mfilt)
    
    def max(self, snr):
        """
        Find the maximum in the generated snr timeseries \rho(t)
        """

        assert repr(snr).find("datavectoropencl") >= 0, "try to call gen_snr_opencl() with wrong type of datavector for snr"
        
        self.__logger.debug("finding maximum of snr series - to be implemented in clayer")            
        return 5.5
        

