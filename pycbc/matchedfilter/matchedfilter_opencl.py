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

from matchedfilter_base import MatchedFilterBase
from matchedfilter_base import GenSnrImplementationBase
from matchedfilter_base import MaxImplementationBase

from matchedfilteropencl import gen_snr_opencl

import logging

class MatchedFilterOpenCl(MatchedFilterBase):

    def __init__(self, length=0):
        self.__logger= logging.getLogger('pycbc.MatchedFilterOpenCl')
        self.__logger.debug("instanciated MatchedFilterOpenCl") 
        
        # Instanciate generate-snr-implementation in base class  
        super(MatchedFilterOpenCl, self).__init__(length, 
              GenSnrImplementationOpenCl, MaxImplementationOpenCl)

class  GenSnrImplementationOpenCl(GenSnrImplementationBase):

    def __init__(self):
        self.__logger= logging.getLogger('pycbc.GenSnrImplementationOpenCl')
        self.__logger.debug("instanciated GenSnrImplementationOpenCl")
        super(GenSnrImplementationOpenCl, self).__init__()
    
    def generate_snr(self, stilde, htilde, snr):
        """
        Process matched filtering by generating snr timeseries \rho(t)
        """
        assert repr(stilde).find("datavectoropencl") >= 0, "try to call gen_snr_opencl() with wrong type of datavector for stilde"
        assert repr(htilde).find("datavectoropencl") >= 0, "try to call gen_snr_opencl() with wrong type of datavector for htilde"
        assert repr(snr).find("datavectoropencl") >= 0, "try to call gen_snr_opencl() with wrong type of datavector for snr"
        
        gen_snr_opencl(stilde, htilde, snr)
        
        return 0
        
class  MaxImplementationOpenCl(MaxImplementationBase):

    def __init__(self):
        self.__logger= logging.getLogger('pycbc.MaxImplementationOpenCl')
        self.__logger.debug("instanciated MaxImplementationOpenCl") 
        super(MaxImplementationOpenCl, self).__init__()
    
    def max(self, snr):
        """
        Find the maximum in the generated snr timeseries \rho(t)
        """

        assert repr(snr).find("datavectoropencl") >= 0, "try to call gen_snr_opencl() with wrong type of datavector for snr"
        
        self.__logger.debug("finding maximum of snr series - to be implemented in clayer")            
        return 5.5
        

