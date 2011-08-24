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

from matchedfilter_base import MatchedFilterBase
from matchedfilter_base import GenSnrImplementationBase
from matchedfilter_base import MaxImplementationBase

# import processing functions from the clayer 
from matchedfiltercpu import gen_snr_cpu

# import member datavectors
from pycbc.datavector.datavectorcpu import complex_vector_single_t

import logging

class MatchedFilterCpu(MatchedFilterBase):

    def __init__(self, length=0, delta_x=1):

        self.__logger= logging.getLogger('pycbc.MatchedFilterCpu')
        self.__logger.debug("instanciated MatchedFilterCpu")

        self.__length= length
        self.__delta_x= delta_x
        
        super(MatchedFilterCpu, self).__init__(self.__length, self.__delta_x,
              GenSnrImplementationCpu, MaxImplementationCpu, 
              complex_vector_single_t, complex_vector_single_t, self)


class  GenSnrImplementationCpu(GenSnrImplementationBase):

    def __init__(self, owner_mfilt):

        self.__logger= logging.getLogger('pycbc.GenSnrImplementationCpu')
        self.__logger.debug("instanciated GenSnrImplementationCpu")

        super(GenSnrImplementationCpu, self).__init__(owner_mfilt)
    
    def generate_snr(self, context, snr, stilde, htilde):
        """
        Process matched filtering by generating snr timeseries \rho(t)
        """
        assert repr(stilde).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for stilde"
        assert repr(htilde).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for htilde"
        assert repr(snr).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for snr"

        gen_snr_cpu(context, snr, stilde, htilde, self._owner_mfilt._q, 
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
        assert repr(snr).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for snr"
        
        self.__logger.debug("finding maximum of snr series - to be implemented in clayer")            

        return 5.5


