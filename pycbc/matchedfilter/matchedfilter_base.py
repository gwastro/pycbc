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
Abstract Base Class (ABC) of matched filters in the pycbc package based on:
<http://www.python.org/dev/peps/pep-3119/>
"""

from abc import ABCMeta, abstractmethod, abstractproperty

import logging
       
class MatchedFilterBase:

    __metaclass__ = ABCMeta

    def __init__(self, length, gen_snr_impl, max_impl, derived_mfilt):
        self.__logger= logging.getLogger('pycbc.MatchedFilterBase')
        self.__logger.debug("instanciated MatchedFilterBase")
        self.__length = length
        self.__gen_snr_impl = gen_snr_impl(derived_mfilt)
        if not isinstance(self.__gen_snr_impl, GenSnrImplementationBase):
            self.__logger.debug("MatchedFilterBase.__init__: gen_snr_impl is not a derivate of GenSnrImplementationBase")
            exit(0)
        self.__max_impl = max_impl()
        if not isinstance(self.__max_impl, MaxImplementationBase):
            self.__logger.debug("MatchedFilterBase.__init__: max_impl is not a derivate of MaxImplementationBase")
            exit(0)

    #-properties----------------------------------------------------------------
 
    @property
    def length(self):
        return self.__length

    #---------------------------------------------------------------------------
            
        
    def perform_generate_snr(self, context, snr, stilde, htilde):
        """
        calls the generate_snr methode of the derived implementation object
        @type  context: Device Context
        @param context: Input: Device Context
        @type  snr:     DataVectorBase
        @param snr:     Output: Signal to noise ratio series
        @type  stilde: DataVectorBase
        @param stilde: Input: Straindata frequency domain
        @type  htilde: DataVectorBase
        @param htilde: Input: Template waveform frequency domain
        """
        self.__gen_snr_impl.generate_snr(context, snr, stilde, htilde)


    def perform_max(self, snr):
        """
        calls the max methode of the derived implementation object
        @rtype:  snr:  DataVectorBase
        @return: snr:  Signal to noise ratio series
        @rtype:  float
        @return: Maximum of snr series
        """
        return self.__max_impl.max(snr)
        

class GenSnrImplementationBase:
    
    __metaclass__ = ABCMeta
    
    def __init__(self, owner_mfilt):
        self.__logger= logging.getLogger('pycbc.GenSnrImplementationBase')
        self.__logger.debug("instanciated GenSnrImplementationBase")
        self._owner_mfilt = owner_mfilt
        
    @abstractmethod
    def generate_snr(self, stilde, htilde ,snr):
        pass
        
class MaxImplementationBase:
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        self.__logger= logging.getLogger('pycbc.MaxImplementationBase')
        self.__logger.debug("instanciated MaxImplementationBase") 
        
    @abstractmethod
    def max(self, snr):
        pass
        
        