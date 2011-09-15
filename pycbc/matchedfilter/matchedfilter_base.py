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

    def __init__(self, context, length, delta_x, gen_snr_impl, max_impl, 
                 snr_vector_t =    None, 
                 qtilde_vector_t = None, 
                 q_vector_t =      None, 
                 derived_mfilt =   None):

        # actually calling init of <arch>ProcessingObj (that's the way super() does)
        super(MatchedFilterBase, self).__init__(context)

        self.__logger= logging.getLogger('pycbc.MatchedFilterBase')
        self.__logger.debug("instanciated MatchedFilterBase")

        self.__context = context
        self.__length = length
        self.__delta_x = delta_x

        # instanciate member datavectors:
        self._snr= snr_vector_t(self.__context, self.__length, self.__delta_x)
        self._qtilde= qtilde_vector_t(self.__length, self.__delta_x)
        self._q= q_vector_t(self.__length, self.__delta_x)

        # instanciate implementation class members of the processing functions
        self.__gen_snr_impl = gen_snr_impl(derived_mfilt)
        if not isinstance(self.__gen_snr_impl, GenSnrImplementationBase):
            self.__logger.debug("MatchedFilterBase.__init__: gen_snr_impl is not a derivate of GenSnrImplementationBase")
            exit(0)

        self.__max_impl = max_impl(derived_mfilt)
        if not isinstance(self.__max_impl, MaxImplementationBase):
            self.__logger.debug("MatchedFilterBase.__init__: max_impl is not a derivate of MaxImplementationBase")
            exit(0)

    #-properties----------------------------------------------------------------
 
    @property
    def length(self):
        return self.__length

    #---------------------------------------------------------------------------
            
        
    def perform_generate_snr(self, stilde, htilde):
        """
        calls the generate_snr methode of the derived implementation object
        @type  context: Device Context
        @param context: Input: Device Context
        @type  stilde: DataVectorBase
        @param stilde: Input: Straindata frequency domain
        @type  htilde: DataVectorBase
        @param htilde: Input: Template waveform frequency domain
        @rtype  snr:   DataVectorBase
        @return snr:   Output: Signal to noise ratio series
        """
        self.__gen_snr_impl.generate_snr(stilde, htilde)
        
        return self._snr


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
    
    def __init__(self, owner_mfilt):

        self.__logger= logging.getLogger('pycbc.MaxImplementationBase')
        self.__logger.debug("instanciated MaxImplementationBase") 
        self._owner_mfilt = owner_mfilt
        
    @abstractmethod
    def max(self, snr):
        pass
        
        
