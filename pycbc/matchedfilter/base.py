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

"""@package docstring
Base class of matched filter
"""
## @package matchedfilter.base
#  Base class of matched filter


from abc import ABCMeta, abstractmethod, abstractproperty

import logging
       
class MatchedFilterBase:

    __metaclass__ = ABCMeta

    def __init__(self, context, 
                 length, delta_x, 
                 snr_vector_t =    None, 
                 qtilde_vector_t = None, 
                 q_vector_t =      None):

        # calling constructor of <arch>ProcessingObj
        super(MatchedFilterBase, self).__init__(context)

        self.__logger= logging.getLogger('pycbc.MatchedFilterBase')
        self.__logger.debug("instanciated MatchedFilterBase")

        self._context = context
        self._length =  length
        self._delta_x = delta_x

        # instanciate member datavectors:
        self._snr= snr_vector_t(self._context, self._length, self._delta_x)
        self._qtilde= qtilde_vector_t(self._context, self._length, self._delta_x)
        self._q= q_vector_t(self._context, self._length, self._delta_x)


    @abstractmethod    
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
        pass        
       
    @abstractmethod
    def perform_find_max(self, snr):
        """
        calls the max methode of the derived implementation object
        @rtype:  snr:  DataVectorBase
        @return: snr:  Signal to noise ratio series
        @rtype:  float
        @return: Maximum of snr series
        """
        pass
        
