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
        
class MatchedFilterBase:

    __metaclass__ = ABCMeta

    def __init__(self, length, gen_snr_impl, max_impl):
        print "instanciated MatchedFilterBase" 
        self.length = length
        self.gen_snr_impl = gen_snr_impl()
        if not isinstance(self.gen_snr_impl, GenSnrImplementationBase):
            print "MatchedFilterBase.__init__: gen_snr_impl is not a derivate of GenSnrImplementationBase "
            exit(0)
        self.max_impl = max_impl()
        if not isinstance(self.max_impl, MaxImplementationBase):
            print "MatchedFilterBase.__init__: max_impl is not a derivate of MaxImplementationBase "
            exit(0)
            
        
    def perform_generate_snr(self, stilde, htilde, snr):
        """
        calls the generate_snr methode of the derived implementation object
        @type  stilde: DataVectorBase
        @param stilde: Straindata frequency domain
        @type  htilde: DataVectorBase
        @param htilde: Template waveform frequency domain
        @rtype:  snr:  DataVectorBase
        @return: snr:  Signal to noise ratio series
        """
        return self.gen_snr_impl.generate_snr(stilde, htilde, snr)

    def perform_max(self, snr):
        """
        calls the max methode of the derived implementation object
        @rtype:  snr:  DataVectorBase
        @return: snr:  Signal to noise ratio series
        @rtype:  float
        @return: Maximum of snr series
        """
        return self.max_impl.max(snr)
        

class GenSnrImplementationBase:
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        print "instanciated GenSnrImplementationBase" 
        
    @abstractmethod
    def generate_snr(self, stilde, htilde ,snr):
        pass
        
class MaxImplementationBase:
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        print "instanciated MaxImplementationBase" 
        
    @abstractmethod
    def max(self, snr):
        pass
        
        