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

    def __init__(self, length, gen_snr_impl):
        print "MatchedFilterBase.__init__ called" 
        self.length = length
        self.gen_snr_impl = gen_snr_impl()

    def perform_generate_snr(self, stilde, htilde):
        """
        Fill this objects internal memory with \rho(t)
        """
        self.gen_snr_impl.generate_snr(stilde, htilde)
        

class GenSnrImplementationBase:
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        print "GenSnrImplementationBase.__init__ called" 

    @abstractmethod
    def generate_snr(self, stilde, htilde):
        pass
        
        
        