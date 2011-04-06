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

from pycbc.datavector.datavector_base import DataVectorBase
from matchedfilter_base import *

class MatchedFilterOpenCl(MatchedFilterBase):

    def __init__(self, length=0):
        print "MatchedFilterOpenCl.__init__ called" 
        self.gen_snr_impl_open_cl = GenSnrImplementationOpenCl()
        super(MatchedFilterOpenCl, self).__init__(length, GenSnrImplementationOpenCl )
        #super(MatchedFilterOpenCl, self).__init__(length, self.gen_snr_impl_open_cl )
        

class  GenSnrImplementationOpenCl(GenSnrImplementationBase):

    def __init__(self):
        print "GenSnrImplementationOpenCl.__init__ called" 
        super(GenSnrImplementationOpenCl, self).__init__()
    
    def generate_snr(self, stilde, htilde):
        """
        Process matched filtering by generating snr timeseries \rho(t)
        """
        isinstance(stilde, DataVectorBase)
        isinstance(htilde, DataVectorBase)
        print " Generate Snr OpenCl Implementation"


