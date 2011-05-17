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

from matchedfiltercpu import gen_snr_cpu

class MatchedFilterCpu(MatchedFilterBase):

    def __init__(self, length=0):
        print "instanciated MatchedFilterCpu" 
        
        # Instanciate generate-snr-implementation in base class  
        super(MatchedFilterCpu, self).__init__(length, 
              GenSnrImplementationCpu, MaxImplementationCpu)

class  GenSnrImplementationCpu(GenSnrImplementationBase):

    def __init__(self):
        print "instanciated GenSnrImplementationCpu" 
        super(GenSnrImplementationCpu, self).__init__()
    
    def generate_snr(self, stilde, htilde, snr):
        """
        Process matched filtering by generating snr timeseries \rho(t)
        """
        assert repr(stilde).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for stilde"
        assert repr(htilde).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for htilde"
        assert repr(snr).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for snr"

        gen_snr_cpu(stilde, htilde, snr)
        
        return 0

class  MaxImplementationCpu(MaxImplementationBase):

    def __init__(self):
        print "instanciated MaxImplementationCpu" 
        super(MaxImplementationCpu, self).__init__()
    
    def max(self, snr):
        """
        Find the maximum in the generated snr timeseries \rho(t)
        """
        assert repr(snr).find("datavectorcpu") >= 0, "try to call gen_snr_cpu() with wrong type of datavector for snr"
        
        print "finding maximum of snr series - to be implemented in clayer"             
        return 5.5


