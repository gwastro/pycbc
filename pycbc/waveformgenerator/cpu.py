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
Waveform Generator Cpu implementation class for the pycbc package
"""

from pycbc.cpu import CpuProcessingObj
from pycbc.waveformgenerator.base import WaveFormGeneratorBase

import logging

class WaveFormGeneratorCpu(WaveFormGeneratorBase, CpuProcessingObj):

    def __init__(self, 
                 context, 
                 waveform_length=0, 
                 waveform_delta_x=1,
                 approximation_model=None):
                 
        self.__logger= logging.getLogger('pycbc.WaveFormGeneratorCpu')
        
        super(WaveFormGeneratorCpu, self).__init__(
                                context,  
                                waveform_length=waveform_length, 
                                waveform_delta_x=waveform_delta_x,
                                approximation_model=approximation_model) 