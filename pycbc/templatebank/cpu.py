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
TemplateBank Cpu implementation class for the pycbc package
"""

from pycbc.cpu import CpuProcessingObj
from pycbc.templatebank.base import TemplateBankBase

from pycbc.datavector.clayer.cpu import complex_vector_single_cpu_t
from pycbc.datavector.clayer.cpu import real_vector_single_cpu_t

from pycbc.waveformgenerator.cpu import WaveFormGeneratorCpu

import logging

class TemplateBankCpu(TemplateBankBase, CpuProcessingObj):

    def __init__(self, 
                 context, 
                 waveform_length=0, 
                 waveform_delta_x=1,
                 template_table_filename=None):
                 
        self.__logger= logging.getLogger('pycbc.TemplateBankCpu')
        self.__logger.debug("instanciate TemplateBankCpu")   
        
        super(TemplateBankCpu, self).__init__(
                context,  
                waveform_length= waveform_length, 
                waveform_delta_x= waveform_delta_x,
                template_table_filename= template_table_filename,
                filter_waveform_vector_t= complex_vector_single_cpu_t,
                pre_condition_vector_t=real_vector_single_cpu_t, # precodition data 
                                                                 # is currently 
                                                                 # fixed to a cpu vector
                waveform_generator_t= WaveFormGeneratorCpu) 