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

from pycbc.waveformgenerator.clayer.cpu import gen_precon_vector_Tf2_from_row

import logging

class WaveFormGeneratorCpu(WaveFormGeneratorBase, CpuProcessingObj):

    def __init__(self, 
                 context, 
                 waveform_length=0, 
                 waveform_delta_x=1,
                 approximation_model=None):
                 
        self.__logger= logging.getLogger('pycbc.WaveFormGeneratorCpu')
        self.__logger.debug("instanciate WaveFormGeneratorCpu")
        
        super(WaveFormGeneratorCpu, self).__init__(
                                context,  
                                waveform_length=waveform_length, 
                                waveform_delta_x=waveform_delta_x,
                                approximation_model=approximation_model) 
        
        # fetch right approx clayer function via functor
        self.gen_precon_vector= GenPreconVec()                                               
    
    
    # implementation of ABC's abstractmethod
    def perform_generate_precondition(self, sngl_insp_tab_row, pre_condition_vector_t):
        #super(WaveFormGeneratorCpu, self).perform_generate_precondition(pre_condition_vector_t)

        self.__logger.debug("called perform_generate_precondition")
        
        # pre instanciate precon vector
        precon_vec= pre_condition_vector_t(self._devicecontext,
                                           self.waveform_length,
                                           self.waveform_delta_x,)
        
        # depending on the approx model and thus on the implementation of
        # self.clayer generate_precondition() would return True or False
        # if false return none if true return precon vec! 
        
        if self.gen_precon_vector(sngl_insp_tab_row, precon_vec):
            return precon_vec
        else:
            return None        
        
        
class GenPreconVec:
    """
    functor definition of GenPreconVec
    """
    
    def __init__(self):
    
        # we know self._approximation_model here
        # if then else approx func= clayerfunction???
        pass
        # self._statevars in clayer regrading this function goes here !
    
    def __call__(self, sngl_insp_tab_row, precon_vec ):
     
            #ToDo have to define clayer typefor : sngl_insp_tab_row, precon_vec
        return gen_precon_vector_Tf2_from_row(precon_vec)
    
