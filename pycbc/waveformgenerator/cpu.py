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

from pycbc.waveformgenerator.clayer.cpu import gen_precon_vector_TaylorF2

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
                                
        
        self._gen_precon_map={"TaylorF2":GenPreconVecTaylorF2(), 
                              "SpinTaylorT4":None}                                                      
        self.__logger.debug("created a map of aproximation-models to " +
              "generate-precondition-functors {0}".format(self._gen_precon_map))
        
        self.gen_precon_vector= self._gen_precon_map[self._approximation_model]
        self.__logger.debug("mapped self.gen_precon_vector functor to " + 
                            "{0}".format(self.gen_precon_vector))
        
        print self.gen_precon_vector
    
    # implementation of ABC's abstractmethod
    def perform_generate_precondition(self, pre_condition_vector_t):

        self.__logger.debug("called perform_generate_precondition")
        
        # pre instanciate precon vector
        precon_vec= pre_condition_vector_t(self._devicecontext,
                                           self.waveform_length,
                                           self.waveform_delta_x,)
        
        # depending on the approx model and thus on the implementation of
        # self.clayer generate_precondition() would return True or False
        # if false return none if true return precon vec! 
        
        if self.gen_precon_vector(precon_vec):
            return precon_vec
        else:
            return None        
        
        
class GenPreconVecTaylorF2:
    """
    functor definition of gen_precon_vector_Tf2_from_row
    """
    
    def __init__(self):
    
        pass
        # status variables etc for/in clayer regrading this function goes here !
    
    def __call__(self, precon_vec ):
     
        #ToDo have to define clayer type for : sngl_insp_tab_row, precon_vec
        return gen_precon_vector_TaylorF2(precon_vec)
        
#ToDo: class GenPreconVecSpinTaylorT4 and so on ...   
