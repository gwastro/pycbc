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
from pycbc.waveformgenerator.clayer.cpu import gen_waveform_filter_TaylorF2

import logging

class WaveFormGeneratorCpu(WaveFormGeneratorBase, CpuProcessingObj):

    def __init__(self, 
                 context, 
                 waveform_length=0, 
                 waveform_delta_x=1,
                 approximation_model=None,
                 waveform_filter=None):
                 
        self.__logger= logging.getLogger('pycbc.WaveFormGeneratorCpu')
        self.__logger.debug("instanciate WaveFormGeneratorCpu")
        
        super(WaveFormGeneratorCpu, self).__init__(
                                context,  
                                waveform_length=waveform_length, 
                                waveform_delta_x=waveform_delta_x,
                                approximation_model=approximation_model,
                                waveform_filter=waveform_filter) 
                                
        # mapping clayer functions according to approximation model
        
        self._gen_precon_map={"TaylorF2":GenPreconVecTaylorF2(), 
                              "SpinTaylorT4":None}                                                      
        self.__logger.debug("created a map of aproximation-models to " +
              "generate-precondition-functors {0}".format(self._gen_precon_map))
        
        self.gen_precon_vector= self._gen_precon_map[self._approximation_model]
        self.__logger.debug("mapped self.gen_precon_vector functor to " + 
                            "{0}".format(self.gen_precon_vector))
        
        self._gen_waveform_filter_map={"TaylorF2":GenWaveformFilterTaylorF2(), 
                                       "SpinTaylorT4":None}                                                      
        self.__logger.debug("created a map of aproximation-models to " +
              "generate-waveform_filter-functors {0}".format(self._gen_waveform_filter_map))
        
        self.gen_waveform_filter= self._gen_waveform_filter_map[self._approximation_model]
        self.__logger.debug("mapped self.gen_waveform_filter functor to " + 
                            "{0}".format(self.gen_waveform_filter))
        
        
    
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

    # implementation of ABC's abstractmethod
    def perform_generate_waveform_filter(self, template):

        self.__logger.debug("called perform_generate_waveform_filter")
        
        self.gen_waveform_filter(template, self._waveform_filter)
        
        return self._waveform_filter         
        
class GenPreconVecTaylorF2:
    """
    functor definition for gen_precon_vector_TaylorF2
    """
    
    def __init__(self):
        # status variables etc for/in clayer regrading this function goes here !
        pass
    
    def __call__(self, precon_vec ):
    
        return gen_precon_vector_TaylorF2(precon_vec)
    
            
class GenWaveformFilterTaylorF2:
    """
    functor definition for gen_waveform_filter_TaylorF2
    """
    
    def __init__(self):
        # status variables etc for/in clayer regrading this function goes here !
        pass
    
    def __call__(self, single_inspiral_table_row, waveform_filter):
        
        # ToDo: fetch proper parameters from single inspiral table according
        # to TaylorF2 approximation (clayer function)
        
        gen_waveform_filter_TaylorF2(waveform_filter, 1.2, 1.8)
 
