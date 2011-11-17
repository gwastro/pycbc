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
Base class of waveform generator
"""


import logging

from abc import ABCMeta, abstractmethod, abstractproperty

class WaveFormGeneratorBase:
    """
    docstring for WaveFormGeneratorBase
    """
    __metaclass__ = ABCMeta

    
    # constructor --------------------------------------------------------------
    def __init__(self, 
                 context, 
                 waveform_length=0, 
                 waveform_delta_x=1,
                 approximation_model=None):
                 
        self.__logger = logging.getLogger('pycbc.WaveFormGeneratorBase')

        # call constructor of <arch>ProcessingObj (2nd parent of WaveFormGeneratorBase 
        # derivative TemplatBank<arch>
        super(WaveFormGeneratorBase, self).__init__(context)
        
        self.waveform_length = waveform_length
        self.waveform_delta_x = waveform_delta_x

 
    def perform_generate_precondition(self, pre_condition_vector_t):
        
        # pre instanciate precon vector
        precon_vec= pre_condition_vector_t(self._devicecontext,
                                           self.waveform_length,
                                           self.waveform_delta_x,)
        
        # CHECK can a methode be handled as a member?
        # call member self.clayer generate_precondition( precon_vec , ...)

        # depending on the approx model and thus on the implementation of
        # self.clayer generate_precondition() would return True or False
        # if false return none if true return precon vec! 
                
        return precon_vec
    
    
    def perform_generate_filterwaveform(self, filter_waveform):
        
        # will be called from the hot loop and reuse the filter_waveform
        # vector for each loop run
                
        pass       
        
        