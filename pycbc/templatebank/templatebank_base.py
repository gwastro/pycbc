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
Base class of templatebanks
"""

from abc import ABCMeta, abstractmethod, abstractproperty
from math import * 
import random


class TemplateBankBase(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self,  n_templates, waveform_length,
                 waveform_frequency_series_t):
        
        # init members
        
        self.__templates= n_templates
        self.__template_index = 0
        self.__template_params= []

        for i in range(self.__templates):
            tmp = [i*1.0, i*0.5] # m1, m2 very prototyping model of parameter space
            self.__template_params.append(tmp)

        print self.__template_params

        self.__waveform_length= waveform_length
        self.__waveform_frequency_series_t = waveform_frequency_series_t
        
        self.__template_params= int(3)
        
        # setup initial data vectors            
        self.__waveform = self.__waveform_frequency_series_t(self.__waveform_length)
        

    #-interface-----------------------------------------------------------------

    @property
    def waveform_length(self):
        return self.__waveform_length

    #---------------------------------------------------------------------------


    def next(self):
        if self.__template_index == self.__templates:
            raise StopIteration
        self.__template_index = self.__template_index + 1
        return self.__template_params[self.__template_index-1]

    # iterate over waveform buffers  
    #def next(self):
    #    if self.__waveform_index == self.__templates:
    #        self.__waveform_index = 0
    #        raise StopIteration
    #    self.__waveform_index = self.__waveform_index + 1
    #    return self.__waveforms[self.__waveform_index-1]

    # define the iterater of TemplateBank. Other access patterns to the data 
    # should be implemented by generators (i.g. reverse())
    def __iter__(self):
        """
        define TemplateBank itself to iterate over it's inherent list of
        waveform series
        """
        return self
        
    def read_parameter_space(self):
        """
        read the parameterspace of the template bank from a file (LAL)
        """
        pass
                
    def perform_generate_waveform(self, template):
        """
        generate waveform on target architecture/device - get's into the 
        ABC scheme later (like in Matchedfilter
        perform_generate_snr ...) 
        """
            
        # call the right approximation model generator with "template"!
        for i in range(self.waveform_length):
            self.__waveform[i] = random.uniform(-1,1)
                
        return self.__waveform
            
         
            
        
        

