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
Base class of template bank
"""

from abc import ABCMeta, abstractmethod, abstractproperty
from math import * 
import random
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import ligolw_add

import logging

class TemplateBank(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self, context, waveform_length, waveform_delta_x,
                 waveform_frequency_series_t):
        self.__logger = logging.getLogger('pycbc.TemplateBank')

        # init members
        self._context = context
        self.__templates_num = 0
        self.__template_index = 0
        self.__template_params = []
        self.__filename = ''
        self.templates = None

        self.__waveform_length = waveform_length
        self.__waveform_delta_x = waveform_delta_x
        
        self.__waveform_frequency_series_t = waveform_frequency_series_t
        
        # setup initial data vectors            
        self.__waveform = self.__waveform_frequency_series_t(self._context,
                                                        self.__waveform_length,
                                                        self.__waveform_delta_x)
    

    #-interface-----------------------------------------------------------------

    @property
    def waveform_length(self):
        return self.__waveform_length

    #---------------------------------------------------------------------------

    # iterate over templates in the parameter space
    # define the iterater of TemplateBank. Other access patterns to the data 
    # should be implemented by generators (e.g., reverse())
    def __iter__(self):
        """
        define the template bank object itself to iterate over it's inherent
        parameter space
        """
        return self

    def next(self):
        if self.__template_index == self.__templates:
            raise StopIteration
        self.__template_index = self.__template_index + 1
        return self.__template_params[self.__template_index-1]

    def read_from_file(self, filename, verbose=False):
        # load the xml file containing the templates and extract the sngl_inspiral table
        self.__logger.debug("loading templatebank: %s"%(filename))
        self.__filename = filename
        xmldoc = ligolw_add.ligolw_add(ligolw.Document(), [filename], verbose=verbose)
        self.templates = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
        self.__logger.debug("obtained %i templates"%(len(self.templates)))
        self.__templates_num = len(self.templates)
        self.__template_index = 0

        self.__logger.debug("extracting templatebank mass parameters")
        m1s = self.templates.get_column('mass1')
        m2s = self.templates.get_column('mass2')
        self.__template_params = zip(m1s,m2s)

    def generate_filter(self, template):
        """
        generate waveform on target architecture/device - get's into the 
        ABC scheme later (like in Matchedfilter
        perform_generate_snr ...) 
        """

        self.__waveform_generator( template)

        return self.__waveform
            
         
            
        
        

