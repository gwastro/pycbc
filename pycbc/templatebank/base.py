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

from math import * 
import random

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import ligolw_add

import logging

from abc import ABCMeta, abstractmethod, abstractproperty

class TemplateBankBase:
    """
    docstring for TemplateBank
    """
    __metaclass__ = ABCMeta

    
    # constructor --------------------------------------------------------------
    def __init__(self, 
                 context, 
                 waveform_length=0, 
                 waveform_delta_x=1,
                 template_table_filename=None,
                 filter_waveform_vector_t=None,
                 pre_condition_vector_t=None,
                 waveform_generator_t=None):

        self.__logger = logging.getLogger('pycbc.TemplateBank')

        # call constructor of <arch>ProcessingObj (2nd parent of TemplateBankBase 
        # derivative TemplatBank<arch>
        super(TemplateBankBase, self).__init__(context) # setting the devicecontext!

        self._sngl_inspiral_table = None        # SnglInspTable from ligolw
        self._sngl_inspiral_table_fname = template_table_filename
        
        self._template_index = 0
        self._template_params = []              # List of tupels of templ params

        self.waveform_length = waveform_length
        self.waveform_delta_x = waveform_delta_x
        
        self.filter_waveform_vector_t = filter_waveform_vector_t
        self.pre_condition_vector_t = pre_condition_vector_t
        
        # setup initial data vectors            
        self.filter_waveform = self.filter_waveform_vector_t(self._devicecontext,
                                                         self.waveform_length,
                                                         self.waveform_delta_x)
        
        
        # try to obtain templates and approximation model from ligolw
        try: 
            if not self._sngl_inspiral_table_fname:
                raise ValueError
                
            self.read_single_inspiral_table(self._sngl_inspiral_table_fname)
        
        except ValueError:
            print "no filename given prototyping simple Tf2 templatebank"
            self._templates_num = 20
            for i in range(self._templates_num):
                tmp = [i*1.0, i*0.5] # m1, m2 very prototyping model of a parameter space              
                self._template_params.append(tmp)
            self.approximation_model= 'Tf2'

        # instanciate waveform generator
        self.waveform_generator= waveform_generator_t(self._devicecontext,
                                                      self.waveform_length,
                                                      self.waveform_delta_x,
                                                      self.approximation_model)

        # get the precondition vector which might be None depending on the 
        # approximation model                                               
        self.precondition_data= self.waveform_generator.perform_generate_precondition(self.pre_condition_vector_t)
        
               
        
                                    
    
    # iterator over the templates in the parameter space -----------------------
    def __iter__(self):
        """
        define the template bank object itself to iterate over it's inherent
        parameter space
        """
        return self

    def next(self):
        if self._template_index == self._templates_num:
            raise StopIteration
        self._template_index = self._template_index + 1
        return self._template_params[self._template_index-1]
   
    # methods ------------------------------------------------------------------
    def read_single_inspiral_table(self, filename, verbose=False):
        """ 
        load the xml file containing the templates and extract the sngl_inspiral table
        """        
        self.__logger.debug("loading templatebank parameterspace: %s"%(filename))
        
        self._sngl_inspiral_table_fname = filename
        xmldoc = ligolw_add.ligolw_add(ligolw.Document(), 
                            [self._sngl_inspiral_table_fname], verbose=verbose)
        
        self._sngl_inspiral_table= table.get_table(xmldoc, 
                                        lsctables.SnglInspiralTable.tableName)
                                        
        self.__logger.debug("obtained %i templates"%(len(self._sngl_inspiral_table)))
        
        self.approximation_model= '' # ToDo - get approximation model from ligolw
        
        self._templates_num = len(self._sngl_inspiral_table)
        self._template_index = 0

        self.__logger.debug("extracting templatebank mass parameters")
        m1s = self._sngl_inspiral_table.get_column('mass1')
        m2s = self._sngl_inspiral_table.get_column('mass2')
        self._template_params = zip(m1s,m2s)


    def precondition_data(self, strain_data):
        """
         
        """
        # if self.precondition_data
        #     strain_data *= self.precondition_data

        pass            
         
            
        
        

