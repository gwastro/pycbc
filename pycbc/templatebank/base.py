# Copyright (C) 2011 Karsten Wiesner, Drew Keppel
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

    # constructor -------------------------------------------------------------
    def __init__(self, context, waveform_length=0, waveform_delta_x=1,
            template_table_filename=None, waveform_approximant=None,
            filter_waveform_vector_t=None, pre_condition_vector_t=None,
            waveform_generator_t=None):

        self.__logger = logging.getLogger('pycbc.TemplateBankBase')
        self.__logger.debug("instanciate TemplateBankBase")

        # call constructor of <arch>ProcessingObj (2nd parent of TemplateBankBase 
        # derivative TemplatBank<arch>
        super(TemplateBankBase, self).__init__(context) # setting the devicecontext!

        self._sngl_inspiral_table = None        # SnglInspTable from ligolw
        self._sngl_inspiral_table_fname = template_table_filename

        self._template_index = 0

        self.waveform_length = waveform_length
        self.waveform_delta_x = waveform_delta_x
        
        self.filter_waveform_vector_t = filter_waveform_vector_t
        self.pre_condition_vector_t = pre_condition_vector_t

        # setup initial data vectors            
        self.waveform_filter = self.filter_waveform_vector_t(self._devicecontext,
            self.waveform_length,
            self.waveform_delta_x)
        

        # try to obtain templates and approximation model from ligolw
        try: 
            if not self._sngl_inspiral_table_fname:
                raise ValueError

            self.read_single_inspiral_table(self._sngl_inspiral_table_fname)

        except ValueError:
            self.__logger.debug("no filename given")
            raise ValueError

        # override template bank approximant if given to constructor
        if waveform_approximant:
            self.approximation_model = waveform_approximant

        # instanciate waveform generator
        self.waveform_generator = waveform_generator_t(self._devicecontext,
            self.waveform_length,
            self.waveform_delta_x,
            self.approximation_model,
            self.waveform_filter)

        # get the precondition vector which might be None depending on the 
        # approximation model   

        self.precondition_factor = \
        self.waveform_generator.perform_generate_precondition( 
            self.pre_condition_vector_t)
                                    

    # iterator over the templates in the parameter space ----------------------
    def __iter__(self):
        """
        define the template bank object itself to iterate over it's inherent
        parameter space
        """
        return self

    def next(self):
        if self._template_index == self._templates_num:
            self._template_index = 0
            raise StopIteration
        self._template_index = self._template_index + 1
        return self._sngl_inspiral_table[self._template_index-1]

    # methods -----------------------------------------------------------------
    def read_single_inspiral_table(self, filename, verbose=False):
        """ 
        load the xml file containing the templates and extract the sngl_inspiral table
        """        
        self.__logger.debug("loading templatebank parameterspace: %s"%(filename))

        self._sngl_inspiral_table_fname = filename
        xmldoc = ligolw_add.ligolw_add(ligolw.Document(), 
            [self._sngl_inspiral_table_fname], verbose=verbose)

        self._sngl_inspiral_table = table.get_table(xmldoc, 
            lsctables.SnglInspiralTable.tableName)

        self.__logger.debug("obtained %i templates"%(len(self._sngl_inspiral_table)))

        process_params = table.get_table(xmldoc, 
            lsctables.ProcessParamsTable.tableName)

        for row in process_params:
            if row.param == '--approximant':
                self.approximation_model = row.value
        if self.approximant is None:
            self.__logger.debug("no approximant found in template bank file")
            raise KeyError

        self._templates_num = len(self._sngl_inspiral_table)
        self._template_index = 0

    def precondition_data(self, strain_data):
        """
        """
        # if self.precondition_data
        #     strain_data *= self.precondition_data

        pass

