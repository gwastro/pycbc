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

    # constructor -------------------------------------------------------------
    def __init__(self, context, approximation_model=None):

        self.__logger = logging.getLogger('pycbc.WaveFormGeneratorBase')
        self.__logger.debug("instanciate WaveFormGeneratorBase")

        # call constructor of <arch>ProcessingObj (2nd parent of WaveFormGeneratorBase 
        # derivative TemplatBank<arch>
        super(WaveFormGeneratorBase, self).__init__(context)

        self._approximation_model= approximation_model

    @abstractmethod
    def perform_generate_precondition_factor(self, length, delta_x, pre_condition_vector_t):

        pass

    @abstractmethod
    def perform_generate_waveform_filter(self, waveform_filter, **kwargs):

        pass

    @abstractmethod
    def perform_generate_waveform_filter_from_row(self, waveform_filter, template):

        pass

