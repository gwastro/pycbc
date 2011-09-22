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
Datavector stimulator for the pycbc package
transfers datavectors from a CPU to the GPU

"""

# import framework related parent
from pycbc.pycbc_opencl import OpenClProcessingObj


import logging

class DataVecStimOpenCl(OpenClProcessingObj):

    def __init__(self, context):
        self.__logger= logging.getLogger('pycbc.DataVecStimOpenCl')
        self.__logger.debug("instanciated DataVecStimOpenCl") 

        super(DataVecStimOpenCl, self).__init__(context)
