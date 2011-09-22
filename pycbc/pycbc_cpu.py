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
pyCBC Cpu processing object - base class and context
"""

from pycbc_base import PyCbcProcessingObj
from pycbccpu import cpu_context_t as CpuContext

from datavector.datavectorcpu import *    

import logging
import re

# ------------------- architecture dependent processing base classes section ---

class CpuProcessingObj(PyCbcProcessingObj):

    def __init__(self, device_context):
    
        self.__logger= logging.getLogger('pycbc.CpuProcessingObj')

        super(CpuProcessingObj, self).__init__(device_context)
        
    def data_in(self, datavector):
        """
        this method has to be called by a pyCBC processing object 
        for/with every incoming input-datavector. In case of an "alien"
        datavector (which does not fit to the architecture of the calling)
        pyCBC processing object) the method instanciates the correct
        datavector (an "aboriginal" datavector), copies the data or transfer it 
        to the GPU and return a reference of the new datavector. The old 
        datavector will be dereferenced so it will be deleted by the 
        garbage collection
        @type  datavector: datavector_<arch>_t
        @param datavector: any datavector
        @rtype  snr:   datavector_cpu_t
        @return snr:   Cpu datavector
        """

        vector_repr = repr(datavector)
        if (vector_repr.find("datavectorcpu") >= 0):
            # aboriginal datavector. just return it as it is
            self.__logger.debug("data_in found aboriginal datavector {0} thus return it".format(vector_repr))
            return datavector

        else:
        
            raise TypeError('data_in cpu has to be implemented to transfer alien datavectors to the cpu')
        
            return none
        


# ------------------- device context section -----------------------------------

class CpuDeviceContext:

    def __init__(self, devicehandle):
        self.__logger= logging.getLogger('pycbc.CpuDeviceContext')
        self.__devicehandle = devicehandle

    def __enter__(self):
        self.__logger.debug("__enter__ called ")
        
        # create Cpu device context
        self.__cpucontext = CpuContext(self.__devicehandle)
        
        self.__logger.debug(" On __enter__ create context for device {0}:".format(self.__devicehandle)) 
        self.__logger.debug( repr(self.__cpucontext) )
        self.__logger.debug( str(self.__cpucontext) )
        
        return self.__cpucontext  # at the "with" statement binding to "as" context
                
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__logger.debug( "__exit__ called " )
        
        # destroy Cpu device context clayer member
        del(self.__cpucontext)
        self.__logger.debug(" On __exit__ destroyed cpucontext of device {0}:".format(self.__devicehandle))
                
