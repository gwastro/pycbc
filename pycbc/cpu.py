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

from pycbc.base import PyCbcProcessingObj
from pycbc.clayer.cpu import cpu_context_t as CpuContext

try:
    from datavector.clayer.opencl import *
except:
    pass
from datavector.clayer.cpu import *    

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



        Currently supporting only OPENCL aliens !!!

        todo -> distinguish by the context!


        """

        vector_repr = repr(datavector)
        print vector_repr
        if (vector_repr.find("datavector.clayer.cpu") >= 0):
            # aboriginal datavector. just return it as it is
            self.__logger.debug("data_in found aboriginal datavector {0} thus return it".format(vector_repr))
            return datavector

        else:

            # cloning the alien datavector (currently only _cpu_t is allowed)
            self.__logger.debug("data_in found alien datavector {0} thus transfer it".format(vector_repr))
            datatype_match= re.search( r'<pycbc.datavector.datavector.clayer.opencl.(.*)(_opencl_t);(.*)', vector_repr)
            tmptype = datatype_match.groups(1)
            datatype = tmptype[0]
            self.__logger.debug("extracted: {0}".format(datatype))

            datavector_tupel= re.split('_+', datatype)

            if datavector_tupel[1] != 'vector':
                raise TypeError('data_in called with something else than a datavector')
        
            if datavector_tupel[0] == 'real':
                if datavector_tupel[2] == 'single':
                    self.__logger.debug("data_in create real_vector_single")
                    self.__logger.debug("data_in create complex_vector_single")
                    new_arch_vector = real_vector_single_cpu_t(len(datavector),
                                                               datavector.get_delta_x())
                    transfer_real_vector_single_to_cpu(self._devicecontext,
                                                         new_arch_vector, datavector)

                elif datavector_tupel[2] == 'double':
                    self.__logger.debug("data_in create real_vector_double")
                    new_arch_vector = real_vector_double_opencl_t(self._devicecontext,
                                                                  len(datavector),
                                                                  datavector.get_delta_x())
                    transfer_real_vector_double_to_cpu(self._devicecontext,
                                                         new_arch_vector, datavector)
                    
                else:
                    raise TypeError('real datavector neither single nor double')

            elif datavector_tupel[0] == 'complex':
                if datavector_tupel[2] == 'single':
                    self.__logger.debug("data_in create complex_vector_single")
                    new_arch_vector = complex_vector_single_opencl_t(self._devicecontext,
                                                                     len(datavector),
                                                                     datavector.get_delta_x())
                    transfer_complex_vector_single_to_cpu(self._devicecontext,
                                                            new_arch_vector, datavector)

                elif datavector_tupel[2] == 'double':
                    self.__logger.debug("data_in create complex_vector_double")
                    new_arch_vector = complex_vector_double_opencl_t(self._devicecontext,
                                                                     len(datavector),
                                                                     datavector.get_delta_x())
                    transfer_complex_vector_double_to_cpu(self._devicecontext,
                                                            new_arch_vector, datavector)

                else:
                    raise TypeError('complex datavector neither single nor double')
            else:
                raise TypeError('datavector neither real nor complex')
    
            return new_arch_vector


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
                
