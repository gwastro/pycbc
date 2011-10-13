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
pyCBC OpenCl processing object - base class and context
"""

from pycbc_base import PyCbcProcessingObj
from pycbcopencl import cl_context_t as OpenClContext

from datavector.datavectoropencl import *    

import logging
import re

# ------------------- architecture dependent processing base classes section ---

class OpenClProcessingObj(PyCbcProcessingObj):

    def __init__(self, device_context):
    
        self.__logger= logging.getLogger('pycbc.OpenClProcessingObj')

        super(OpenClProcessingObj, self).__init__(device_context)
        
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
        @rtype  snr:   datavector_opencl_t
        @return snr:   OpenCl datavector
        """

        vector_repr = repr(datavector)
        if (vector_repr.find("datavectoropencl") >= 0):
            # aboriginal datavector. just return it as it is
            self.__logger.debug("data_in found aboriginal datavector {0} thus return it".format(vector_repr))
            return datavector

        else:
            # cloning the alien datavector (currently only _cpu_t is allowed)
            self.__logger.debug("data_in found alien datavector {0} thus transfer it".format(vector_repr))
            datatype_match= re.search( r'<pycbc.datavector.datavectorcpu.(.*)(_cpu_t);(.*)', vector_repr)
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
                    new_arch_vector = real_vector_single_opencl_t(self._devicecontext,
                                                                  len(datavector),
                                                                  datavector.get_delta_x())
                    transfer_real_vector_single_from_cpu(self._devicecontext,
                                                         new_arch_vector, datavector)

                elif datavector_tupel[2] == 'double':
                    self.__logger.debug("data_in create real_vector_double")
                    new_arch_vector = real_vector_double_opencl_t(self._devicecontext,
                                                                  len(datavector),
                                                                  datavector.get_delta_x())
                    transfer_real_vector_double_from_cpu(self._devicecontext,
                                                         new_arch_vector, datavector)
                    
                else:
                    raise TypeError('real datavector neither single nor double')

            elif datavector_tupel[0] == 'complex':
                if datavector_tupel[2] == 'single':
                    self.__logger.debug("data_in create complex_vector_single")
                    new_arch_vector = complex_vector_single_opencl_t(self._devicecontext,
                                                                     len(datavector),
                                                                     datavector.get_delta_x())
                    transfer_complex_vector_single_from_cpu(self._devicecontext,
                                                            new_arch_vector, datavector)

                elif datavector_tupel[2] == 'double':
                    self.__logger.debug("data_in create complex_vector_double")
                    new_arch_vector = complex_vector_double_opencl_t(self._devicecontext,
                                                                     len(datavector),
                                                                     datavector.get_delta_x())
                    transfer_complex_vector_double_from_cpu(self._devicecontext,
                                                            new_arch_vector, datavector)

                else:
                    raise TypeError('complex datavector neither single nor double')
            else:
                raise TypeError('datavector neither real nor complex')
    
            return new_arch_vector

# ------------------- device context section -----------------------------------
        
class OpenClDeviceContext:

    def __init__(self, devicehandle):
        self.__logger= logging.getLogger('pycbc.OpenClDeviceContext')
        self.__devicehandle = devicehandle

    def __enter__(self):
        self.__logger.debug("__enter__ called ")
        
        # create OpenCl device context
        self.__openclcontext = OpenClContext(self.__devicehandle)
        
        self.__logger.debug(" On __enter__ create context for device {0}:".format(self.__devicehandle)) 
        self.__logger.debug( repr(self.__openclcontext) )
        self.__logger.debug( str(self.__openclcontext) )
        
        return self.__openclcontext # at the "with" statement binding to "as" context
                
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__logger.debug( "__exit__ called " )
        
        # destroy OpenCl device context clayer member
        del(self.__openclcontext)
        self.__logger.debug(" On __exit__ destroyed openclcontext of device {0}:".format(self.__devicehandle))
                
