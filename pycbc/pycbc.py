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
pycbc management and tools
"""

from pycbcopencl import cl_context_t as OpenClContext
from pycbccpu import cpu_context_t as CpuContext

from datavector.datavectorcpu import *
from datavector.datavectoropencl import *    

from abc import ABCMeta, abstractmethod, abstractproperty

import logging
import re

# ------------------- pycbc processing base classes section --------------------

# data_in prototyping and generic ProcessingObj inheritance showcase
# 

# All processing objects have to inherit from this base class via ...
class PyCbcProcessingObj:

    __metaclass__ = ABCMeta

    def __init__(self, device_context):
        self._logger= logging.getLogger('pycbc.pycbc')
        self._devicecontext = device_context

    @abstractmethod
    def data_in(self, datavector):
        pass

# ... their correct derivative according to their processing architecture:
class CpuProcessingObj(PyCbcProcessingObj):

    def __init__(self, device_context):
    
        super(CpuProcessingObj, self).__init__(device_context)

    # data_in() is to be called for every input datavector. 
    # in case of an alien datavector (does not fit to self-architecture) data_in
    # create the new datavector, copies the data by calling the proper transfer 
    # function in the C layer and set new_datavector = old_datavector 
    # (thus destroy the old datavector)

    def data_in(self, datavector):

        #print 'data_in of ' + repr(self) + ' called'
        #print 'with ' + repr(datavector)
        
        if repr(datavector).find("datavectorcpu") >= 0:
            # it is one of us
            return datavector

        else:
            #print 'aliendatavector found:'
            alien_repr_str= repr(datavector)
            #print alien_repr_str
            # find correct new datatype. by parsing alien_repr_str.
            # and instanciate the correct thing, " cloning " from the alien
            new_arch_vector = real_vector_single_cpu_t(len(datavector), datavector.get_delta_x())
            
            # call the transfer function in the C layer
            # prototyping it here:
            for i in range(len(datavector)):
                pass
               # new_arch_vector[i] = datavector[i]
               # fix opencl datavector probs first
            
            return new_arch_vector

class OpenClProcessingObj(PyCbcProcessingObj):


    #
    #
    #  TODO   split into separated modules cpu/cuda/opencl 
    #  to avoid name conflicts with
    #  Ex.: transfer_real_vector_double_from_cpu
    #
    #


    def __init__(self, device_context):
    
        self.__logger= logging.getLogger('pycbc.OpenClProcessingObj')

        super(OpenClProcessingObj, self).__init__(device_context)


    # data_in() is to be called for every input datavector. 
    # in case of an alien datavector (does not fit to self-architecture) data_in
    # create the new datavector, copies the data by calling the proper transfer 
    # function in the C layer and set new_datavector = old_datavector 
    # (thus destroy the old datavector)
    def data_in(self, datavector):

        vector_repr = repr(datavector)
        if (vector_repr.find("datavectoropencl") >= 0):
            self.__logger.debug("data_in found aboriginal datavector {0} thus return".format(vector_repr))
            return datavector

        else:
            # cloning the alien datavector
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

class CpuDeviceContext:

    def __init__(self, devicehandle):
        self.__logger= logging.getLogger('pycbc.pycbc')
        self.__devicehandle = devicehandle

    def __enter__(self):
        self.__logger.debug("__enter__ called ")
        
        # create Cpu device context
        self.__cpucontext = CpuContext(self.__devicehandle)
        
        self.__logger.debug(" On __enter__ create context with {0}:".format(self.__devicehandle)) 
        self.__logger.debug( repr(self.__cpucontext) )
        self.__logger.debug( str(self.__cpucontext) )
        
        return self.__cpucontext  # at with statement binding to as ###
                
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__logger.debug( "__exit__ called " )
        
        # destroy Cpu device context
        del(self.__cpucontext)
         
        
class OpenClDeviceContext:

    def __init__(self, devicehandle):
        self.__logger= logging.getLogger('pycbc.pycbc')
        self.__devicehandle = devicehandle

    def __enter__(self):
        self.__logger.debug("__enter__ called ")
        
        # create OpenCl device context
        self.__openclcontext = OpenClContext(self.__devicehandle)
        
        self.__logger.debug(" On __enter__ create context with {0}:".format(self.__devicehandle)) 
        self.__logger.debug( repr(self.__openclcontext) )
        self.__logger.debug( str(self.__openclcontext) )
        
        return self.__openclcontext  # at with statement binding to as ###
                
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__logger.debug( "__exit__ called " )
        
        # destroy OpenCl device context
        del(self.__openclcontext)
                

class CudaDeviceContext:

    def __init__(self, devicehandle):
        self.__logger= logging.getLogger('pycbc.pycbc')
        self.__devicehandle = devicehandle
        self.__logger.debug("instanciated CudaDeviceContext {0}".format(self.__devicehandle))

    def __enter__(self):
        self.__logger.debug("__enter__ called ")
        
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__logger.debug("__exit__ called ")
        

        
