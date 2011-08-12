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

import logging

class CpuDeviceContext:

    def __init__(self, devicehandle):
        self.__logger= logging.getLogger('pycbc.pycbc')
        self.__devicehandle = devicehandle
        self.__logger.debug("instanciated CpuDeviceContext {0}".format(self.__devicehandle))


    def __enter__(self):
        self.__logger.debug("__enter__ called ")
        
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__logger.debug("__exit__ called ")


class OpenClDeviceContext:

    def __init__(self, devicehandle):
        self.__logger= logging.getLogger('pycbc.pycbc')
        self.__devicehandle = devicehandle
        self.__logger.debug("instanciated OpenClDeviceContext {0}".format(self.__devicehandle))


    def __enter__(self):
        self.__logger.debug("__enter__ called ")
        
        # create OpenCl device context
        self.__openclcontext = OpenClContext(1)
        self.__logger.debug( repr(self.__openclcontext) )
        self.__logger.debug( str(self.__openclcontext) )
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        print "__exit__ called for OpenClDeviceContext "
        self.__logger.debug("__exit__ called ")
        
        # destroy physical device context
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
        

        