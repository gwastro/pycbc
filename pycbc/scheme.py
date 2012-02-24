# Copyright (C) 2012  Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
This modules provides python contexts that set the default behavior for PyCBC
objects. 
"""
import pycbc

class CPUScheme(object):
    """Context that sets PyCBC objects to use CPU processing. """
    def __enter__(self):
        pass
    def __exit__(self,type,value,traceback):
        pass

class _DeviceScheme(object):
    def __enter__(self):
        global current_context 
        if type(current_context) is not CPUScheme:
            raise RuntimeError("Nesting processing schemes are not supported.")
        
        current_context = self
        return self

    def __exit__(self,type,value,traceback):
        global current_context
        current_context = CPUScheme()

class CUDAScheme(_DeviceScheme):
    """Context that sets PyCBC objects to use a CUDA processing scheme. """
    def __init__(self,device_num=0):
        if not pycbc.HAVE_CUDA:
            raise RuntimeError("Install PyCUDA to use CUDA processing")    
        import pycuda.driver   
        pycuda.driver.init()
        self.device = pycuda.driver.Device(device_num)
        self.context = self.device.make_context()
        self.context.pop()
        
    def __enter__(self):
        _DeviceScheme.__enter__(self)
        self.context.push()
        
    def __exit__(self,type,value,traceback):
        _DeviceScheme.__exit__(self,type,value,traceback)
        self.context.pop()
           
class OpenCLScheme(_DeviceScheme):
    """Context that sets PyCBC objects to use a OpenCL processing scheme. """
    def __init__(self,platform_id=0,device_num=0):
        if not pycbc.HAVE_OPENCL:
            raise RuntimeError("Install PyOpenCL to use OpenCL processing")   
        import pyopencl  
        self.platform = pyopencl.get_platforms()[platform_id]
        self.device =  self.platform.get_devices()[device_num]
        self.device_context = pyopencl.Context(self.platform.get_devices())
        self.queue = pyopencl.CommandQueue(self.device_context)     


#Set the default Processing Context to be the CPU
current_context=CPUScheme()

