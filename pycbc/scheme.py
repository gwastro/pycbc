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

class _SchemeManager(object):

    _single = None

    def __init__(self):

        if _SchemeManager._single is not None: 
            raise RuntimeError("SchemeManager is a private class")
        _SchemeManager._single= self

        self.state= None  
        self._lock= False
        
    def lock(self):
        self._lock= True 
        
    def unlock(self):
        self._lock= False
    
    def shift_to(self, state):
        if self._lock is False:
            self.state = state
        else:
            raise RuntimeError("The state is locked, cannot shift schemes")

# Create the global processing scheme manager
mgr = _SchemeManager()


class DefaultScheme(object):
    """Context that sets PyCBC objects to use CPU processing. """
    _single = None
    def __init__(self):
        if DefaultScheme._single is not None:
            raise RuntimeError("Only one processing scheme can be used")
        DefaultScheme._single = True
    def __enter__(self):
        if type(self) is DefaultScheme:
            mgr.shift_to(None)
        else:
            mgr.shift_to(self)
        mgr.lock()
    def __exit__(self,type,value,traceback):
        mgr.unlock()
        mgr.shift_to(None)
        
    def __del__(self):
        DefaultScheme._single = None

_cuda_cleanup_list=[]

def register_clean_cuda(function):
    _cuda_cleanup_list.append(function)

def clean_cuda(context):
    #Before cuda context is destroyed, all item destructions dependent on cuda must take place
    #This calls all functions that have been registered with _register_clean_cuda() in reverse order
    #So the last one registered, is the first one cleaned
    _cuda_cleanup_list.reverse()
    for func in _cuda_cleanup_list:
        func()
    
    context.pop()
    from pycuda.tools import clear_context_caches
    clear_context_caches()

class CUDAScheme(DefaultScheme):
    """Context that sets PyCBC objects to use a CUDA processing scheme. """
    def __init__(self,device_num=0):
        DefaultScheme.__init__(self)
        if not pycbc.HAVE_CUDA:
            raise RuntimeError("Install PyCUDA to use CUDA processing")    
        import pycuda.driver   
        pycuda.driver.init()
        self.device = pycuda.driver.Device(device_num)
        self.context = self.device.make_context()
        import atexit
        atexit.register(clean_cuda,self.context)
        


class OpenCLScheme(DefaultScheme):
    """Context that sets PyCBC objects to use a OpenCL processing scheme. """
    _single = None
    def __init__(self,platform_name=None,device_num=0):
        DefaultScheme.__init__(self)
        if not pycbc.HAVE_OPENCL:
            raise RuntimeError("Install PyOpenCL to use OpenCL processing")   
        import pyopencl  
        
        #If no platform is given, use the first one
        if platform_name is None:
            platform_id = 0
        elif platform_name is not None:
            for platform in pyopencl.get_platforms():
                if platform.name == platform_name:
                    platform_id = pyopencl.get_platforms().index(platform)
        
        self.platform = pyopencl.get_platforms()[platform_id]
        self.device =  self.platform.get_devices()[device_num]
        self.context = pyopencl.Context([self.device])
        self.queue = pyopencl.CommandQueue(self.context)     



