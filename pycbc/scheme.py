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


import pycbc 


current_context = None

class CPUScheme(object):
    def __enter__(self):
        global current_context 
        
        if type(current_context) is not CPUScheme:
            raise RuntimeError("Nesting Contexts is not allowed")
        
        current_context=self
        return self

    def __exit__(self,type,value,tracebook):
        pass


class CUDAScheme(object):
    def __init__(self,device_id=None):
        self.device_context = None
        self.device = None
        self.device_id=device_id
        import pycuda.autoinit      
                
    
    def __enter__(self):
        global current_context
    
        if type(current_context) is not CPUScheme:
            raise RuntimeError("Nesting is not allowed")
        
        current_context=self     
        return self

    def __exit__(self,type,value,tracebook):
        global current_context
        current_context = CPUScheme()


class OpenCLScheme(object):

    def __init__(self,platform_id=0,device_id=0):
        import pyopencl
        
        self.platform = pyopencl.get_platforms()[platform_id]
        self.device =  self.platform.get_devices()[device_id]
    
        self.device_context = pyopencl.Context(self.platform.get_devices())
        self.queue = pyopencl.CommandQueue(self.device_context)   

    def __enter__(self):
        global current_context
        if type(current_context) is not CPUScheme:
            raise RuntimeError("Nesting is not allowed")

        current_context=self
        
        return self

    def __exit__(self,type,value,tracebook):
        global current_context
        current_context = CPUScheme()
    
#Set the default Processing Context to be the CPU
current_context=CPUScheme()

