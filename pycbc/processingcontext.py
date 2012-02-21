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

class CPUContext(object):
    def __init__(self):
        self.device_context = None
        self.prior_context = None

    def __enter__(self):
        global current_context 
        
        if type(current_context) is not CPUContext:
            raise RuntimeError("Nesting Contexts is not allowed")
        
        current_context=self
        return self

    def __exit__(self,type,value,tracebook):
        pass


class CUDAContext(object):
    def __init__(self,device_id=None):
        self.device_context = None
        self.device = None
        self.device_id=device_id
        import pycuda.autoinit
                
    
    def __enter__(self):
    
        global current_context
        if type(current_context) is not CPUContext:
            raise RuntimeError("Nesting Contexts is not allowed")
        
        current_context=self 
        return self

    def __exit__(self,type,value,tracebook):

    
        global current_context
        current_context = CPUContext()

class OpenCLContext(object):

    def __init__(self,device_id=None):
        self.device_context = None
        self.device_id=device_id
        import pyopencl
        self.device_context = pyopencl.create_some_context()
        self.queue = pyopencl.CommandQueue(self.device_context)   
        pass

    def __enter__(self): 


        global current_context
        current_context=self
        
        return self

    def __exit__(self,type,value,tracebook):
        global current_context
        current_context = CPUContext()
    
#Set the default Processing Context to be the CPU
current_context=CPUContext()

