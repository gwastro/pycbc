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
This modules provides a device independent Array class based on PyCUDA, PyOpenCL, and 
Numpy. 
"""

import numpy
import pycbc.processingcontext

__author__ = "Alex Nitz <alex.nitz@ligo.org>"

# Check for optional components
have_cuda=False
have_opencl=False

try:
    import pycuda
    import pycuda.gpuarray
    import pycuda.autoinit
    have_cuda=True
except ImportError:
    have_cuda=False
    
try:
    import pyopencl
    import pyopencl.array
    have_opencl=True
except ImportError:
    have_opencl=False
    

class Array(object):
    def __init__(self,object, dtype=None, copy=True,context=None):
    
        #The private members of the class. No guarantee is made about the behavior
        # of these variables
        self._context=None
        self._data=None

        
        #Determine the correct context for the new Array instance
        if context:
            self._context=context
        else:
            self._context=pycbc.processingcontext.current_context

        if copy is False:
            if type(object) is Array:
                self._data=object._data
            elif type(object) is numpy.ndarray:
                self._data = object
            elif have_cuda and type(object) is pycuda.gpuarray.GPUArray:
                self._data = object
            elif have_opencl and type(object) is pyopencl.array.Array:
                self._data = object
            else:
                raise TypeError, str(type(object))+' is not supported'
                
        if copy is True:
        
            input_data = None
            if type(object) is Array:
                input_data = object._data
            else:
                input_data=object

            if type(self._context) is pycbc.processingcontext.CPUContext:
                if have_cuda and type(input_data) is pycuda.gpuarray.GPUArray:
                    self._data = input_data.get()
                elif have_opencl and type(input_data) is pyopencl.array.Array:
                    self._data = input_data.get() 
                else:
                    self._data = numpy.array(input_data,dtype = dtype)           
                    
            elif have_cuda and type(self._context) is pycbc.processingcontext.CUDAContext:
                if type(input_data) is pycuda.gpuarray.GPUArray:
                    self._data = pycuda.gpuarray.zeros(input_data.size,input_data.dtype)
                    pycuda.driver.memcpy_dtod(self._data.gpudata,input_data.gpudata,self._data.nbytes)
                else:
                    self._data = pycuda.gpuarray.to_gpu(numpy.array(input_data,dtype=dtype))
                    
            elif have_opencl and type(self._context) is pycbc.processingcontext.OpenCLContext:
                if type(input_data) is pyopencl.array.Array:
                    self._data = pyopencl.array.zeros(self._context.queue,input_data.size,input_data.dtype)
                    pyopencl.enqueue_copy(self._context.queue,self._data.data,input_data.data)
                else:
                    self._data = pyopencl.array.to_device(self._context.queue,numpy.array(input_data,dtype=dtype))
            
            else:
                raise RuntimeError   
                                    
        
    def _returnarray(fn):
        def wrapped(*args):
            return Array(fn(*args),copy=False)
        return wrapped
    
    def _checkother(fn):
        def checked(self,other):
            # TODO Check if other is compatible array (dtype)
            # TODO Check if other is compatible scaler type
            if type(other) is Array:
                other=other._data
            return fn(self,other)
        return checked
            
    @_returnarray
    @_checkother 
    def __mul__(self,other):
        return self._data * other
    
    @_checkother
    @_returnarray
    def __rmul__(self,other):
        return self._data * other

    @_checkother
    def __imul__(self,other):
        self._data *= other
        return self
    
    @_checkother
    @_returnarray
    def __add__(self,other):
        return self._data + other
        
    @_checkother
    @_returnarray
    def __radd__(self,other):
        return self._data + other
    
    @_checkother 
    def __iadd__(self,other):
        self._data += other
        return self
        
    @_checkother
    @_returnarray    
    def __div__(self,other):
        return self._data / other
        
    @_checkother
    @_returnarray
    def __rdiv__(self,other):
        return other / self._data
        
    @_checkother
    def __idiv__(self,other):
        self._data /= other
        return self
        
    @_checkother
    @_returnarray
    def __sub__(self,other):
        return self._data - other 
        
    @_checkother
    @_returnarray
    def __rsub__(self,other):
        return other - self._data
        
    @_checkother
    def __isub__(self,other):
        self._data -= other
        return self
        
    @_checkother
    @_returnarray
    def __pow__(self,other):
        return self._data ** other
    
    @_returnarray  
    def __abs__(self):
        return abs(self._data)
        
    def __len__(self):
        return len(self._data)
        
    def __str__(self):
        return str(self._data)
    
    @_returnarray
    def real(self):
        return Array(self._data.real,copy=True,context=self._context)

    @_returnarray
    def imag(self):
        return Array(self._data.imag,copy=True,context=self._context)
    
    @_returnarray
    def conj(self):
        return self._data.conj()
        
    def norm(self):
        pass
    
    def innerprod(self,other):
        pass
        
    def ptr(self):
        if type(_data) is numpy.ndarray:
            return _data.ctypes.get_data()
            
        if have_cuda and type(_data) is pycuda.gpuarray.GPUArray:
            return _data.ptr
    
        if have_opencl and type(_data) is pyopencl.array.Array:
            return _data.data
        
    
    
    
    
__all__ = ["Array"]    
