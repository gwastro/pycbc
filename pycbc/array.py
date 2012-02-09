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

import numpy as _numpy
import pycbc as _pycbc
import pycbc.processingcontext as _processingcontext
import functools as _functools


__author__ = "Alex Nitz <alex.nitz@ligo.org>"

if _pycbc.have_cuda:
    import pycuda as _pycuda
    import pycuda.gpuarray as _cudaarray

if _pycbc.have_opencl:
    import pyopencl as _pyopencl
    import pyopencl.array as _openclarray
    


# Supported Types 
_allowed_dtypes = [_numpy.float32,_numpy.float64,_numpy.complex64,_numpy.complex128]
_allowed_scalars = [int,long, float]


from numpy import float32,float64,complex64,complex128

class Array(object):
    """Array used to do numeric calculations on a various compute devices. It is a 
    convience wrapper around _numpy, _pyopencl, and _pycuda. 

    object: An array-like object as specified by NumPy, this also includes instances 
    of an underlying data type as described in section 3 or an instance of the PYCBC
    Array class itself. This object is used to populate the data of the array.

    dtype: A NumPy style dtype that describes the type of 
    encapsulated data (float32,compex64, etc)

    copy: This defines whether the object is copied to instantiate the array 
    or is simply referenced. If copy is false, new data is not created, and so the 
    context is ignored. The default is to copy the given object.

    context: This one of the contexts found in _processingcontext
    """
    def __init__(self,object, dtype=None, copy=True,context=None):
    
        #The private members of the class. No guarantee is made about the behavior
        # of these variables
        self._context=None
        self._data=None

        
        #Determine the correct context for the new Array instance
        if context:
            self._context=context
        else:
            self._context=_processingcontext.current_context

        if copy is False:
            if type(object) is Array:
                self._data=object._data
            elif type(object) is _numpy.ndarray:
                self._data = object
            elif _pycbc.have_cuda and type(object) is _cudaarray.GPUArray:
                self._data = object
            elif _pycbc.have_opencl and type(object) is _openclarray.Array:
                self._data = object
            else:
                raise TypeError, str(type(object))+' is not supported'
                
            # Check that the dtype is supported
            if self._data.dtype not in _allowed_dtypes:
                raise TypeError, str(self._data.dtype) + ' is not supported' 
            
                
        if copy is True:
        
            #Check that a valid dtype was given
            if dtype: 
                if dtype not in _allowed_dtypes:
                    raise TypeError, str(dtype) + ' is not supported'   
        
            #Unwrap object    
            input_data = None
            if type(object) is Array:
                input_data = object._data
            else:
                input_data=object

            #Create new instance with input_data as initialization
            if type(self._context) is _processingcontext.CPUContext:
                if _pycbc.have_cuda and type(input_data) is _cudaarray.GPUArray:
                    self._data = input_data.get()
                elif _pycbc.have_opencl and type(input_data) is _openclarray.Array:
                    self._data = input_data.get()
                else:
                    self._data = _numpy.array(input_data,dtype=dtype)     
  
                    
            elif _pycbc.have_cuda and type(self._context) is _processingcontext.CUDAContext:
                if type(input_data) is _cudaarray.GPUArray:
                    input_data = input_data.get()

                self._data = _cudaarray.to_gpu(_numpy.array(input_data,dtype=dtype))
                    
            elif _pycbc.have_opencl and type(self._context) is _processingcontext.OpenCLContext:
                if type(input_data) is _openclarray.Array:
                    input_data = input_data.get()
                    
                self._data = _openclarray.to_device(self._context.queue,_numpy.array(input_data,dtype=dtype))
            
            else:
                raise TypeError, ' Invalid Processing Context Type'  
                
                   
            #If no dtype was not given, default to Double, and Double Complex 
            if not dtype:
                if self._data.dtype.kind == 'c':
                    self._data = self._data.astype(_numpy.complex128)
                else:
                    self._data = self._data.astype(_numpy.float64)
                
        
    def _returnarray(fn):
        """ Wrapper for method functions to return a PyCBC Array class """
        @_functools.wraps(fn)
        def wrapped(*args):
            return Array(fn(*args),copy=False)
        return wrapped
    
    def _checkother(fn):
        """ Checks the input to method functions """
        @_functools.wraps(fn)
        def checked(self,other):

            if type(other) not in _allowed_scalars and not Array:
                raise TypeError,str( type(other) ) + ' is incompatible with ' + str( type(self)  ) 
                
            if type(other) is Array:
                if type(self._context) is not type(other._context):
                    raise TypeError, "Incompatible Contexts"
                
                if self._data.dtype is not other._data.dtype:
                    raise TypeError, "dtypes do not match"
 
                other=other._data
            return fn(self,other)
        return checked
            
    @_returnarray
    @_checkother 
    def __mul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        return self._data * other
    
    @_checkother
    @_returnarray
    def __rmul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        return self._data * other

    @_checkother
    def __imul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        self._data *= other
        return self
    
    @_checkother
    @_returnarray
    def __add__(self,other):
        """ Add Array to Array or scalar and return an Array. """
        return self._data + other
        
    @_checkother
    @_returnarray
    def __radd__(self,other):
        """ Add Array to Array or scalar and return an Array. """
        return self._data + other
    
    @_checkother 
    def __iadd__(self,other):
        """ Add Array to Array or scalar and return an Array. """
        self._data += other
        return self
        
    @_checkother
    @_returnarray    
    def __div__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        return self._data / other
        
    @_checkother
    @_returnarray
    def __rdiv__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        return other / self._data
        
    @_checkother
    def __idiv__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        self._data /= other
        return self
        
    @_checkother
    @_returnarray
    def __sub__(self,other):
        """ Subtract Array or scalar from Array and return an Array. """
        return self._data - other 
        
    @_checkother
    @_returnarray
    def __rsub__(self,other):
        """ Subtract Array or scalar from Array and return an Array. """
        return other - self._data
        
    @_checkother
    def __isub__(self,other):
        """ Subtract Array or scalar from Array and return an Array. """
        self._data -= other
        return self
        
    @_checkother
    @_returnarray
    def __pow__(self,other):
        """ Exponentiate Arrary by scalar """
        return self._data ** other
    
    @_returnarray  
    def __abs__(self):
        """ Return absolute value of Array """
        return abs(self._data)
        
    def __len__(self):
        """ Return length of Array """
        return len(self._data)
        
    def __str__(self):
        return str(self._data)
    
    @_returnarray
    def real(self):
        """ Return real part of Array """
        return Array(self._data.real,copy=True,context=self._context)

    @_returnarray
    def imag(self):
        """ Return imaginary part of Array """
        return Array(self._data.imag,copy=True,context=self._context)
    
    @_returnarray
    def conj(self):
        """ Return complex conjugate of Array """ 
        return self._data.conj()
        
    def norm(self):
        pass
    
    def innerprod(self,other):
        pass
    
    def __getitem__(self, index):
        if isinstance(index, slice):
            return Array(self._data[index],copy=False)
        else:
            if type(self._data) is _numpy.ndarray:
                return self._data[index]
            elif _pycbc.have_cuda and type(self._data) is _cudaarray.GPUArray:
                return self._data.get()[index]
            elif _pycbc.have_opencl and type(self._data) is _openclarray.Array:
                return self._data.get()[index]
        
    def ptr(self):
        """ Return pointer to memory """
        if type(self._data) is _numpy.ndarray:
            return self._data.ctypes.get_data()
            
        if _pycbc.have_cuda and type(self._data) is _cudaarray.GPUArray:
            return self._data.ptr
    
        if _pycbc.have_opencl and type(self._data) is _openclarray.Array:
            return self._data.data
        
   
