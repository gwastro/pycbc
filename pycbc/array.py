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
This modules provides a device independent Array class based on PyCUDA,
PyOpen CL, and Numpy. 
"""

import functools as _functools

import swiglal as _swiglal
import numpy as _numpy
from numpy import float32,float64,complex64,complex128

import pycbc as _pycbc
import pycbc.scheme as _scheme

if _pycbc.HAVE_CUDA:
    import pycuda as _pycuda
    import pycuda.gpuarray as _cudaarray
if _pycbc.HAVE_OPENCL:
    import pyopencl as _pyopencl
    import pyopencl.array as _openclarray 

_ALLOWED_DTYPES = [_numpy.float32,_numpy.float64,_numpy.complex64,
                   _numpy.complex128]
_ALLOWED_SCALARS = [int,long, float]


def _convert_to_scheme(ary):
    if ary._scheme is not _scheme.mgr.state: 
        converted_array = Array(ary,dtype=ary._data.dtype) 
        ary._data = converted_array._data
        ary._scheme = _scheme.mgr.state


class Array(object):
    """Array used to do numeric calculations on a various compute
    devices. It is a convience wrapper around _numpy, _pyopencl, and
    _pycuda. 
    """

    def __init__(self,initial_array, dtype=None, copy=True):
        """ initial_array: An array-like object as specified by NumPy, this 
        also includes instances of an underlying data type as described in 
        section 3 or an instance of the PYCBC Array class itself. This
        object is used to populate the data of the array.

        dtype: A NumPy style dtype that describes the type of 
        encapsulated data (float32,compex64, etc)

        copy: This defines whether the initial_array is copied to instantiate 
        the array or is simply referenced. If copy is false, new data is not 
        created, and so all arguments that would force a copy are ignored. 
        The default is to copy the given object. 
        """
        self._scheme=_scheme.mgr.state
        self._data = None
        
        if not copy:
            if isinstance(initial_array,Array):
                if self._scheme is not initial_array._scheme:
                    raise TypeError("Cannot avoid a copy of this Array")
                self._data = initial_array._data
            elif type(initial_array) is _numpy.ndarray:
                self._data = initial_array
            elif type(self._scheme) is _scheme.CUDAScheme  and \
            type(initial_array) is _cudaarray.GPUArray:
                self._data = initial_array
            elif type(self._scheme) is _scheme.OpenCLScheme and \
            type(initial_array) is _openclarray.Array:
                self._data = initial_array
            else:
                raise TypeError(str(type(initial_array))+' is not supported')
                
            # Check that the dtype is supported.
            if self._data.dtype not in _ALLOWED_DTYPES:
                raise TypeError(str(self._data.dtype) + ' is not supported')
            
                
        if copy:                  
            #Check that a valid dtype was given.
            if dtype is not None: 
                if dtype not in _ALLOWED_DTYPES:
                    raise TypeError(str(dtype) + ' is not supported')   
        
            #Unwrap initial_array
            input_data = None
            if isinstance(initial_array,Array):
                input_data = initial_array._data
            else:
                input_data = initial_array

            #Create new instance with input_data as initialization.
            if self._scheme is None:
                if _pycbc.HAVE_CUDA and (type(input_data) is
                                         _cudaarray.GPUArray):
                    self._data = input_data.get()
                elif _pycbc.HAVE_OPENCL and (type(input_data) is
                                             _openclarray.Array):
                    self._data = input_data.get()
                else:
                    self._data = _numpy.array(input_data,dtype=dtype)          
            elif _pycbc.HAVE_CUDA and (type(self._scheme) is 
                                       _scheme.CUDAScheme):
                if type(input_data) is _cudaarray.GPUArray:
                    input_data = input_data.get()
                self._data = _cudaarray.to_gpu(_numpy.array(input_data,
                                                            dtype=dtype))       
            elif _pycbc.HAVE_OPENCL and (type(self._scheme) is 
                                         _scheme.OpenCLScheme):
                if type(input_data) is _openclarray.Array:
                    input_data = input_data.get()
                self._data = _openclarray.to_device(self._scheme.queue,
                                                    _numpy.array(input_data,
                                                                 dtype=dtype))
            else:
                raise TypeError('Invalid Processing scheme Type')

            #If no dtype was given, default to Double, and Double Complex. 
            if dtype is None:
                if self._data.dtype.kind == 'c':
                    self._data = self._data.astype(_numpy.complex128)
                else:
                    self._data = self._data.astype(_numpy.float64)            
        
    def _returnarray(fn):
        # Wrapper for method functions to return a PyCBC Array class 
        @_functools.wraps(fn)
        def wrapped(*args):
            return Array(fn(*args),copy=False)
        return wrapped
    

    
    def _checkother(fn):
        # Checks the input to method functions 
        @_functools.wraps(fn)
        def checked(self,other):
            if type(other) not in _ALLOWED_SCALARS and not Array:
                raise TypeError(str(type(other)) + ' is incompatible with ' +
                                str(type(self)))          
            if isinstance(other,Array):
                _convert_to_scheme(other)    
                other = other._data
            return fn(self,other)
        return checked  
            
    def _convert(fn):
        # Convert this array to the current processing scheme
        @_functools.wraps(fn)
        def converted(self,*args):
            _convert_to_scheme(self)
            return fn(self,*args)
        return converted
            
    @_returnarray
    @_convert
    @_checkother 
    def __mul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        return self._data * other
    
    @_returnarray
    @_convert
    @_checkother
    def __rmul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        return self._data * other

    @_convert
    @_checkother
    def __imul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        self._data *= other
        return self
    
    @_returnarray
    @_convert
    @_checkother
    def __add__(self,other):
        """ Add Array to Array or scalar and return an Array. """
        return self._data + other
        
    @_returnarray
    @_convert
    @_checkother
    def __radd__(self,other):
        """ Add Array to Array or scalar and return an Array. """
        return self._data + other
    
    @_convert
    @_checkother
    def __iadd__(self,other):
        """ Add Array to Array or scalar and return an Array. """
        self._data += other
        return self
        
    @_convert
    @_checkother
    @_returnarray    
    def __div__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        return self._data / other
        
    @_returnarray
    @_convert
    @_checkother
    def __rdiv__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        return other / self._data
        
    @_convert
    @_checkother
    def __idiv__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        self._data /= other
        return self
        
    @_returnarray
    @_convert
    @_checkother
    def __sub__(self,other):
        """ Subtract Array or scalar from Array and return an Array. """
        return self._data - other 
        
    @_returnarray
    @_convert
    @_checkother
    def __rsub__(self,other):
        """ Subtract Array or scalar from Array and return an Array. """
        return other - self._data
        
    @_convert   
    @_checkother
    def __isub__(self,other):
        """ Subtract Array or scalar from Array and return an Array. """
        self._data -= other
        return self
        
    @_returnarray
    @_convert
    @_checkother
    def __pow__(self,other):
        """ Exponentiate Arrary by scalar """
        return self._data ** other
    
    @_returnarray  
    @_convert
    def __abs__(self):
        """ Return absolute value of Array """
        return abs(self._data)
        
    def __len__(self):
        """ Return length of Array """
        return len(self._data)
        
    def __str__(self):
        return str(self._data)
    
    @_returnarray
    @_convert
    def real(self):
        """ Return real part of Array """
        return Array(self._data.real,copy=True)

    @_returnarray
    @_convert
    def imag(self):
        """ Return imaginary part of Array """
        return Array(self._data.imag,copy=True)
    
    @_returnarray
    @_convert
    def conj(self):
        """ Return complex conjugate of Array. """ 
        return self._data.conj()
        
    @_convert
    def sum(self):
        """ Return the sum of the the array. """
        if type(self._data) is _numpy.ndarray:
            return self._data.sum()
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return _pycuda.gpuarray.sum(self._data).get().max()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return _pyopencl.array.sum(self._data).get().max()
    
    @_convert
    @_checkother
    def dot(self,other):
        """ Return the dot product"""
        if type(self._data) is _numpy.ndarray:
            return self._data.dot(other)
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return _pycuda.gpuarray.dot(self._data,other).get().max()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return _pyopencl.array.dot(self._data,other).get().max()
    
    @_convert
    def __getitem__(self, index):
        if isinstance(index, slice):
            return Array(self._data[index],copy=False)
        else:
            if type(self._data) is _numpy.ndarray:
                return self._data[index]
            elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
                return self._data.get()[index]
            elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
                return self._data.get()[index]

    @_convert
    def ptr(self):
        """ Returns a pointer to the memory of this array """
        if type(self._data) is _numpy.ndarray:
            raise TypeError("Please use swiglal for CPU objects")
        if _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return self._data.ptr
        if _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return self._data.data
        
        
    @_convert
    def  lal(self):
        """ Returns a LAL Object that contains this data """
        lal_data = None
        if type(self._data) is not _numpy.ndarray:
            raise TypeError("Cannot return lal type from the GPU") 
        elif self._data.dtype == float32:
            lal_data = _swiglal.XLALCreateREAL4Vector(len(self))
        elif self._data.dtype == float64:
            lal_data = _swiglal.XLALCreateREAL8Vector(len(self))
        elif self._data.dtype == complex64:
            lal_data = _swiglal.XLALCreateCOMPLEX8Vector(len(self))
        elif self._data.dtype == complex128:
            lal_data = _swiglal.XLALCreateCOMPLEX16Vector(len(self))    

        lal_data.data = self._data
        return lal_data

def zeros(length,dtype=None):
	return Array(_numpy.zeros(length),dtype=None)



