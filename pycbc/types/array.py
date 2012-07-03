# Copyright (C) 2012  Alex Nitz, Josh Willis, Andrew Miller
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
PyOpenCL, and Numpy.
"""

import functools as _functools

import swiglal as _swiglal
import numpy as _numpy
from numpy import float32,float64,complex64,complex128

import pycbc as _pycbc
import pycbc.scheme as _scheme

if _pycbc.HAVE_CUDA:
    import pycuda as _pycuda
    import pycuda.driver as _cudriver
    import pycuda.gpuarray as _cudaarray    
    from array_cuda import squared_norm as cuda_squared_norm
    
if _pycbc.HAVE_OPENCL:
    import pyopencl as _pyopencl
    import pyopencl.array as _openclarray

_ALLOWED_DTYPES = [_numpy.float32,_numpy.float64,_numpy.complex64,
                   _numpy.complex128]
_ALLOWED_SCALARS = [int,long, float, complex]+_ALLOWED_DTYPES


def _convert_to_scheme(ary):
    if ary._scheme is not _scheme.mgr.state:
        converted_array = Array(ary,dtype=ary._data.dtype)
        ary._data = converted_array._data
        ary._scheme = _scheme.mgr.state
        
def _convert(fn):
    # Convert this array to the current processing scheme
    @_functools.wraps(fn)
    def converted(self,*args):
        _convert_to_scheme(self)
        return fn(self,*args)
    return converted

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
            if self._data.dtype == float32 or self._data.dtype == float64:
                self.kind = 'real'
            else:
                self.kind = 'complex'
            if self._data.dtype == float32 or self._data.dtype == complex64:
                self.precision = 'single'
            else:
                self.precision = 'double'


        if copy:
            #Check that a valid dtype was given.
            if dtype is not None:
                if dtype not in _ALLOWED_DTYPES:
                    raise TypeError(str(dtype) + ' is not supported')
            elif hasattr(initial_array,'dtype'):
                dtype = initial_array.dtype
        
            #Unwrap initial_array
            input_data = None
            if isinstance(initial_array,Array):
                input_data = initial_array._data
            else:
                input_data = initial_array

            #If no dtype was given, default to Double, and Double Complex.
            if dtype is None:
                if type(input_data[0]) == complex:
                    dtype = complex128
                else:
                    dtype = float64

            if dtype == float32 or dtype == float64:
                self.kind = 'real'
            else:
                self.kind = 'complex'
            if dtype == float32 or dtype == complex64:
                self.precision = 'single'
            else:
                self.precision = 'double'

            #Create new instance with input_data as initialization.
            if self._scheme is None:
                if _pycbc.HAVE_CUDA and (type(input_data) is
                                         _cudaarray.GPUArray):
                    self._data = input_data.get().astype(dtype)
                elif _pycbc.HAVE_OPENCL and (type(input_data) is
                                             _openclarray.Array):
                    self._data = input_data.get().astype(dtype)
                else:
                    self._data = _numpy.array(input_data,dtype=dtype)
            elif _pycbc.HAVE_CUDA and (type(self._scheme) is
                                       _scheme.CUDAScheme):
                if type(input_data) is _cudaarray.GPUArray:
                    if input_data.dtype == dtype:
                        self._data = _cudaarray.GPUArray((input_data.size),dtype)
                        _cudriver.memcpy_dtod(self._data.gpudata,input_data.gpudata,
                                                input_data.nbytes)
                    else:
                        self._data = input_data.astype(dtype)
                else:
                    self._data = _cudaarray.to_gpu(_numpy.array(input_data,
                                                            dtype=dtype))
            elif _pycbc.HAVE_OPENCL and (type(self._scheme) is
                                         _scheme.OpenCLScheme):
                if type(input_data) is _openclarray.Array:
                    if input_data.dtype == dtype:
                    #Unsure how to copy memory directly with OpenCL, these two lines are 
                    #a temporary workaround, still probably better than going to the host and back though
                        self._data = _openclarray.zeros_like(input_data)
                        self._data += input_data
                    else:
                        self._data = input_data.astype(dtype)
                else:
                    self._data = _openclarray.to_device(self._scheme.queue,
                                                    _numpy.array(input_data,
                                                                 dtype=dtype))
            else:
                raise TypeError('Invalid Processing scheme Type')

    def _returntype(fn):
        """ Wrapper for method functions to return a PyCBC Array class """
        @_functools.wraps(fn)
        def wrapped(self,*args):
            ary = fn(self,*args)
            if ary is NotImplemented:
                return NotImplemented
            return self._return(ary)
        return wrapped

    def _return(self,ary):
        """Wrap the ary to return an Array type """
        return Array(ary, copy=False)

    def _checkother(fn):
        """ Checks the input to method functions """
        @_functools.wraps(fn)
        def checked(self,other):
            self._typecheck(other) 
    
            if type(other) in _ALLOWED_SCALARS:
                pass
            elif isinstance(other,Array):
                if other.precision == self.precision:
                    _convert_to_scheme(other)
                    other = other._data
                else:
                    raise TypeError('dtypes do not match')
            else:
                return NotImplemented

            return fn(self,other)
        return checked

    def _typecheck(self,other):
        """ Additional typechecking for other. Stops Array from avoiding
        checks of derived types
        """
        pass

    @_returntype
    @_convert
    @_checkother
    def __mul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        return self._data * other

    @_returntype
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

    @_returntype
    @_convert
    @_checkother
    def __add__(self,other):
        """ Add Array to Array or scalar and return an Array. """
        return self._data + other
        
    def fill(self,value):
        self._data.fill(value)

    @_returntype
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
    @_returntype
    def __div__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        return self._data / other

    @_returntype
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

    @_returntype
    @_convert
    @_checkother
    def __sub__(self,other):
        """ Subtract Array or scalar from Array and return an Array. """
        return self._data - other

    @_returntype
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

    @_returntype
    @_convert
    @_checkother
    def __pow__(self,other):
        """ Exponentiate Arrary by scalar """
        return self._data ** other

    @_returntype
    @_convert
    def __abs__(self):
        """ Return absolute value of Array """
        return abs(self._data)

    def __len__(self):
        """ Return length of Array """
        return len(self._data)

    def __str__(self):
        return str(self._data)

    @_returntype
    @_convert
    def real(self):
        """ Return real part of Array """
        return Array(self._data.real,copy=True)

    @_returntype
    @_convert
    def imag(self):
        """ Return imaginary part of Array """
        return Array(self._data.imag,copy=True)

    @_returntype
    @_convert
    def conj(self):
        """ Return complex conjugate of Array. """
        return self._data.conj()
        
    @_returntype
    @_convert
    def squared_norm(self):
        """ Return the squared norm of the array """
        if type(self._data) is _numpy.ndarray:
            return (self.data * self.data.conj()).real
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            tmp = zeros(len(self),dtype=real_same_precision_as(self))
            cuda_squared_norm[self.precision](self.data,tmp.data)
            return tmp
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            raise NotImplementedError         

    @_convert
    def sum(self):
        """ Return the sum of the the array. """
        if type(self._data) is _numpy.ndarray:
            if self.kind == 'real':
                return _numpy.sum(self._data,dtype=float64)
            else:
                return _numpy.sum(self._data,dtype=complex128)
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return _pycuda.gpuarray.sum(self._data).get().max()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return _pyopencl.array.sum(self._data).get().max()
           
    def max(self):
        """ Return the maximum value in the array. """
        if type(self._data) is _numpy.ndarray:
            return self._data.max()
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return _pycuda.gpuarray.sum(self._data).get().max()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return _pyopencl.array.sum(self._data).get().max()        
     
    @_convert
    def max(self):
        """ Return the maximum value in the array. """
        if type(self._data) is _numpy.ndarray:
            return self._data.max()
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return _pycuda.gpuarray.max(self._data).get().max()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return _pyopencl.array.max(self._data).get().max()  
            
    @_convert
    def min(self):
        """ Return the maximum value in the array. """
        if type(self._data) is _numpy.ndarray:
            return self._data.min()
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return _pycuda.gpuarray.min(self._data).get().max()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return _pyopencl.array.min(self._data).get().max()         

    @_convert
    @_checkother
    def dot(self,other):
        """ Return the dot product"""
        if type(self._data) is _numpy.ndarray:
            return _numpy.dot(self._data,other)
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return _pycuda.gpuarray.dot(self._data,other).get().max()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return _pyopencl.array.dot(self._data,other).get().max()

    @_convert
    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._return(self._data[index])
        else:
            if type(self._data) is _numpy.ndarray:
                return self._data[index]
            elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
                return self._data.get()[index]
            elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
                return self._data.get()[index]
                
    @_convert
    def __setitem__(self,index,other):
        if isinstance(other,Array):
            if type(self._data) is _numpy.ndarray:
                self._data[index]=other.data
            elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
                if (other._data[index].nbytes <= self._data[index].nbytes) :
                    _cudriver.memcpy_dtod_async(self[index].ptr,other.ptr,other._data.nbytes)
                else:
                    raise RuntimeError("The arrays must the same length")
            elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
                raise NotImplementedError
        else:
            raise TypeError('Can only copy data from another Array')
                

    @property
    @_convert
    def data(self):
        """Returns the internal python array """
        return self._data

    @data.setter
    def data(self,other):
        dtype = None
        if hasattr(other,'dtype'):
            dtype = other.dtype
        temp = Array(other,dtype=dtype)
        self._data = temp._data


    @property
    @_convert
    def ptr(self):
        """ Returns a pointer to the memory of this array """
        if type(self._data) is _numpy.ndarray:
            raise TypeError("Please use swiglal for CPU objects")
        if _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return self._data.ptr
        if _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return self._data.data

    @property
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

        lal_data.data[:] = self._data
        self._data = lal_data.data

        return lal_data

    @property
    def dtype(self):
        return self._data.dtype

# Convenience functions for determining dtypes            

def real_same_precision_as(data):
    if data.precision is 'single':
        return float32
    elif data.precision is 'double':
        return float64
        
def complex_same_precision_as(data):
    if data.precision is 'single':
        return complex64
    elif data.precision is 'double':
        return complex128            

def zeros(length,dtype=None):
    if _scheme.mgr.state is None:
        return Array(_numpy.zeros(length,dtype=dtype),copy=False)
    if type(_scheme.mgr.state) is _scheme.CUDAScheme:
        return Array(_cudaarray.zeros(length,dtype),copy=False)
    if type(_scheme.mgr.state) is _scheme.OpenCLScheme:
        return Array(_openclarray.zeros(_scheme.mgr.state.queue,length,dtype),copy=False)

