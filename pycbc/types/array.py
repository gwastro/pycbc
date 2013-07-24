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

BACKEND_PREFIX="pycbc.types.array_"

import functools as _functools
from decorator import decorator

import lal as _lal
import numpy as _numpy
from numpy import float32, float64, complex64, complex128

import pycbc as _pycbc
import pycbc.scheme as _scheme
from pycbc.scheme import schemed, cpuonly

if _pycbc.HAVE_CUDA:
    import pycuda as _pycuda
    import pycuda.driver as _cudriver
    import pycuda.gpuarray as _cudaarray    
    
if _pycbc.HAVE_OPENCL:
    import pyopencl as _pyopencl
    import pyopencl.array as _openclarray

_ALLOWED_DTYPES = [_numpy.float32,_numpy.float64,_numpy.complex64,
                   _numpy.complex128]
_ALLOWED_SCALARS = [int,long, float, complex]+_ALLOWED_DTYPES

def _convert_to_scheme(ary):
    if ary._scheme is not _scheme.mgr.state:
        converted_array = Array(ary, dtype=ary._data.dtype)
        ary._data = converted_array._data
        ary._scheme = _scheme.mgr.state
        
def _convert(fn):
    # Convert this array to the current processing scheme
    @_functools.wraps(fn)
    def converted(self,*args):
        _convert_to_scheme(self)
        return fn(self,*args)
    return converted

def force_precision_to_match(scalar, precision):
    if _numpy.iscomplex(scalar):
        if precision is 'single':
            return _numpy.complex64(scalar)
        else:
            return _numpy.complex128(scalar)
    else:
        if precision is 'single':
            return _numpy.float32(scalar)
        else:
            return _numpy.float64(scalar)
        

def common_kind(*dtypes):
    for dtype in dtypes:
        if dtype.kind is 'c':
            return dtype
    return dtypes[0]

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
            if isinstance(initial_array, Array):
                if self._scheme is not initial_array._scheme:
                    raise TypeError("Cannot avoid a copy of this Array")
                self._data = initial_array._data
            elif type(initial_array) is _numpy.ndarray:
                if type(self._scheme) is not _scheme.CPUScheme:
                    raise TypeError("Cannot avoid a copy of this Array")
                else:
                    self._data = initial_array
            elif _pycbc.HAVE_CUDA and type(initial_array) is _cudaarray.GPUArray:
                if type(self._scheme) is not _scheme.CUDAScheme:
                    raise TypeError("Cannot avoid a copy of this Array")
                else:
                    self._data = initial_array
            elif _pycbc.HAVE_OPENCL and type(initial_array) is _openclarray.Array:
                if type(self._scheme) is not _scheme.OpenCLScheme:
                    raise TypeError("Cannot avoid a copy of this Array")
                else:
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

            if dtype and dtype != self._data.dtype:
                raise TypeError("Cannot set dtype when not copying")


        if copy:
            # First we will check the dtype that we are given
            initdtype = None
            if not hasattr(initial_array,'dtype'):
                initial_array = _numpy.array(initial_array)

            if initial_array.dtype in _ALLOWED_DTYPES:
                initdtype = initial_array.dtype
            else:
                if initial_array.dtype.kind == 'c':
                    initdtype = complex128
                else:
                    initdtype = float64
                        
            # Now that we know the dtype of the data, we can determine whether the specified dtype
            # is valid. If the data is complex, and a real dtype has been specified, this should
            # raise an error.
            if dtype is not None:
                if dtype not in _ALLOWED_DTYPES:
                    raise TypeError(str(dtype) + ' is not supported')
                if _numpy.dtype(dtype).kind != 'c' and _numpy.dtype(initdtype).kind == 'c':
                    raise TypeError(str(initdtype) + ' cannot be cast as ' + str(dtype))
            else:
                dtype = initdtype  
                       
            if dtype == float32 or dtype == float64:
                self.kind = 'real'
            else:
                self.kind = 'complex'
            if dtype == float32 or dtype == complex64:
                self.precision = 'single'
            else:
                self.precision = 'double'                        
            #Unwrap initial_array
            input_data = None
            if isinstance(initial_array,Array):
                input_data = initial_array._data
            else:
                input_data = initial_array

            #Create new instance with input_data as initialization.
            if type(self._scheme) is _scheme.CPUScheme:
                if _pycbc.HAVE_CUDA and (type(input_data) is
                                         _cudaarray.GPUArray):
                    self._data = input_data.get().astype(dtype)
                elif _pycbc.HAVE_OPENCL and (type(input_data) is
                                             _openclarray.Array):
                    self._data = input_data.get().astype(dtype)
                else:
                    self._data = _numpy.array(input_data,dtype=dtype,ndmin=1)
            elif _pycbc.HAVE_CUDA and (type(self._scheme) is
                                       _scheme.CUDAScheme):
                if type(input_data) is _cudaarray.GPUArray:
                    if input_data.dtype == dtype:
                        self._data = _cudaarray.GPUArray((input_data.size),dtype)
                        if len(input_data) > 0:
                            _cudriver.memcpy_dtod(self._data.gpudata,input_data.gpudata,
                                                    input_data.nbytes)
                    else:
                        self._data = input_data.astype(dtype)
                else:
                    self._data = _cudaarray.to_gpu(_numpy.array(input_data,
                                                            dtype=dtype,ndmin=1))
            elif _pycbc.HAVE_OPENCL and (type(self._scheme) is
                                         _scheme.OpenCLScheme):
                if type(input_data) is _openclarray.Array:
                    if input_data.dtype == dtype:
                    #Unsure how to copy memory directly with OpenCL, these two lines are 
                    #a temporary workaround, still probably better than going to the host and back though
                        # This function doesn't behave nicely when the array is empty
                        # because it tries to fill it with zeros
                        if len(input_data) > 0:
                            self._data = _openclarray.zeros_like(input_data)
                            self._data += input_data
                        else:
                            self._data = _openclarray.Array(self._scheme.queue,(0,),dtype)
                    else:
                        self._data = input_data.astype(dtype)
                else:
                    self._data = _openclarray.to_device(self._scheme.queue,
                                                    _numpy.array(input_data,
                                                                 dtype=dtype,ndmin=1))
            else:
                raise TypeError('Invalid Processing scheme Type')

    @decorator
    def _returnarray(fn, self, *args):
        return Array(fn(self, *args), copy=False)

    @decorator
    def _returntype(fn, self, *args):
        ary = fn(self,*args)
        if ary is NotImplemented:
            return NotImplemented
        return self._return(ary)
        
    def _return(self,ary):
        """Wrap the ary to return an Array type """
        return Array(ary, copy=False)

    @decorator
    def _checkother(fn, self,*args):
        nargs = ()
        for other in args:
            self._typecheck(other)  
            if type(other) in _ALLOWED_SCALARS:
                other = force_precision_to_match(other, self.precision)
                nargs +=(other,)
            elif isinstance(other, type(self)) or type(other) is Array:
                if len(other) != len(self):
                    raise ValueError('lengths do not match')
                if other.precision == self.precision:
                    _convert_to_scheme(other)
                    nargs += (other._data,)
                else:
                    raise TypeError('precisions do not match')
            else:
                return NotImplemented

        return fn(self,*nargs)
    
    @decorator  
    def _vcheckother(fn, self,*args):
        nargs = ()
        for other in args:
            self._typecheck(other)  
            if isinstance(other, type(self)) or type(other) is Array:
                if len(other) != len(self):
                    raise ValueError('lengths do not match')
                if other.precision == self.precision:
                    _convert_to_scheme(other)
                    nargs += (other._data,)
                else:
                    raise TypeError('precisions do not match')
            else:
                raise TypeError('array argument required')                    

        return fn(self,*nargs)

        
    def _icheckother(fn):
        """ Checks the input to in-place operations """
        @_functools.wraps(fn)
        def checked(self,other):
            self._typecheck(other) 
    
            if type(other) in _ALLOWED_SCALARS:
                if self.kind == 'real' and type(other) == complex:
                    raise TypeError('dtypes are incompatible')
                other = force_precision_to_match(other, self.precision)
            elif isinstance(other, type(self)) or type(other) is Array:
                if len(other) != len(self):
                    raise ValueError('lengths do not match')
                if self.kind == 'real' and other.kind == 'complex':
                    raise TypeError('dtypes are incompatible')
                if other.precision == self.precision:
                    _convert_to_scheme(other)
                    other = other._data
                else:
                    raise TypeError('precisions do not match')
            else:
                return NotImplemented

            return fn(self,other)
        return checked

    def _typecheck(self,other):
        """ Additional typechecking for other. Placeholder for use by derived
        types. 
        """
        pass

    @_returntype
    @_convert
    @_checkother
    def __mul__(self,other):
        """ Multiply by an Array or a scalar and return an Array. """
        return self._data * other

    __rmul__ = __mul__

    @_convert
    @_icheckother
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

    __radd__ = __add__
       
    def fill(self,value):
        self._data.fill(value)

    @_convert
    @_icheckother
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
        return self._data.__rdiv__(other)

    @_convert
    @_icheckother
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
        return self._data.__rsub__(other)

    @_convert
    @_icheckother
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

    def __eq__(self,other):
        """
        This is the Python special method invoked whenever the '=='
        comparison is used.  We want it to return true if the data
        of our PyCBC array-like type are identical, and all of the
        numeric meta-data are identical, irrespective of whether or
        not the two instances live in the same memory (for that
        comparison, the Python statement 'a is b' should be used
        instead).

        Thus, this method returns 'True' if the types of both 'self'
        and 'other' are identical, as well as their dtypes and the
        data in the arrays, element by element.  It will always do
        the comparison on the CPU, but will *not* move either object
        to the CPU if it is not already there.  This may make it
        slower than best-possible when both objects live on the GPU.
        If this is an issue the implementation of this method may be
        revisited later, but for now we do not see many use cases
        beyond unit-testing, where speed is not such a concern.
        Element-by-element comparisons will always be slow when done
        on a GPU, even in a reduction kernel.

        Note in particular that this function returns a single boolean,
        and not an array of booleans as Numpy does.  If the numpy
        behavior is desired it can be obtained using the numpy() method
        of the PyCBC type.

        Parameters
        ----------
        other: another Python object, that should be tested for
            equality with 'self'.

        Returns
        -------
        boolean: 'True' if the types, data, and meta-data (except
            scheme) of the two objects are each identical.
        """

        # Writing the first test as below allows this method to be safely
        # called from subclasses.
        if type(self) != type(other):
            return False
        if self.dtype != other.dtype:
            return False
        if len(self) != len(other):
            return False

        # Now we've checked meta-data, so look at the actual data itself:
        # The numpy() method call will put a copy of GPU data onto a CPU
        # array, and could therefore be slow.  As noted in the help for
        # this function we don't worry about that.

        sary = self.numpy()
        oary = other.numpy()

        # Now we know that both sary and oary are numpy arrays. The
        # '==' statement returns an array of booleans, and the all()
        # method of that array returns 'True' only if every element
        # of that array of booleans is True.
        return bool((sary == oary).all())

    @_returntype
    @_convert
    def real(self):
        """ Return real part of Array """
        return Array(self._data.real, copy=True)

    @_returntype
    @_convert
    def imag(self):
        """ Return imaginary part of Array """
        return Array(self._data.imag, copy=True)

    @_returntype
    @_convert
    def conj(self):
        """ Return complex conjugate of Array. """
        return self._data.conj()
        
    @_returntype
    @_convert
    @schemed(BACKEND_PREFIX)
    def squared_norm(self):
        """ Return the elementwise squared norm of the array """

    @_vcheckother
    @_convert
    @schemed(BACKEND_PREFIX)
    def inner(self, other):
        """ Return the inner product of the array with complex conjugation.
        """

    @_convert
    def clear(self):
        """ Clear out the values of the array. """
        if type(self._scheme) is _scheme.CPUScheme:
            self[:] = 0 
        if type(self._scheme) is _scheme.CUDAScheme:
            n32 = self.data.nbytes / 4
            _cudriver.memset_d32(self.data.gpudata, 0, n32)
        if (self._scheme) is _scheme.OpenCLScheme:
            self[:] = 0

    @_vcheckother
    @_convert
    @schemed(BACKEND_PREFIX)
    def weighted_inner(self, other, weight):
        """ Return the inner product of the array with complex conjugation.
        """

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

    @_returntype
    @_convert
    @schemed(BACKEND_PREFIX)
    def cumsum(self):
        """ Return the cumulative sum of the the array. """
     
    @_convert
    @schemed(BACKEND_PREFIX)
    def max(self):
        """ Return the maximum value in the array. """
            
    @_convert
    @schemed(BACKEND_PREFIX)
    def max_loc(self):
        """Return the maximum value in the array along with the index location """

    @_convert
    @schemed(BACKEND_PREFIX)
    def abs_max_loc(self):
        """Return the maximum elementwise norm in the array along with the index location"""

    @_convert
    @schemed(BACKEND_PREFIX)
    def min(self):
        """ Return the maximum value in the array. """ 
        
    @_returnarray
    @_convert
    @schemed(BACKEND_PREFIX)
    def take(self, indices):
        """ Return the values at the given indices. """                           

    @_convert
    @_vcheckother
    @schemed(BACKEND_PREFIX)
    def dot(self, other):
        """ Return the dot product"""
            
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
    def resize(self, new_size):
        """Resize self to new_size
        """
        if new_size == len(self):
            return
        else:
            new_arr = zeros(new_size, dtype=self.dtype)
            if len(self) <= new_size:
                new_arr[0:len(self)] = self
            else:
                new_arr[:] = self[0:new_size]
                
            self._data = new_arr._data

    @_convert
    def roll(self, shift):
        """shift vector
        """
        new_arr = zeros(len(self), dtype=self.dtype)

        if shift == 0:
            return
        if shift < 0:
            shift=len(self) + shift

        new_arr[0:shift] = self[len(self)-shift: len(self)]
        new_arr[shift:len(self)] = self[0:len(self)-shift]
            
        self._data = new_arr._data

    @_returntype
    @_convert
    def astype(self, dtype):
        if self.dtype is dtype:
            return self
        else:
            return self._data.astype(dtype)
                
    @_convert
    def __setitem__(self, index, other):
        if isinstance(other,Array):
            _convert_to_scheme(other)

            if self.kind is 'real' and other.kind is 'complex':
                raise ValueError('Cannot set real value with complex')

            if isinstance(index,slice):          
                self_ref = self._data[index]
                other_ref = other._data
            else:
                self_ref = self._data[index:index+1]
                other_ref = other._data

            if type(self._data) is _numpy.ndarray:
                self_ref[:] = other_ref[:]

            elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
                if (len(other_ref) <= len(self_ref)) :
                    from pycuda.elementwise import get_copy_kernel
                    func = get_copy_kernel(self.dtype, other.dtype)
                    func.prepared_async_call(self_ref._grid, self_ref._block, None,
                            self_ref.gpudata, other_ref.gpudata,
                            self_ref.mem_size)
                else:
                    raise RuntimeError("The arrays must the same length")

            elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
                if (len(other_ref) <= len(self_ref)) :
                    self_ref._copy(self_ref, other_ref)
                else:
                    raise RuntimeError("The arrays must the same length")

        
        elif type(other) in _ALLOWED_SCALARS:
            if isinstance(index, slice):
                self[index].fill(other)
            else:
                self[index:index+1].fill(other)
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
    @schemed(BACKEND_PREFIX)
    def ptr(self):
        """ Returns a pointer to the memory of this array """        

    @property
    @cpuonly
    @_convert
    def  _swighelper(self):
        """ Used internally by SWIG typemaps to ensure @_convert 
            is called and scheme is correct  
        """
        return self;

    def  numpy(self):
        """ Returns a Numpy Array that contains this data """
        if type(self._data) is _numpy.ndarray:
            return self._data
        elif _pycbc.HAVE_CUDA and type(self._data) is _cudaarray.GPUArray:
            return self._data.get()
        elif _pycbc.HAVE_OPENCL and type(self._data) is _openclarray.Array:
            return self._data.get()       
    
    @cpuonly
    @_convert
    def  lal(self):
        """ Returns a LAL Object that contains this data """

        lal_data = None
        if self._data.dtype == float32:
            lal_data = _lal.CreateREAL4Vector(len(self))
        elif self._data.dtype == float64:
            lal_data = _lal.CreateREAL8Vector(len(self))
        elif self._data.dtype == complex64:
            lal_data = _lal.CreateCOMPLEX8Vector(len(self))
        elif self._data.dtype == complex128:
            lal_data = _lal.CreateCOMPLEX16Vector(len(self))

        lal_data.data[:] = self._data

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

@decorator
def _return_array(fn, *args, **kwds):
    return Array(fn(*args, **kwds), copy=False)


@_return_array
@schemed(BACKEND_PREFIX)
def zeros(length, dtype=float64):
    """ Return an Array filled with zeros.
    """
    pass



