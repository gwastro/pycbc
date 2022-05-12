# Copyright (C) 2012  Alex Nitz, Josh Willis, Andrew Miller, Tito Dal Canton
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
This modules provides a device independent Array class based on PyCUDA and Numpy.
"""

BACKEND_PREFIX="pycbc.types.array_"

import h5py
import os as _os

from functools import wraps

import lal as _lal
import numpy as _numpy
from numpy import float32, float64, complex64, complex128, ones
from numpy.linalg import norm

import pycbc.scheme as _scheme
from pycbc.scheme import schemed, cpuonly
from pycbc.opt import LimitedSizeDict

#! FIXME: the uint32 datatype has not been fully tested,
# we should restrict any functions that do not allow an
# array of uint32 integers
_ALLOWED_DTYPES = [_numpy.float32, _numpy.float64, _numpy.complex64,
                   _numpy.complex128, _numpy.uint32, _numpy.int32, int]
try:
    _ALLOWED_SCALARS = [int, long, float, complex] + _ALLOWED_DTYPES
except NameError:
    _ALLOWED_SCALARS = [int, float, complex] + _ALLOWED_DTYPES

def _convert_to_scheme(ary):
    if not isinstance(ary._scheme, _scheme.mgr.state.__class__):
        converted_array = Array(ary, dtype=ary._data.dtype)
        ary._data = converted_array._data
        ary._scheme = _scheme.mgr.state
      
def _convert(func):
    @wraps(func)
    def convert(self, *args, **kwargs):
        _convert_to_scheme(self)
        return func(self, *args, **kwargs)
    return convert
    
def _nocomplex(func):
    @wraps(func)
    def nocomplex(self, *args, **kwargs):
        if self.kind == 'real':
            return func(self, *args, **kwargs)
        else:
            raise TypeError( func.__name__ + " does not support complex types")
    return nocomplex

def _noreal(func):
    @wraps(func)
    def noreal(self, *args, **kwargs):
        if self.kind == 'complex':
            return func(self, *args, **kwargs)
        else:
            raise TypeError( func.__name__ + " does not support real types")
    return noreal

def force_precision_to_match(scalar, precision):
    if _numpy.iscomplexobj(scalar):
        if precision == 'single':
            return _numpy.complex64(scalar)
        else:
            return _numpy.complex128(scalar)
    else:
        if precision == 'single':
            return _numpy.float32(scalar)
        else:
            return _numpy.float64(scalar)

def common_kind(*dtypes):
    for dtype in dtypes:
        if dtype.kind == 'c':
            return dtype
    return dtypes[0]
   
@schemed(BACKEND_PREFIX) 
def _to_device(array):
    """ Move input to device """
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)
    
@schemed(BACKEND_PREFIX)
def _copy_base_array(array):
    """ Copy a backend array"""
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

@schemed(BACKEND_PREFIX)
def _scheme_matches_base_array(array):
    """ Check that input matches array type for scheme """
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

def check_same_len_precision(a, b):
    """Check that the two arguments have the same length and precision.
    Raises ValueError if they do not.
    """
    if len(a) != len(b):
        msg = 'lengths do not match ({} vs {})'.format(
                len(a), len(b))
        raise ValueError(msg)
    if a.precision != b.precision:
        msg = 'precisions do not match ({} vs {})'.format(
                a.precision, b.precision)
        raise TypeError(msg)

class Array(object):
    """Array used to do numeric calculations on a various compute
    devices. It is a convience wrapper around numpy, and
    pycuda.
    """

    def __init__(self, initial_array, dtype=None, copy=True):
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
        self._saved = LimitedSizeDict(size_limit=2**5)
        
        #Unwrap initial_array
        if isinstance(initial_array, Array):
            initial_array = initial_array._data

        if not copy:
            if not _scheme_matches_base_array(initial_array):
                raise TypeError("Cannot avoid a copy of this array")
            else:
                self._data = initial_array

            # Check that the dtype is supported.
            if self._data.dtype not in _ALLOWED_DTYPES:
                raise TypeError(str(self._data.dtype) + ' is not supported')

            if dtype and dtype != self._data.dtype:
                raise TypeError("Can only set dtype when allowed to copy data")


        if copy:
            # First we will check the dtype that we are given
            if not hasattr(initial_array, 'dtype'):
                initial_array = _numpy.array(initial_array)

            # Determine the dtype to use
            if dtype is not None:  
                dtype = _numpy.dtype(dtype)
                if dtype not in _ALLOWED_DTYPES:
                    raise TypeError(str(dtype) + ' is not supported')
                if dtype.kind != 'c' and initial_array.dtype.kind == 'c':
                    raise TypeError(str(initial_array.dtype) + ' cannot be cast as ' + str(dtype))          
            elif initial_array.dtype in _ALLOWED_DTYPES:
                dtype = initial_array.dtype
            else:
                if initial_array.dtype.kind == 'c':
                    dtype = complex128
                else:
                    dtype = float64
                     
            # Cast to the final dtype if needed
            if initial_array.dtype != dtype:
                initial_array = initial_array.astype(dtype)
                                              
            #Create new instance with initial_array as initialization.
            if issubclass(type(self._scheme), _scheme.CPUScheme):
                if hasattr(initial_array, 'get'):
                    self._data = _numpy.array(initial_array.get())
                else:
                    self._data = _numpy.array(initial_array, dtype=dtype, ndmin=1)
            elif _scheme_matches_base_array(initial_array):
                self._data = _copy_base_array(initial_array) # pylint:disable=assignment-from-no-return
            else:
                initial_array = _numpy.array(initial_array, dtype=dtype, ndmin=1)
                self._data = _to_device(initial_array) # pylint:disable=assignment-from-no-return

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        inputs = [i.numpy() if isinstance(i, Array) else i for i in inputs]
        ret = getattr(ufunc, method)(*inputs, **kwargs)
        if hasattr(ret, 'shape') and ret.shape == self.shape:
            ret = self._return(ret)
        return ret

    def __array__(self, dtype=None):
        arr = self.numpy()
        if dtype is not None:
            arr = arr.astype(dtype)
        return arr

    @property
    def shape(self):
        return self._data.shape
     
    def _memoize_single(func):
        @wraps(func)
        def memoize_single(self, arg):
            badh = str(arg)

            if badh in self._saved:
                return self._saved[badh]

            res = func(self, arg) # pylint:disable=not-callable
            self._saved[badh] = res
            return res
        return memoize_single

    def _returnarray(func):
        @wraps(func)
        def returnarray(self, *args, **kwargs):
            return Array(func(self, *args, **kwargs), copy=False) # pylint:disable=not-callable
        return returnarray

    def _returntype(func):
        @wraps(func)
        def returntype(self, *args, **kwargs):
            ary = func(self, *args, **kwargs) # pylint:disable=not-callable
            if ary is NotImplemented:
                return NotImplemented
            return self._return(ary)
        return returntype
        
    def _return(self, ary):
        """Wrap the ary to return an Array type """
        if isinstance(ary, Array):
            return ary
        return Array(ary, copy=False)

    def _checkother(func):
        @wraps(func)
        def checkother(self, *args):
            nargs = ()
            for other in args:
                self._typecheck(other)
                if type(other) in _ALLOWED_SCALARS:
                    other = force_precision_to_match(other, self.precision)
                    nargs +=(other,)
                elif isinstance(other, type(self)) or type(other) is Array:
                    check_same_len_precision(self, other)
                    _convert_to_scheme(other)
                    nargs += (other._data,)
                else:
                    return NotImplemented

            return func(self, *nargs) # pylint:disable=not-callable
        return checkother

    def _vcheckother(func):
        @wraps(func)
        def vcheckother(self, *args):
            nargs = ()
            for other in args:
                self._typecheck(other)
                if isinstance(other, type(self)) or type(other) is Array:
                    check_same_len_precision(self, other)
                    _convert_to_scheme(other)
                    nargs += (other._data,)
                else:
                    raise TypeError('array argument required')

            return func(self, *nargs) # pylint:disable=not-callable
        return vcheckother
        
    def _vrcheckother(func):
        @wraps(func)
        def vrcheckother(self, *args):
            nargs = ()
            for other in args:
                if isinstance(other, type(self)) or type(other) is Array:
                    check_same_len_precision(self, other)
                    _convert_to_scheme(other)
                    nargs += (other._data,)
                else:
                    raise TypeError('array argument required')

            return func(self, *nargs) # pylint:disable=not-callable
        return vrcheckother

    def _icheckother(func):
        @wraps(func)
        def icheckother(self, other):
            """ Checks the input to in-place operations """
            self._typecheck(other)
            if type(other) in _ALLOWED_SCALARS:
                if self.kind == 'real' and type(other) == complex:
                    raise TypeError('dtypes are incompatible')
                other = force_precision_to_match(other, self.precision)
            elif isinstance(other, type(self)) or type(other) is Array:
                check_same_len_precision(self, other)
                if self.kind == 'real' and other.kind == 'complex':
                    raise TypeError('dtypes are incompatible')
                _convert_to_scheme(other)
                other = other._data
            else:
                return NotImplemented

            return func(self, other) # pylint:disable=not-callable
        return icheckother

    def _typecheck(self, other):
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
       
    def fill(self, value):
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
    def __truediv__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        return self._data / other

    @_returntype
    @_convert
    @_checkother
    def __rtruediv__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        return self._data.__rtruediv__(other)

    @_convert
    @_icheckother
    def __itruediv__(self,other):
        """ Divide Array by Array or scalar and return an Array. """
        self._data /= other
        return self
        
    __div__ = __truediv__
    __idiv__ = __itruediv__
    __rdiv__ = __rtruediv__

    @_returntype
    @_convert
    def __neg__(self):
        """ Return negation of self """
        return - self._data

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
        """ Exponentiate Array by scalar """
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
        
    @property
    def ndim(self):
        return self._data.ndim

    def __eq__(self,other):
        """
        This is the Python special method invoked whenever the '=='
        comparison is used.  It will return true if the data of two
        PyCBC arrays are identical, and all of the numeric meta-data
        are identical, irrespective of whether or not the two
        instances live in the same memory (for that comparison, the
        Python statement 'a is b' should be used instead).

        Thus, this method returns 'True' if the types of both 'self'
        and 'other' are identical, as well as their lengths, dtypes
        and the data in the arrays, element by element.  It will always
        do the comparison on the CPU, but will *not* move either object
        to the CPU if it is not already there, nor change the scheme of
        either object. It is possible to compare a CPU object to a GPU
        object, and the comparison should be true if the data and
        meta-data of the two objects are the same.

        Note in particular that this function returns a single boolean,
        and not an array of booleans as Numpy does.  If the numpy
        behavior is instead desired it can be obtained using the numpy()
        method of the PyCBC type to get a numpy instance from each
        object, and invoking '==' on those two instances.

        Parameters
        ----------
        other: another Python object, that should be tested for equality
            with 'self'.

        Returns
        -------
        boolean: 'True' if the types, dtypes, lengths, and data of the
            two objects are each identical.
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
        return (sary == oary).all()

    def almost_equal_elem(self,other,tol,relative=True):
        """
        Compare whether two array types are almost equal, element
        by element.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
        abs(self[i]-other[i]) <= tol*abs(self[i])
        for all elements of the array.

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
        abs(self[i]-other[i]) <= tol
        for all elements of the array.

        Other meta-data (type, dtype, and length) must be exactly equal.
        If either object's memory lives on the GPU it will be copied to
        the CPU for the comparison, which may be slow.  But the original
        object itself will not have its memory relocated nor scheme
        changed.

        Parameters
        ----------
        other
            Another Python object, that should be tested for
            almost-equality with 'self', element-by-element.
        tol
            A non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative
            A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).

        Returns
        -------
        boolean 
            'True' if the data agree within the tolerance, as
            interpreted by the 'relative' keyword, and if the types,
            lengths, and dtypes are exactly the same.
        """
        # Check that the tolerance is non-negative and raise an
        # exception otherwise.
        if (tol<0):
            raise ValueError("Tolerance cannot be negative")
        # Check that the meta-data agree; the type check is written in
        # this way so that this method may be safely called from
        # subclasses as well.
        if type(other) != type(self):
            return False
        if self.dtype != other.dtype:
            return False
        if len(self) != len(other):
            return False

        # The numpy() method will move any GPU memory onto the CPU.
        # Slow, but the user was warned.

        diff = abs(self.numpy()-other.numpy())
        if relative:
            cmpary = tol*abs(self.numpy())
        else:
            cmpary = tol*ones(len(self),dtype=self.dtype)

        return (diff<=cmpary).all()

    def almost_equal_norm(self,other,tol,relative=True):
        """
        Compare whether two array types are almost equal, normwise.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
        abs(norm(self-other)) <= tol*abs(norm(self)).

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
        abs(norm(self-other)) <= tol

        Other meta-data (type, dtype, and length) must be exactly equal.
        If either object's memory lives on the GPU it will be copied to
        the CPU for the comparison, which may be slow.  But the original
        object itself will not have its memory relocated nor scheme
        changed.

        Parameters
        ----------
        other
            another Python object, that should be tested for
            almost-equality with 'self', based on their norms.
        tol 
            a non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative
            A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).

        Returns
        -------
        boolean
            'True' if the data agree within the tolerance, as
            interpreted by the 'relative' keyword, and if the types,
            lengths, and dtypes are exactly the same.
        """
        # Check that the tolerance is non-negative and raise an
        # exception otherwise.
        if (tol<0):
            raise ValueError("Tolerance cannot be negative")
        # Check that the meta-data agree; the type check is written in
        # this way so that this method may be safely called from
        # subclasses as well.
        if type(other) != type(self):
            return False
        if self.dtype != other.dtype:
            return False
        if len(self) != len(other):
            return False

        # The numpy() method will move any GPU memory onto the CPU.
        # Slow, but the user was warned.

        diff = self.numpy()-other.numpy()
        dnorm = norm(diff)
        if relative:
            return (dnorm <= tol*norm(self))
        else:
            return (dnorm <= tol)

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
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_returntype
    @_checkother
    @_convert
    @schemed(BACKEND_PREFIX)
    def multiply_and_add(self, other, mult_fac):
        """ Return other multiplied by mult_fac and with self added.
        Self is modified in place and returned as output.
        Precisions of inputs must match.
        """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_vrcheckother
    @_convert
    @schemed(BACKEND_PREFIX)
    def inner(self, other):
        """ Return the inner product of the array with complex conjugation.
        """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_vrcheckother
    @_convert
    @schemed(BACKEND_PREFIX)
    def vdot(self, other):
        """ Return the inner product of the array with complex conjugation.
        """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_convert
    @schemed(BACKEND_PREFIX)
    def clear(self): 
        """ Clear out the values of the array. """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_vrcheckother
    @_convert
    @schemed(BACKEND_PREFIX)
    def weighted_inner(self, other, weight):
        """ Return the inner product of the array with complex conjugation.
        """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_convert
    @schemed(BACKEND_PREFIX)
    def sum(self):
        """ Return the sum of the the array. """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_returntype
    @_convert
    @schemed(BACKEND_PREFIX)
    def cumsum(self):
        """ Return the cumulative sum of the the array. """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)
     
    @_convert
    @_nocomplex
    @schemed(BACKEND_PREFIX)
    def max(self):
        """ Return the maximum value in the array. """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)
            
    @_convert
    @_nocomplex
    @schemed(BACKEND_PREFIX)
    def max_loc(self):
        """Return the maximum value in the array along with the index location """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_convert
    @schemed(BACKEND_PREFIX)
    def abs_arg_max(self):
        """ Return location of the maximum argument max """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_convert
    @schemed(BACKEND_PREFIX)
    def abs_max_loc(self):
        """Return the maximum elementwise norm in the array along with the index location"""
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_convert
    @_nocomplex
    @schemed(BACKEND_PREFIX)
    def min(self):
        """ Return the maximum value in the array. """ 
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)
        
    @_returnarray
    @_convert
    @schemed(BACKEND_PREFIX)
    def take(self, indices):
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_convert
    @_vcheckother
    @schemed(BACKEND_PREFIX)
    def dot(self, other):
        """ Return the dot product"""
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)
    
    @schemed(BACKEND_PREFIX)
    def _getvalue(self, index):
        """Helper function to return a single value from an array. May be very
           slow if the memory is on a gpu.
        """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)

    @_memoize_single
    @_returntype
    def _getslice(self, index):
        return self._return(self._data[index])
    
    @_convert
    def __getitem__(self, index):
        """ Return items from the Array. This not guaranteed to be fast for
            returning single values. 
        """
        if isinstance(index, slice):
            return self._getslice(index)
        else:
            return self._getvalue(index)

    @_convert
    def resize(self, new_size):
        """Resize self to new_size
        """
        if new_size == len(self):
            return
        else:
            self._saved = LimitedSizeDict(size_limit=2**5)
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

        if shift < 0:
            shift = shift - len(self) * (shift // len(self))
        
        if shift == 0:
            return
        
        new_arr[0:shift] = self[len(self)-shift: len(self)]
        new_arr[shift:len(self)] = self[0:len(self)-shift]
        
        self._saved = LimitedSizeDict(size_limit=2**5)
        
        self._data = new_arr._data

    @_returntype
    @_convert
    def astype(self, dtype):
        if _numpy.dtype(self.dtype) == _numpy.dtype(dtype):
            return self
        else:
            return self._data.astype(dtype)
    
    @schemed(BACKEND_PREFIX)
    def _copy(self, self_ref, other_ref):
        """Helper function to copy between two arrays. The arrays references
           should be bare array types and not `Array` class instances. 
        """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)
                
    @_convert
    def __setitem__(self, index, other):
        if isinstance(other,Array):
            _convert_to_scheme(other)

            if self.kind == 'real' and other.kind == 'complex':
                raise ValueError('Cannot set real value with complex')

            if isinstance(index,slice):          
                self_ref = self._data[index]
                other_ref = other._data
            else:
                self_ref = self._data[index:index+1]
                other_ref = other._data

            self._copy(self_ref, other_ref)

        elif type(other) in _ALLOWED_SCALARS:
            if isinstance(index, slice):
                self[index].fill(other)
            else:
                self[index:index+1].fill(other)
        else:
            raise TypeError('Can only copy data from another Array')

    @property
    def precision(self):
        if self.dtype == float32 or self.dtype == complex64:
            return 'single'
        else:
            return 'double'        
                
    @property
    def kind(self):
        if self.dtype == float32 or self.dtype == float64:
            return 'real'
        elif self.dtype == complex64 or self.dtype == complex128:
            return 'complex'
        else:
            return 'unknown'

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
        temp = Array(other, dtype=dtype)
        self._data = temp._data

    @property
    @_convert
    @schemed(BACKEND_PREFIX)
    def ptr(self):
        """ Returns a pointer to the memory of this array """
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)
        
    @property
    def itemsize(self):
        return self.dtype.itemsize
    
    @property
    def nbytes(self):
        return len(self.data) * self.itemsize

    @property
    @cpuonly
    @_convert
    def _swighelper(self):
        """ Used internally by SWIG typemaps to ensure @_convert 
            is called and scheme is correct  
        """
        return self;

    @_convert
    @schemed(BACKEND_PREFIX)
    def numpy(self):
        """ Returns a Numpy Array that contains this data """     
        err_msg = "This function is a stub that should be overridden using "
        err_msg += "the scheme. You shouldn't be seeing this error!"
        raise ValueError(err_msg)
    
    @_convert
    def lal(self):
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

        lal_data.data[:] = self.numpy()

        return lal_data

    @property
    def dtype(self):
        return self._data.dtype
    
    def save(self, path, group=None):
        """
        Save array to a Numpy .npy, hdf, or text file. When saving a complex array as
        text, the real and imaginary parts are saved as the first and second
        column respectively. When using hdf format, the data is stored
        as a single vector, along with relevant attributes.

        Parameters
        ----------
        path: string
            Destination file path. Must end with either .hdf, .npy or .txt.
            
        group: string 
            Additional name for internal storage use. Ex. hdf storage uses
            this as the key value.

        Raises
        ------
        ValueError
            If path does not end in .npy or .txt.
        """

        ext = _os.path.splitext(path)[1]
        if ext == '.npy':
            _numpy.save(path, self.numpy())
        elif ext == '.txt':
            if self.kind == 'real':
                _numpy.savetxt(path, self.numpy())
            elif self.kind == 'complex':
                output = _numpy.vstack((self.numpy().real,
                                        self.numpy().imag)).T
                _numpy.savetxt(path, output)
        elif ext == '.hdf':
            key = 'data' if group is None else group
            with h5py.File(path, 'a') as f:
                f.create_dataset(key, data=self.numpy(), compression='gzip',
                                 compression_opts=9, shuffle=True)
        else:
            raise ValueError('Path must end with .npy, .txt, or .hdf')
           
    @_convert 
    def trim_zeros(self):
        """Remove the leading and trailing zeros.
        """      
        tmp = self.numpy()
        f = len(self)-len(_numpy.trim_zeros(tmp, trim='f'))
        b = len(self)-len(_numpy.trim_zeros(tmp, trim='b'))
        return self[f:len(self)-b]

    @_returntype
    @_convert
    def view(self, dtype):
        """
        Return a 'view' of the array with its bytes now interpreted according
        to 'dtype'. The location in memory is unchanged and changing elements
        in a view of an array will also change the original array.

        Parameters
        ----------
        dtype : numpy dtype (one of float32, float64, complex64 or complex128)
            The new dtype that should be used to interpret the bytes of self
        """
        return self._data.view(dtype)

    def copy(self):
        """ Return copy of this array """
        return self._return(self.data.copy())
        
    def __lt__(self, other):
        return self.numpy().__lt__(other)
        
    def __le__(self, other):
        return self.numpy().__le__(other)
        
    def __ne__(self, other):
        return self.numpy().__ne__(other)
        
    def __gt__(self, other):
        return self.numpy().__gt__(other)
        
    def __ge__(self, other):
        return self.numpy().__ge__(other)
            
# Convenience functions for determining dtypes
def real_same_precision_as(data):
    if data.precision == 'single':
        return float32
    elif data.precision == 'double':
        return float64

def complex_same_precision_as(data):
    if data.precision == 'single':
        return complex64
    elif data.precision == 'double':
        return complex128

def _return_array(func):
    @wraps(func)
    def return_array(*args, **kwds):
        return Array(func(*args, **kwds), copy=False)
    return return_array

@_return_array
@schemed(BACKEND_PREFIX)
def zeros(length, dtype=float64):
    """ Return an Array filled with zeros.
    """
    err_msg = "This function is a stub that should be overridden using "
    err_msg += "the scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

@_return_array
@schemed(BACKEND_PREFIX)
def empty(length, dtype=float64):
    """ Return an empty Array (no initialization)
    """
    err_msg = "This function is a stub that should be overridden using "
    err_msg += "the scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

def load_array(path, group=None):
    """Load an Array from an HDF5, ASCII or Numpy file. The file type is
    inferred from the file extension, which must be `.hdf`, `.txt` or `.npy`.

    For ASCII and Numpy files with a single column, a real array is returned.
    For files with two columns, the columns are assumed to contain the real
    and imaginary parts of a complex array respectively.

    The default data types will be double precision floating point.

    Parameters
    ----------
    path : string
        Input file path. Must end with either `.npy`, `.txt` or `.hdf`.

    group: string
        Additional name for internal storage use. When reading HDF files, this
        is the path to the HDF dataset to read.

    Raises
    ------
    ValueError
        If path does not end with a supported extension. For Numpy and ASCII
        input files, this is also raised if the array does not have 1 or 2
        dimensions.
    """
    ext = _os.path.splitext(path)[1]
    if ext == '.npy':
        data = _numpy.load(path)
    elif ext == '.txt':
        data = _numpy.loadtxt(path)
    elif ext == '.hdf':
        key = 'data' if group is None else group
        with h5py.File(path, 'r') as f:
            array = Array(f[key])
        return array
    else:
        raise ValueError('Path must end with .npy, .hdf, or .txt')

    if data.ndim == 1:
        return Array(data)
    elif data.ndim == 2:
        return Array(data[:,0] + 1j*data[:,1])

    raise ValueError('File has %s dimensions, cannot convert to Array, \
                      must be 1 (real) or 2 (complex)' % data.ndim)
