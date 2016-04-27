# Copyright (C) 2012  Alex Nitz
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
"""Numpy based CPU backend for PyCBC Array
"""
import numpy as _np
from array import common_kind, complex128, float64
import aligned as _algn
from scipy.linalg import blas

def zeros(length, dtype=_np.float64):
    return _algn.zeros(length, dtype=dtype)

def empty(length, dtype=_np.float64):
    return _algn.empty(length, dtype=dtype)

def ptr(self):
    return self.data.ctypes.data
    
def dot(self, other):
    return _np.dot(self._data,other)     

def min(self):
    return self.data.min()  

def abs_max_loc(self):
    tmp = abs(self.data)
    ind = _np.argmax(tmp)
    return tmp[ind], ind

def cumsum(self):
    return self.data.cumsum()

def max(self):
    return self.data.max()

def max_loc(self):
    ind = _np.argmax(self.data)
    return self.data[ind], ind
    
def take(self, indices):
    return self.data.take(indices)
    
def weighted_inner(self, other, weight):
    """ Return the inner product of the array with complex conjugation.
    """
    if weight is None:
        return self.inner(other)
        
    cdtype = common_kind(self.dtype, other.dtype)
    if cdtype.kind == 'c':
        acum_dtype = complex128
    else:
        acum_dtype = float64

    return _np.sum(self.data.conj() * other / weight, dtype=acum_dtype)
    
def inner(self, other):
    """ Return the inner product of the array with complex conjugation.
    """
    cdtype = common_kind(self.dtype, other.dtype)
    if cdtype.kind == 'c':
        acum_dtype = complex128
    else:
        acum_dtype = float64
    return _np.sum(self.data.conj() * other, dtype=acum_dtype)
    #return _np.vdot(self.data, other)

def vdot(self, other):
    """ Return the inner product of the array with complex conjugation.
    """
    return _np.vdot(self.data, other)

def squared_norm(self):
    """ Return the elementwise squared norm of the array """
    return (self.data.real**2 + self.data.imag**2)

_blas_mandadd_funcs = {}
_blas_mandadd_funcs[_np.float32] = blas.saxpy
_blas_mandadd_funcs[_np.float64] = blas.daxpy
_blas_mandadd_funcs[_np.complex64] = blas.caxpy
_blas_mandadd_funcs[_np.complex128] = blas.zaxpy

def multiply_and_add(self, other, mult_fac):
    """
    Return other multiplied by mult_fac and with self added.
    Self will be modified in place. This requires all inputs to be of the same
    precision.
    """
    # Sanity checking should have already be done. But we don't know if
    # mult_fac and add_fac are arrays or scalars.
    inpt = _np.array(self.data, copy=False)
    N = len(inpt)
    # For some reason, _checkother decorator returns other.data so we don't
    # take .data here
    other = _np.array(other, copy=False)

    assert(inpt.dtype == other.dtype)

    blas_fnc = _blas_mandadd_funcs[inpt.dtype.type]
    return blas_fnc(other, inpt, a=mult_fac)

def numpy(self):
    return self._data
    
def _copy(self, self_ref, other_ref):
    self_ref[:] = other_ref[:]
    
def _getvalue(self, index):
    return self._data[index]

def sum(self):
    if self.kind == 'real':
        return _np.sum(self._data,dtype=float64)
    else:
        return _np.sum(self._data,dtype=complex128)   
        
def clear(self): 
    self[:] = 0 
    
def _scheme_matches_base_array(array):
    # Since ArrayWithAligned is a subclass of ndarray,
    # and since converting to ArrayWithAligned will
    # *not* copy 'array', the following is the way to go:
    if isinstance(array, _np.ndarray):
        return True
    else:
        return False
