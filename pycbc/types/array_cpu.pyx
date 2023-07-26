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
from pycbc.types.array import common_kind, complex128, float64
from . import aligned as _algn
from scipy.linalg import blas
from pycbc.types import real_same_precision_as
cimport cython, numpy

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
    if self.kind == 'real':
        tmp = abs(self.data)
        ind = _np.argmax(tmp)
        return tmp[ind], ind
    else:
        tmp = self.data.real ** 2.0
        tmp += self.data.imag ** 2.0
        ind = _np.argmax(tmp)
        return tmp[ind] ** 0.5, ind

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

ctypedef fused REALTYPE:
    float
    double

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_real(numpy.ndarray [REALTYPE, ndim = 1] a, numpy.ndarray [REALTYPE, ndim = 1] b):
    cdef double total = 0
    cdef unsigned int xmax = a.shape[0]
    cdef unsigned int i

    cdef REALTYPE* x = &a[0]
    cdef REALTYPE* y = &b[0]

    for i in range(xmax):
        total += x[i] * y[i]
    return total

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

def abs_arg_max_complex(numpy.ndarray [COMPLEXTYPE, ndim=1] a):
    cdef unsigned int xmax = a.shape[0]
    cdef double mag
    cdef double magmax = 0
    cdef unsigned int idx = 0

    for i in range(xmax):
        mag = a[i].real * a[i].real + a[i].imag * a[i].imag
        if mag > magmax:
            magmax = mag
            idx = i

    return idx

def abs_arg_max(self):
    if self.dtype == _np.float32 or self.dtype == _np.float64:
        return _np.argmax(abs(self.data))
    else:
        return abs_arg_max_complex(self._data)

def inner(self, other):
    """ Return the inner product of the array with complex conjugation.
    """
    cdtype = common_kind(self.dtype, other.dtype)
    if cdtype.kind == 'c':
        return _np.sum(self.data.conj() * other, dtype=complex128)
    else:
        return inner_real(self.data, other)

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
    if isinstance(array, _np.ndarray):
        return True
    else:
        return False
