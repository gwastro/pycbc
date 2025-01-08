# Copyright (C) 2024 Y Ddraig Goch
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
"""Cupy based CPU backend for PyCBC Array
"""
import cupy as cp
from pycbc.types.array import common_kind, complex128, float64

def zeros(length, dtype=cp.float64):
    return cp.zeros(length, dtype=dtype)

def empty(length, dtype=cp.float64):
    return cp.empty(length, dtype=dtype)

def ptr(self):
    return self.data.data.mem.ptr

def dot(self, other):
    return cp.dot(self._data,other)

def min(self):
    return self.data.min()

def abs_max_loc(self):
    if self.kind == 'real':
        tmp = abs(self.data)
        ind = cp.argmax(tmp)
        return tmp[ind], ind
    else:
        tmp = self.data.real ** 2.0
        tmp += self.data.imag ** 2.0
        ind = cp.argmax(tmp)
        return tmp[ind] ** 0.5, ind

def cumsum(self):
    return self.data.cumsum()

def max(self):
    return self.data.max()

def max_loc(self):
    ind = cp.argmax(self.data)
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

    return cp.sum(self.data.conj() * other / weight, dtype=acum_dtype)

def abs_arg_max(self):
    if self.dtype == cp.float32 or self.dtype == cp.float64:
        return cp.argmax(abs(self.data))
    else:
        return abs_arg_max_complex(self._data)

def inner(self, other):
    """ Return the inner product of the array with complex conjugation.
    """
    cdtype = common_kind(self.dtype, other.dtype)
    if cdtype.kind == 'c':
        return cp.sum(self.data.conj() * other, dtype=complex128)
    else:
        return inner_real(self.data, other)

def vdot(self, other):
    """ Return the inner product of the array with complex conjugation.
    """
    return cp.vdot(self.data, other)

def squared_norm(self):
    """ Return the elementwise squared norm of the array """
    return (self.data.real**2 + self.data.imag**2)

def numpy(self):
    return cp.asnumpy(self.data)

def _copy(self, self_ref, other_ref):
    self_ref[:] = other_ref[:]

def _getvalue(self, index):
    return self._data[index]

def sum(self):
    if self.kind == 'real':
        return cp.sum(self._data,dtype=float64)
    else:
        return cp.sum(self._data,dtype=complex128)

def clear(self):
    self[:] = 0

def _scheme_matches_base_array(array):
    if isinstance(array, cp.ndarray):
        return True
    else:
        return False

def _to_device(array):
    return cp.asarray(array)

def numpy(self):
    return cp.asnumpy(self._data)

def _copy_base_array(array):
    return array.copy()

