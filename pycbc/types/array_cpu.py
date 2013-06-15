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
import numpy
from array import common_kind, complex128, float64

def zeros(length, dtype=numpy.float64):
    return numpy.zeros(length, dtype=dtype)

def ptr(self):
    raise TypeError("Please use lal for CPU objects")

def dot(self, other):
    return numpy.dot(self._data,other)     

def min(self):
    return self.data.min()  

def abs_max_loc(self):
    tmp = abs(self.data)
    ind = numpy.argmax(tmp)
    return tmp[ind], ind

def cumsum(self):
    return self.data.cumsum()

def max(self):
    return self.data.max()

def max_loc(self):
    return self.data.max(), numpy.argmax(self._data)
    
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

    return numpy.sum(self.data.conj() * other / weight, dtype=acum_dtype)

    
