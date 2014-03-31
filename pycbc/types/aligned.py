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
This module provides a clas derived from numpy.ndarray that also indicates
whether or not its memory is aligned.  It further provides functions for
creating zeros and empty (unitialized) arrays with this class.
"""
import numpy as _np
from pycbc import PYCBC_ALIGNMENT

def check_aligned(ndarr):
    return ((ndarr.__array_interface__['data'][0] % PYCBC_ALIGNMENT) == 0)

class ArrayWithAligned(_np.ndarray):
    def __new__(cls, input_array):
        obj = _np.asarray(input_array).view(cls)
        obj.isaligned = check_aligned(input_array)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.isaligned = check_aligned(self)

def zeros(n, dtype):
    d = _np.dtype(dtype)
    nbytes = (d.itemsize)*n
    #print "nbytes = {0}".format(nbytes)                                                                      
    tmp = _np.zeros(nbytes+PYCBC_ALIGNMENT, dtype=_np.uint8)
    address = tmp.__array_interface__['data'][0]
    offset = (PYCBC_ALIGNMENT - address%PYCBC_ALIGNMENT)%PYCBC_ALIGNMENT
    return ArrayWithAligned(tmp[offset:offset+nbytes].view(dtype=_np.dtype(dtype)))

def empty(n, dtype):
    d = _np.dtype(dtype)
    nbytes = (d.itemsize)*n
    tmp = _np.empty(nbytes+PYCBC_ALIGNMENT, dtype=_np.uint8)
    address = tmp.__array_interface__['data'][0]
    offset = (PYCBC_ALIGNMENT - address%PYCBC_ALIGNMENT)%PYCBC_ALIGNMENT
    return ArrayWithAligned(tmp[offset:offset+nbytes].view(dtype=d))
