# Copyright (C) 2012  Josh Willis, Andrew Miller
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
This package provides a front-end to various fast Fourier transform
implementations within PyCBC.
"""

from .backend_support import get_backend

def _fft_factory(invec, outvec, nbatch=1, size=None):
    backend = get_backend()
    cls = getattr(backend, 'FFT')
    return cls

def _ifft_factory(invec, outvec, nbatch=1, size=None):
    backend = get_backend()
    cls = getattr(backend, 'IFFT')
    return cls

class FFT(object):
    """ Create a forward FFT  engine

    Parameters
    ----------
    invec : complex64 or float32
      Input pycbc.types.Array (or subclass); its FFT will be computed
    outvec : complex64 
      Output pycbc.types.Array (or subclass); it will hold the FFT of invec
    nbatch : int (default 1)
      When not one, specifies that invec and outvec should each be interpreted
      as nbatch distinct vectors. The total length of invec and outvec should
      then be that appropriate to a single vector, multiplied by nbatch
    size : int (default None)
      When nbatch is not 1, this parameter gives the logical size of each
      transform.  If nbatch is 1 (the default) this can be None, and the
      logical size is the length of invec.

    The addresses in memory of both vectors should be divisible by
    pycbc.PYCBC_ALIGNMENT.
    """
    def __new__(cls, *args, **kwargs):
        real_cls = _fft_factory(*args, **kwargs)
        return real_cls(*args, **kwargs)

class IFFT(object):
    """ Create a reverse FFT  engine

    Parameters
    ----------
    invec : complex64
      Input pycbc.types.Array (or subclass); its IFFT will be computed
    outvec : complex64 or float32 
      Output pycbc.types.Array (or subclass); it will hold the IFFT of invec
    nbatch : int (default 1)
      When not one, specifies that invec and outvec should each be interpreted
      as nbatch distinct vectors. The total length of invec and outvec should
      then be that appropriate to a single vector, multiplied by nbatch
    size : int (default None)
      When nbatch is not 1, this parameter gives the logical size of each
      transform.  If nbatch is 1 (the default) this can be None, and the
      logical size is the length of outvec.

    The addresses in memory of both vectors should be divisible by
    pycbc.PYCBC_ALIGNMENT.
    """
    def __new__(cls, *args, **kwargs):
        real_cls = _ifft_factory(*args, **kwargs)
        return real_cls(*args, **kwargs)

