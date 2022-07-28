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


def create_memory_and_engine_for_class_based_fft(
    npoints_time,
    dtype,
    delta_t=1,
    ifft=False
):
    """ Create memory and engine for class-based FFT/IFFT

    Currently only supports R2C FFT / C2R IFFTs, but this could be expanded
    if use-cases arise.

    Parameters
    ----------
    npoints_time : int
        Number of time samples of the real input vector (or real output vector
        if doing an IFFT).
    dtype : np.dtype
        The dtype for the real input vector (or real output vector if doing an
        IFFT). np.float32 or np.float64 I think in all cases.
    delta_t : float (default 1)
        delta_t of the real vector. If not given this will be set to 1, and we
        will assume it is not needed in the returned TimeSeries/FrequencySeries
    ifft : boolean (default False)
        By default will use the FFT class, set to true to use IFFT.
    """
    from pycbc.types import FrequencySeries, TimeSeries, zeros
    from pycbc.types import complex_same_precision_as

    npoints_freq = npoints_time // 2 + 1
    delta_f_tmp = 1.0 / (npoints_time * delta_t)
    vec = TimeSeries(
        zeros(
            npoints_time,
            dtype=dtype
        ),
        delta_t=delta_t,
        copy=False
    )
    vectilde = FrequencySeries(
        zeros(
            npoints_freq,
            dtype=complex_same_precision_as(vec)
        ),
        delta_f=delta_f_tmp,
        copy=False
    )
    if ifft:
        fft_class = pycbc.fft.IFFT(vectilde, vec)
        invec = vectilde
        outvec = vec
    else:
        fft_class = pycbc.fft.FFT(vec, vectilde)
        invec = vec
        outvec = vectilde

    return invec, outvec, fft_class

    
