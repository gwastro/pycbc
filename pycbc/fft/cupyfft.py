# Copyright (C) 2024 Y Ddraig Goch
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
This module provides the cupy backend of the fast Fourier transform
for the PyCBC package.
"""

from pycbc.types import Array
import logging
import cupy.fft
from .core import _check_fft_args
from .core import _BaseFFT, _BaseIFFT

_INV_FFT_MSG = ("I cannot perform an {} between data with an input type of "
                "{} and an output type of {}")

def fft(invec, outvec, _, itype, otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("cupy backend of pycbc.fft does not "
                                  "support in-place transforms")
    if itype == 'complex' and otype == 'complex':
        outvec.data[:] = cupy.asarray(cupy.fft.fft(invec.data),
                                       dtype=outvec.dtype)
    elif itype == 'real' and otype == 'complex':
        outvec.data[:] = cupy.asarray(cupy.fft.rfft(invec.data),
                                       dtype=outvec.dtype)
    else:
        raise ValueError(_INV_FFT_MSG.format("FFT", itype, otype))


def ifft(invec, outvec, _, itype, otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("cupy backend of pycbc.fft does not "
                                  "support in-place transforms")
    if itype == 'complex' and otype == 'complex':
        outvec.data[:] = cupy.asarray(cupy.fft.ifft(invec.data),
                                       dtype=outvec.dtype)
        outvec *= len(outvec)
    elif itype == 'complex' and otype == 'real':
        outvec.data[:] = cupy.asarray(cupy.fft.irfft(invec.data,len(outvec)),
                                       dtype=outvec.dtype)
        outvec *= len(outvec)
    else:
        raise ValueError(_INV_FFT_MSG.format("IFFT", itype, otype))


class FFT(_BaseFFT):
    """
    Class for performing FFTs via the cupy interface.
    """
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(FFT, self).__init__(invec, outvec, nbatch, size)
        self.prec, self.itype, self.otype = _check_fft_args(invec, outvec)

    def execute(self):
        fft(self.invec, self.outvec, self.prec, self.itype, self.otype)


class IFFT(_BaseIFFT):
    """
    Class for performing IFFTs via the cupy interface.
    """
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(IFFT, self).__init__(invec, outvec, nbatch, size)
        self.prec, self.itype, self.otype = _check_fft_args(invec, outvec)

    def execute(self):
        ifft(self.invec, self.outvec, self.prec, self.itype, self.otype)

def batch_fft(invecs, outvecs, _, itype, otype):
    """Batched FFT operation for multiple templates"""
    if itype == 'complex' and otype == 'complex':
        # Reshape input for batch operation
        batch_size = len(invecs)
        fft_size = len(invecs[0])
        batch_data = cupy.stack([v.data for v in invecs])
        
        # Perform batched FFT
        result = cupy.fft.fft(batch_data)
        
        # Copy results back to output vectors
        for i, outvec in enumerate(outvecs):
            outvec.data[:] = result[i]
            
    elif itype == 'real' and otype == 'complex':
        batch_size = len(invecs)
        fft_size = len(invecs[0])
        batch_data = cupy.stack([v.data for v in invecs])
        
        result = cupy.fft.rfft(batch_data)
        
        for i, outvec in enumerate(outvecs):
            outvec.data[:] = result[i]
    else:
        raise ValueError(_INV_FFT_MSG.format("FFT", itype, otype))

def batch_ifft(invecs, outvecs, _, itype, otype):
    """Batched IFFT operation for multiple templates"""
    if itype == 'complex' and otype == 'complex':
        # Stack the input arrays directly 
        batch_data = cupy.stack([v.data for v in invecs])
        
        # Perform batch IFFT
        result = cupy.fft.ifft(batch_data)
        
        # Copy results back efficiently using cupy.copyto
        for i, outvec in enumerate(outvecs):
            cupy.copyto(outvec.data, result[i])
            outvec *= len(outvec)
            
    elif itype == 'complex' and otype == 'real':
        batch_data = cupy.stack([v.data for v in invecs])
        result = cupy.fft.irfft(batch_data)
        
        for i, outvec in enumerate(outvecs):
            cupy.copyto(outvec.data, result[i])
            outvec *= len(outvec)
    else:
        raise ValueError(_INV_FFT_MSG.format("IFFT", itype, otype))

class BatchFFT(_BaseFFT):
    """Class for performing batched FFTs via the cupy interface"""
    def __init__(self, invecs, outvecs, batch_size):
        self.invecs = invecs
        self.outvecs = outvecs
        self.batch_size = batch_size
        self.prec, self.itype, self.otype = _check_fft_args(invecs[0], outvecs[0])

    def execute(self):
        batch_fft(self.invecs, self.outvecs, self.prec, self.itype, self.otype)

class BatchIFFT(_BaseIFFT):
    """Class for performing batched IFFTs via the cupy interface"""
    def __init__(self, invecs, outvecs, batch_size):
        self.invecs = invecs  
        self.outvecs = outvecs
        self.batch_size = batch_size
        self.prec, self.itype, self.otype = _check_fft_args(Array(invecs[0]), Array(outvecs[0]))

    def execute(self):
        batch_ifft(self.invecs, self.outvecs, self.prec, self.itype, self.otype)