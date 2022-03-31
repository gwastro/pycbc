# Copyright (C) 2012  Josh Willis
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
This module provides the numpy backend of the fast Fourier transform
for the PyCBC package.
"""

import logging
import numpy.fft
from .core import _check_fft_args
from .core import _BaseFFT, _BaseIFFT

_INV_FFT_MSG = ("I cannot perform an {} between data with an input type of "
                "{} and an output type of {}")

def fft(invec, outvec, _, itype, otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("numpy backend of pycbc.fft does not "
                                  "support in-place transforms")
    if itype == 'complex' and otype == 'complex':
        outvec.data[:] = numpy.asarray(numpy.fft.fft(invec.data),
                                       dtype=outvec.dtype)
    elif itype == 'real' and otype == 'complex':
        outvec.data[:] = numpy.asarray(numpy.fft.rfft(invec.data),
                                       dtype=outvec.dtype)
    else:
        raise ValueError(_INV_FFT_MSG.format("FFT", itype, otype))


def ifft(invec, outvec, _, itype, otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("numpy backend of pycbc.fft does not "
                                  "support in-place transforms")
    if itype == 'complex' and otype == 'complex':
        outvec.data[:] = numpy.asarray(numpy.fft.ifft(invec.data),
                                       dtype=outvec.dtype)
        outvec *= len(outvec)
    elif itype == 'complex' and otype == 'real':
        outvec.data[:] = numpy.asarray(numpy.fft.irfft(invec.data,len(outvec)),
                                       dtype=outvec.dtype)
        outvec *= len(outvec)
    else:
        raise ValueError(_INV_FFT_MSG.format("IFFT", itype, otype))


WARN_MSG = ("You are using the class-based PyCBC FFT API, with the numpy "
            "backed. This is provided for convenience only. If performance is "
            "important use the class-based API with one of the other backends "
            "(for e.g. MKL or FFTW)")


class FFT(_BaseFFT):
    """
    Class for performing FFTs via the numpy interface.
    """
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(FFT, self).__init__(invec, outvec, nbatch, size)
        logging.warning(WARN_MSG)
        self.prec, self.itype, self.otype = _check_fft_args(invec, outvec)

    def execute(self):
        fft(self.invec, self.outvec, self.prec, self.itype, self.otype)


class IFFT(_BaseIFFT):
    """
    Class for performing IFFTs via the numpy interface.
    """
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(IFFT, self).__init__(invec, outvec, nbatch, size)
        logging.warning(WARN_MSG)
        self.prec, self.itype, self.otype = _check_fft_args(invec, outvec)

    def execute(self):
        ifft(self.invec, self.outvec, self.prec, self.itype, self.otype)
