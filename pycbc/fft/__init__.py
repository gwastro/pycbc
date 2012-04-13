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
This package provides a front-end to various fast Fourier transform
implementations within PyCBC.
"""

import pycbc
import pycbc.scheme
import pycbc.array
from numpy import dtype

# Helper function to import a possible module and update available.

def _backend_update(key,possible,available):
    mod = __import__('pycbc.fft.'+possible[key],
                     fromlist=['pycbc.fft'])
    available.update({key:mod})

# Package initialization: loop over *possible* backends, trying to
# import them, and adding to the list of available backends. Also set
# the default backend for each of the three schemes, and it is an error
# if we cannot import that.

# CPU backends:

import npfft as _cpu_default

# The next part is a little tricky.  There are two issues:
#  (1) The logical name for what the user specifies as the backend
#      may not be something we can call the module (e.g., 'numpy')
#  (2) We envision a scenario in which several backends supported
#      in principle may not be available on all platforms (even
#      all CPU platforms)
_cpu_possible_backends = {'numpy':'npfft',
                          'lal':'swiglalfft'}
_cpu_backends = {None: _cpu_default}

# NOTE: Syntax below for iteration over dict keys should change in
# Python 3!
for backend in _cpu_possible_backends.iterkeys():
    try:
        _backend_update(backend,_cpu_possible_backends,_cpu_backends)
    except ImportError:
        pass

# CUDA backends; blank for now
if pycbc.HAVE_CUDA:
    import cufft as _cuda_default
    _cuda_possible_backends = {'cuda' : 'cufft'}
    _cuda_backends = {None : _cuda_default}

   # NOTE: Syntax below for iteration over dict keys should change in
   # Python 3!
    for backend in _cuda_possible_backends.iterkeys():
        try:
        _backend_update(backend,_cuda_possible_backends,_cuda_backends)
        except ImportError:
            pass

# OpenCL backends; blank for now
if pycbc.HAVE_OPENCL:
    import cldefault as _opencl_default
    _opencl_possible_backends = {'opencl' : 'cldefault'}
    _opencl_backends = {None: _opencl_default}

   # NOTE: Syntax below for iteration over dict keys should change in
   # Python 3!
    for backend in _opencl_possible_backends.iterkeys():
        try:
        _backend_update(backend,_opencl_possible_backends,
                        _opencl_backends)
        except ImportError:
            pass

# Now create a dict-of-dicts of backends

_fft_backends = {None: _cpu_backends}
if pycbc.HAVE_CUDA:
    _fft_backends.update({pycbc.scheme.CUDAScheme:
                          _cuda_backends})
if pycbc.HAVE_OPENCL:
    _fft_backends.update({pycbc.scheme.OpenCLScheme:
                              _opencl_backends})


# The main purpose of the top-level module is to present a
# uniform interface for a forward and reverse FFT, independent of
# the underlying backend.  We perform sanity checking here, at the
# top-level, and then don't worry about it in submodules. To
# facilitate this checking, we define dicts mapping the numpy dtype
# to the corresponding precisions and types.

_type_dict = {dtype('float32'):'real',
              dtype('float64'):'real',
              dtype('complex64'):'complex',
              dtype('complex128'):'complex'}
_prec_dict = {dtype('float32'):'single',
              dtype('float64'):'double',
              dtype('complex64'):'single',
              dtype('complex128'):'double'}

def _check_fft_args(invec,outvec):
    if not isinstance(invec,pycbc.array.Array):
        raise TypeError("Input is not a PyCBC Array")
    if not isinstance(outvec,pycbc.array.Array):
        raise TypeError("Output is not a PyCBC Array")
    # Later we should add checks for TimeSeries and
    # FrequencySeries: input and output should map one
    # to the other.
    iprec = _prec_dict[invec.dtype]
    oprec = _prec_dict[outvec.dtype]
    if iprec is not oprec:
        raise TypeError("Input and output precisions must agree")
    else:
        prec = iprec
    itype = _type_dict[invec.dtype]
    otype = _type_dict[outvec.dtype]
    return [prec,itype,otype]

def fft(invec,outvec,backend=None):
    [prec,itype,otype] = _check_fft_args(invec,outvec)
    if itype is 'complex' and otype is 'complex':
        if len(invec) is not len(outvec):
            raise ValueError(
                "Lengths of input and output vectors must agree")
    elif itype is 'real' and otype is 'complex':
        if len(outvec) is not (len(invec)/2+1):
            raise ValueError(
                "Output length of R2HC must be half input length plus one")
    else:
        raise ValueError("Inconsistent types for forward FFT")
    thescheme = pycbc.scheme.mgr.state
    thebackend = _fft_backends[thescheme][backend]
    thebackend.fft(invec,outvec,prec,itype,otype)

def ifft(invec,outvec,backend=None):
    [prec,itype,otype] = _check_fft_args(invec,outvec)
    if itype is 'complex' and otype is 'complex':
        if len(invec) is not len(outvec):
            raise ValueError(
                "Lengths of input and output vectors must agree")
    elif itype is 'complex' and otype is 'real':
        if len(invec) is not (len(outvec)/2+1):
            raise ValueError(
                "Input length of R2HC@r must be half output length plus one")
    else:
        raise ValueError("Inconsistent types for reverse FFT")
    thescheme = pycbc.scheme.mgr.state
    thebackend = _fft_backends[thescheme][backend]
    thebackend.ifft(invec,outvec,prec,itype,otype)
