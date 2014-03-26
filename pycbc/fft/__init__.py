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

import pycbc
import pycbc.scheme
import pycbc.types
from pycbc.types import Array as _Array
from pycbc.types import TimeSeries as _TimeSeries
from pycbc.types import FrequencySeries as _FrequencySeries

from numpy import dtype

# Package initialization: loop over *possible* backends, trying to
# import them, and adding to the list of available backends. Also set
# the default backend for each of the three schemes, and it is an error
# if we cannot import that.

# Helper function to import a possible module and update available.

def _backend_update(key,possible,available):
    mod = __import__('pycbc.fft.'+possible[key],
                     fromlist=['pycbc.fft'])
    available.update({key:mod})

# The next part is a little tricky.  There are two issues:
#  (1) The logical name for what the user specifies as the backend
#      may not be something we can call the module (e.g., 'numpy')
#  (2) We envision a scenario in which several backends supported
#      in principle may not be available on all platforms (even
#      all CPU platforms)
#
# We also define a module-level function, call that, and delete
# it, so that we don't pollute the module-level namespace with all
# of the variables and functions used in setting this up.

# The global lists, intended for users to check which backends are
# available.  They should be empty if the corresponding scheme
# is not available
cpu_backends = []
cuda_backends = []
opencl_backends = []

# This is a private global variable, a dict of dicts indexed
# by scheme and then by backend
_fft_backends = {}

def _setup_fft():
    # CPU backends; PyCBC requires LAL so it's an error if the following
    # fails:
    import fftw as _cpu_default
    _cpu_possible_backends = {'numpy':'npfft',
                              'lal':'lalfft',
                              'fftw':'fftw'}
    _cpu_backends = {'Default': _cpu_default}
    cpu_backends.append('Default')

    # NOTE: Syntax below for iteration over dict keys should change in
    # Python 3!
    for backend in _cpu_possible_backends.iterkeys():
        try:
            _backend_update(backend,_cpu_possible_backends,_cpu_backends)
            cpu_backends.append(backend)
        except ImportError:
            pass

    # CUDA backends;
    if pycbc.HAVE_CUDA:
        import cufft as _cuda_default
        _cuda_possible_backends = {'cuda' : 'cufft',
                                   'pyfft':'cuda_pyfft'}
        _cuda_backends = {'Default' : _cuda_default}
        cuda_backends.append('Default')

        # NOTE: Syntax below for iteration over dict keys should change in
        # Python 3!
        for backend in _cuda_possible_backends.iterkeys():
            try:
                _backend_update(backend,_cuda_possible_backends,_cuda_backends)
                cuda_backends.append(backend)
            except ImportError:
                pass

    # OpenCL backends; blank for now
    if pycbc.HAVE_OPENCL:
        import cl_pyfft as _opencl_default
        _opencl_possible_backends = {'pyfft' : 'cl_pyfft'}
        _opencl_backends = {'Default': _opencl_default}
        opencl_backends.append('Default')

        # NOTE: Syntax below for iteration over dict keys should change in
        # Python 3!
        for backend in _opencl_possible_backends.iterkeys():
            try:
                _backend_update(backend,_opencl_possible_backends,
                                _opencl_backends)
                opencl_backends.append(backend)
            except ImportError:
                pass

    # Finally, we update our global dict-of-dicts:
    _fft_backends.update({pycbc.scheme.CPUScheme: _cpu_backends})

    if pycbc.HAVE_CUDA:
        _fft_backends.update({pycbc.scheme.CUDAScheme:
                                  _cuda_backends})
    if pycbc.HAVE_OPENCL:
        _fft_backends.update({pycbc.scheme.OpenCLScheme:
                                  _opencl_backends})

# We're finally at the end of setup_fft(), so call it and delete it:
_setup_fft()
del _setup_fft
del _backend_update

# The main purpose of the top-level module is to present a
# uniform interface for a forward and reverse FFT, independent of
# the underlying backend.  We perform sanity checking here, at the
# top-level, and then don't worry about it in submodules. To
# facilitate this checking, we define dicts mapping the numpy dtype
# to the corresponding precisions and types.

def _check_fft_args(invec,outvec):
    if not isinstance(invec,_Array):
        raise TypeError("Input is not a PyCBC Array")
    if not isinstance(outvec,_Array):
        raise TypeError("Output is not a PyCBC Array")

    if isinstance(invec,_TimeSeries) and not isinstance(
        outvec,_FrequencySeries):
        raise TypeError(
            "When input is TimeSeries output must be FrequencySeries")
    if isinstance(outvec,_TimeSeries) and not isinstance(
        invec,_FrequencySeries):
        raise TypeError(
            "When output is TimeSeries input must be FrequencySeries")
    if isinstance(invec,_FrequencySeries) and not isinstance(
        outvec,_TimeSeries):
        raise TypeError(
            "When input is FrequencySeries output must be TimeSeries")
    if isinstance(outvec,_FrequencySeries) and not isinstance(
        invec,_TimeSeries):
        raise TypeError(
            "When output is FrequencySeries input must be TimeSeries")

    iprec = invec.precision
    oprec = outvec.precision
    if iprec != oprec:
        raise ValueError("Input and output precisions must agree")

    itype = invec.kind
    otype = outvec.kind
    return [iprec,itype,otype]

def fft(invec,outvec,backend='Default'):
    """ Fourier transform from invec to outvec.

    Perform a fourier transform. The type of transform is determined
    by the dtype of invec and outvec.

    Parameters
    ----------
    invec : TimeSeries or FrequencySeries
        The input vector.
    outvec : TimeSeries or FrequencySeries
        The output.
    """
    [prec,itype,otype] = _check_fft_args(invec,outvec)
    if itype == 'complex' and otype == 'complex':
        if len(invec) != len(outvec):
            raise ValueError(
                "Lengths of input and output vectors must agree")
    elif itype == 'real' and otype == 'complex':
        if len(outvec) != (len(invec)/2+1):
            raise ValueError(
                "Output length of R2HC must be half input length plus one")
    else:
        raise ValueError("Inconsistent dtypes for forward FFT")
    thescheme = pycbc.scheme.mgr.state.__class__
    thebackend = _fft_backends[thescheme][backend]
    thebackend.fft(invec,outvec,prec,itype,otype)
    # For a forward FFT, the length of the *input* vector is the length
    # we should divide by, whether C2C or R2HC transform
    if isinstance(invec,_TimeSeries):
        outvec._epoch = invec._epoch
        outvec._delta_f = 1.0/(invec._delta_t*len(invec))
        outvec *= invec._delta_t
    elif isinstance(invec,_FrequencySeries):
        outvec._epoch = invec._epoch
        outvec._delta_t = 1.0/(invec._delta_f*len(invec))
        outvec *= invec._delta_f

def ifft(invec, outvec, backend='Default'):
    """ Inverse fourier transform from invec to outvec.

    Perform an inverse fourier transform. The type of transform is determined
    by the dtype of invec and outvec.

    Parameters
    ----------
    invec : TimeSeries or FrequencySeries
        The input vector.
    outvec : TimeSeries or FrequencySeries
        The output.
    """
    [prec,itype,otype] = _check_fft_args(invec,outvec)
    if itype == 'complex' and otype == 'complex':
        if len(invec) != len(outvec):
            raise ValueError(
                "Lengths of input and output vectors must agree")
    elif itype == 'complex' and otype == 'real':
        if len(invec) != (len(outvec)/2+1):
            raise ValueError(
                "Input length of R2HC@r must be half output length plus one")
    else:
        raise ValueError("Inconsistent dtypes for reverse FFT")
    thescheme = pycbc.scheme.mgr.state.__class__
    thebackend = _fft_backends[thescheme][backend]
    thebackend.ifft(invec,outvec,prec,itype,otype)
    # For an inverse FFT, the length of the *output* vector is the length
    # we should divide by, whether C2C or HC2R transform
    if isinstance(invec,_TimeSeries):
        outvec._epoch = invec._epoch
        outvec._delta_f = 1.0/(invec._delta_t*len(outvec))
        outvec *= invec._delta_t
    elif isinstance(invec,_FrequencySeries):
        outvec._epoch = invec._epoch
        outvec._delta_t = 1.0/(invec._delta_f*len(outvec))
        outvec *= invec._delta_f
