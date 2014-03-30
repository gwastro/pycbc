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

from optparse import OptionGroup

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

# This will be indexed by scheme, and return a list of strings
# naming available backends for that scheme
_backend_names = {}

# Ordinarily, the user will not pass a backend to the fft/ifft
# functions.  Instead, he or she will set the default backend
# for a particular scheme, and rely on that.  Below is the
# machinery to handle that.  Initial values of the defaults
# are here.
_default_backends = {pycbc.scheme.CPUScheme: 'fftw'}

if pycbc.HAVE_CUDA:
    _default_backends.update({pycbc.scheme.CUDAScheme: 'cuda'})

if pycbc.HAVE_OPENCL:
    _default_backends.update({pycbc.scheme.OpenCLScheme: 'pyfft'})

def get_default_backend(scheme):
    return _default_backends[scheme]

def set_default_backend(scheme,backend_name):
    global _default_backends
    if backend_name not in _backend_names[scheme]:
        raise ValueError("{0} is an invalid choice for default FFT backend".format(backend_name))
    _default_backends.update({scheme:backend_name})

# This is a private global variable, a dict of dicts indexed
# by scheme and then by backend.  This double indexing yields
# actual imported modules that can be called
_fft_backends = {}

def _setup_fft():
    # CPU backends
    _cpu_possible_backends = {'numpy':'npfft',
                              'lal':'lalfft',
                              'mkl':'mkl',
                              'fftw':'fftw',
                             }
    _cpu_backends = {}

    # NOTE: Syntax below for iteration over dict keys should change in
    # Python 3!
    for backend in _cpu_possible_backends.iterkeys():
        try:
            _backend_update(backend,_cpu_possible_backends,_cpu_backends)
            cpu_backends.append(backend)
        except (ImportError, OSError):
            pass

    _backend_names.update({pycbc.scheme.CPUScheme: cpu_backends})

    # CUDA backends;
    if pycbc.HAVE_CUDA:
        _cuda_possible_backends = {'cuda' : 'cufft',
                                   'pyfft':'cuda_pyfft'}
        _cuda_backends = {}

        # NOTE: Syntax below for iteration over dict keys should change in
        # Python 3!
        for backend in _cuda_possible_backends.iterkeys():
            try:
                _backend_update(backend,_cuda_possible_backends,_cuda_backends)
                cuda_backends.append(backend)
            except (ImportError, OSError):
                pass

        _backend_names.update({pycbc.scheme.CUDAScheme: cuda_backends})

    # OpenCL backends; blank for now
    if pycbc.HAVE_OPENCL:
        _opencl_possible_backends = {'pyfft' : 'cl_pyfft'}
        _opencl_backends = {}

        # NOTE: Syntax below for iteration over dict keys should change in
        # Python 3!
        for backend in _opencl_possible_backends.iterkeys():
            try:
                _backend_update(backend,_opencl_possible_backends,
                                _opencl_backends)
                opencl_backends.append(backend)
            except (ImportError, OSError):
                pass

        _backend_names.update({pycbc.scheme.OpenCLScheme: opencl_backends})

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

def fft(invec,outvec,backend=None):
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
    if backend is None:
        backend = get_default_backend(thescheme)
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

def ifft(invec, outvec, backend=None):
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
    if backend is None:
        backend = get_default_backend(thescheme)
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

def insert_fft_option_group(parser,scheme):
    """
    Adds the options used to choose an FFT backend. This should be used
    if your program supports the ability to select the FFT backend; otherwise
    you may simply call the fft and ifft functions and rely on default
    choices.  This function must be passed a scheme (so that option must be
    parsed first) since the possible choices depend on that scheme.
    
    Parameters
    ----------
    parser : object
        OptionParser instance
    scheme: object
        pycbc.scheme scheme
    """
    fft_group = parser.add_argument_group("Options for selecting the"
                                          " FFT backend in this program.")   
    fft_group.add_argument("--fft-backend", 
                      help="The choice of FFT backend. "
                           "Choices are " + str(_backend_names[scheme]), 
                      choices=_backend_names[scheme], 
                      default=_default_backends[scheme])                                                          

def from_cli(opt,scheme):
    """Parses the command line options and returns the FFT backend.
    Note that the return value is the actual module (which may be
    further queried for its own options). Also note that the default
    backend for provided scheme will also be set from this command
    line option.

    Parameters
    ----------
    opt: object
        Result of parsing the CLI with OptionParser, or any object with
        the required attributes.

    scheme: Scheme
        A pycbc.scheme type whose backend will be set and returned by
        this cli
        
    Returns
    -------
    backend: module
        The submodule of pycbc.fft that is the FFT backend for the
        specified scheme.
    """
    thebackend = _fft_backends[scheme][opt.fft_backend]
    set_fft_backend(scheme,opt.fft_backend)
    return thebackend
    
def verify_fft_options(opt, parser, scheme):
    """Parses the  processing scheme options and verifies that they are 
       reasonable. 
       
  
    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes.
    parser : object
        OptionParser instance.
    scheme: object
        A pycbc.scheme type whose backend fft option will be verified
    """
    if opt.fft_backend not in _backend_names[scheme]:
        parser.error("{0} is not a valid FFT backend.".format(opt.fft_backend))
