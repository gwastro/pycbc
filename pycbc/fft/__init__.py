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

# Ordinarily, the user will not pass a backend to the fft/ifft
# functions.  Instead, he or she will set the default backend
# for a particular scheme, and rely on that.  Below is the
# machinery to handle that.  Initial values of the defaults
# are here.
_default_backends = {pycbc.scheme.CPUScheme: 'fftw'}
# Note that the following dict is filled in when
# _setup_fft runs
_backend_names = {}

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

# The function that will set all of the lists and dicts of what is really 
# available (could be successfully imported).
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

    _backend_names.update({pycbc.scheme.CPUScheme:cpu_backends})

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

        _backend_names.update({pycbc.scheme.CUDAScheme:cuda_backends})

    # OpenCL backends
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

        _backend_names.update({pycbc.scheme.OpenCLScheme:opencl_backends})


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

# Next we add all of the machinery to set backends and their options
# from the command line.  There's a single function that will define
# an option for setting the backend for each of the possible schemes,
# since a given script might need to setup running under each of them.
# That function also attempts to set up any options that any of the
# backends themselves define.  Likewise, another function validates
# all of the options.  The function to actually act from the command
# line, called from_cli(), expects in addition to a parsed cli to be
# given a scheme.  It therefore expects to be called *after* the scheme
# has been determined, via the pycbc.scheme module's CLI interface.

def insert_fft_option_group(parser):
    """
    Adds the options used to choose an FFT backend. This should be used
    if your program supports the ability to select the FFT backend; otherwise
    you may simply call the fft and ifft functions and rely on default
    choices.  This function will also attempt to add any options exported
    by available backends through a function called insert_fft_options.
    These submodule functions should take the fft_group object as argument.
    
    Parameters
    ----------
    parser : object
        OptionParser instance
    """
    fft_group = parser.add_argument_group("Options for selecting the"
                                          " FFT backend and controlling its performance"
                                          " in this program.")   
    # We have options for each of the backends, because a given script
    # may not know on which platform it will run and may want to set a
    # default backend and options for that backend for more than one
    fft_group.add_argument("--fft-cpu-backend", 
                      help="Determines the FFT CPU backend. "
                           "Choices are: \n" + str(cpu_backends), 
                      choices=cpu_backends, 
                      default=_default_backends[pycbc.scheme.CPUScheme])

    for backend in _fft_backends[pycbc.scheme.CPUScheme].values():
        try:
            backend.insert_fft_options(fft_group)
        except AttributeError:
            pass

    if pycbc.HAVE_CUDA:
        fft_group.add_argument("--fft-cuda-backend", 
                               help="Determines the FFT CUDA backend. "
                               "Choices are: \n" + str(cuda_backends), 
                               choices=cuda_backends, 
                               default=_default_backends[pycbc.scheme.CUDAScheme])

        for backend in _fft_backends[pycbc.scheme.CUDAScheme].values():
            try:
                backend.insert_fft_options(fft_group)
            except AttributeError:
                pass

    if pycbc.HAVE_OPENCL:
        fft_group.add_argument("--fft-opencl-backend", 
                               help="Determines the FFT OpenCL backend. "
                               "Choices are: \n" + str(opencl_backends), 
                               choices=opencl_backends, 
                               default=_default_backends[pycbc.scheme.OpenCLScheme])

        for backend in _fft_backends[pycbc.scheme.OpenCLScheme].values():
            try:
                backend.insert_fft_options(fft_group)
            except AttributeError:
                pass

def verify_fft_options(opt, parser):
    """Parses the  processing scheme options and verifies that they are 
       reasonable. 
       
  
    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes.
    parser : object
        OptionParser instance.
    """
    try:
        if opt.fft_cpu_backend not in cpu_backends:
            parser.error("{0} is not a valid CPU FFT backend.".format(opt.fft_cpu_backend))
    except AttributeError:
        pass

    for backend in _fft_backends[pycbc.scheme.CPUScheme].values():
        try:
            backend.verify_fft_options(opt,parser)
        except AttributeError:
            pass

    if pycbc.HAVE_CUDA:
        try:
            if opt.fft_cuda_backend not in cuda_backends:
                parser.error("{0} is not a valid CUDA FFT backend.".format(opt.fft_cuda_backend))
        except AttributeError:
            pass

        for backend in _fft_backends[pycbc.scheme.CUDAScheme].values():
            try:
                backend.verify_fft_options(opt,parser)
            except AttributeError:
                pass

    if pycbc.HAVE_OPENCL:
        try:
            if opt.fft_opencl_backend not in opencl_backends:
                parser.error("{0} is not a valid OpenCL FFT backend.".format(opt.fft_opencl_backend))
        except AttributeError:
            pass

        for backend in _fft_backends[pycbc.scheme.OpenCLScheme].values():
            try:
                backend.verify_fft_options(opt,parser)
            except AttributeError:
                pass

def from_cli(opt,ctx):
    """Parses the command line options and sets the FFT backend
    for the provided context.  This function may be called more
    than once with different instances of ctx to set the backend
    an options for multiple schemes (if needed). Aside from setting
    the default backed for this context, this function will also
    call (if it exists) the from_cli function of the specified
    backend (that function should only take opt as an argument).

    Parameters
    ----------
    opt: object
        Result of parsing the CLI with OptionParser, or any object with
        the required attributes.

    ctx: Scheme
        An instance of a pycbc.scheme type, whose backend will be set
        and the corresponding options of the backend also parsed from
        the cli.
        
    Returns
    -------
    kwdrets: dict
        A dictionary containing keyword/value pairs returned by the
        call to backend.from_cli.  If that backend has no return
        values, this dict will be empty (i.e., the top level never
        returns anything; it just sets the default backend).
    """
    kwdrets = {}

    thescheme = type(ctx)
    thebackend = _fft_backends[thescheme][opt.fft_backend]
    set_fft_backend(thescheme,opt.fft_backend)

    try:
        tmpdict = thebackend.from_cli(opt)
        if tmpdict is not None:
            kwdrets.update(tmpdict)
    except AttributeError:
        pass

    return kwdrets
    
