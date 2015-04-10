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
from pkg_resources import VersionConflict

import pycbc
import pycbc.scheme
import pycbc.types
from pycbc.types import Array as _Array
from pycbc.types import TimeSeries as _TimeSeries
from pycbc.types import FrequencySeries as _FrequencySeries

from numpy import dtype

from optparse import OptionGroup

# These are global variables, that are modified by the various scheme-
# dependent submodules, to maintain a list of all possible backends
# for all possible schemes that are available at runtime.  This list
# and dict are then used when parsing command-line options.

_all_backends_list = []
_all_backends_dict = {}

# The following is empty by default, and only used if the command line
# option parsing sets it

_default_backends_list = []

# The following helper function is in this top-level module because it
# is used by the scheme-dependent files to write their version of the
# _available_backends() function.

def _list_available(possible_list,possible_dict):
    # It possibly is strange that we have both a list and a dict.
    # The reason for this is that the name the user specfies for a
    # backend, e.g. 'numpy', may be something that the fft submodule
    # cannot be called, so we need a dict mapping those names to the
    # actual names of the modules.  However when we iterate, it must
    # be in the defined order, because that represents a preference
    # for which backends are most likely to be preferrable.  As a
    # dict is unordered, we cannot simply use its keys for this purpose.
    available_list = []
    available_dict = {}
    for backend in possible_list:
        try:
            mod = __import__('pycbc.fft.'+possible_dict[backend],fromlist=['pycbc.fft'])
            available_dict.update({backend:mod})
            available_list.append(backend)
        except (ImportError, OSError, VersionConflict):
            pass
    return available_list, available_dict

# The following is the function called by each scheme's setup to add whatever new
# backends may have been found to the global list.  Since some backends may be
# shared, we must first check to make sure that the item in the list is not already
# in the global list, and we assume that the keys to the dict are in one-to-one
# correspondence with the items in the list.

def _update_global_available(new_list,new_dict,global_list,global_dict):
    for item in new_list:
        if item not in global_list:
            global_list.append(item)
            global_dict.update({item:new_dict[item]})

# We import the three modules for the underlying architectures; each of
# these files may in priniciple provide implementations of the scheme-
# dependent functions for several schemes. We don't need to *use* these
# sub-modules directly, but must import them for the side-effects they
# have on the schemes themselves (adding the availble backends and
# importing those modules)

import fft_cpu as _fft_cpu
import fft_cuda as _fft_cuda
import fft_opencl as _fft_opencl

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

def fft(invec,outvec,backends=[]):
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
    backends_dict = pycbc.scheme.mgr.state.fft_backends_dict
    backends_list = pycbc.scheme.mgr.state.fft_backends_list
    if len(backends) > 0:
        for backend in backends:
            if backend in backends_list:
                thebackend = backends_dict[backend]
                break
        else:
            raise RuntimeError("None of specified FFT backends available")
    elif len(_default_backends_list) > 0:
        for backend in _default_backends_list:
            if backend in backends_list:
                thebackend = backends_dict[backend]
                break
        else:
            raise RuntimeError("None of the command-line specified backends available")
    else:
        # This will fail if the backends_list is empty, but then so it should
        thebackend = backends_dict[backends_list[0]]
    # The following line is where all the work is done:
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

def ifft(invec, outvec, backends=[]):
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
    backends_dict = pycbc.scheme.mgr.state.fft_backends_dict
    backends_list = pycbc.scheme.mgr.state.fft_backends_list
    if len(backends) > 0:
        for backend in backends:
            if backend in backends_list:
                thebackend = backends_dict[backend]
                break
        else:
            raise RuntimeError("None of specified FFT backends available")
    elif len(_default_backends_list) > 0:
        for backend in _default_backends_list:
            if backend in backends_list:
                thebackend = backends_dict[backend]
                break
        else:
            raise RuntimeError("None of the command-line specified backends available")
    else:
        # This will fail if the backends_list is empty, but then so it should
        thebackend = backends_dict[backends_list[0]]
    # The following line is where all the work is done:
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
# from the command line.

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
    # We have one argument to specify the backends.  This becomes the default list used
    # if none is specified for a particular call of fft() of ifft().  Note that this
    # argument expects a *list* of inputs, as indicated by the nargs='*'.
    fft_group.add_argument("--fft-backends",
                      help="Preference list of the FFT backends. "
                           "Choices are: \n" + str(_all_backends_list),
                      nargs='*',default=[])

    for backend in _all_backends_dict.values():
        try:
            backend.insert_fft_options(fft_group)
        except AttributeError:
            pass

def verify_fft_options(opt, parser):
    """Parses the FFT options and verifies that they are
       reasonable.


    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes.
    parser : object
        OptionParser instance.
    """

    if len(opt.fft_backends) > 0:
        for backend in opt.fft_backends:
            if backend not in _all_backends_list:
                parser.error("Backend {0} is not available".format(backend))

    for backend in _all_backends_dict.values():
        try:
            backend.verify_fft_options(opt,parser)
        except AttributeError:
            pass

def from_cli(opt):
    """Parses the command line options and sets the FFT backend
    for each (available) scheme. Aside from setting the default
    backed for this context, this function will also call (if
    it exists) the from_cli function of the specified backends in
    the *current* scheme; typically one would only call this function
    once inside of a scheme context manager, but if it is desired
    to perform FFTs both inside and outside of a context, then
    this function would need to be called again.

    Parameters
    ----------
    opt: object
        Result of parsing the CLI with OptionParser, or any object with
        the required attributes.

    Returns
    """
    global _default_backends_list

    if len(opt.fft_backends) > 0:
        _default_backends_list = opt.fft_backends

    backends_dict = pycbc.scheme.mgr.state.fft_backends_dict
    backends_list = pycbc.scheme.mgr.state.fft_backends_list

    for backend in _default_backends_list:
        if backend in backends_list:
            try:
                backends_dict[backend].from_cli(opt)
            except AttributeError:
                pass
