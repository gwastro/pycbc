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

from .backend_support import get_backend_modules, get_backend_names
from .backend_support import set_backend, get_backend

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
                           "Choices are: \n" + str(get_backend_names()),
                      nargs='*', default=[])

    for backend in get_backend_modules():
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
        _all_backends = get_backend_names()
        for backend in opt.fft_backends:
            if backend not in _all_backends:
                parser.error("Backend {0} is not available".format(backend))

    for backend in get_backend_modules():
        try:
            backend.verify_fft_options(opt, parser)
        except AttributeError:
            pass

# The following function is the only one that is designed
# only to work with the active scheme.  We'd like to fix that,
# eventually, but it's non-trivial because of how poorly MKL
# and FFTW cooperate.

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

    set_backend(opt.fft_backends)

    # Eventually, we need to be able to parse command lines
    # from more than just the current scheme's preference. But
    # the big problem is that calling from_cli for more than one
    # backend could cause interference; apparently, FFTW and MKL
    # don't play nice unless FFTW has been compiled and linked
    # with icc (and possibly numpy, scipy, and/or Python as well?)

    backend = get_backend()
    try:
        backend.from_cli(opt)
    except AttributeError:
        pass
