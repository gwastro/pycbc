# Copyright (C) 2012  Alex Nitz, Josh Willis
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
"""PyCBC contains a toolkit for CBC gravitational wave analysis
"""
from __future__ import (absolute_import, print_function)
import subprocess, os, sys, tempfile, signal, warnings

# Filter annoying Cython warnings that serve no good purpose.
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import logging
import random
import string

try:
    # This will fail when pycbc is imported during the build process,
    # before version.py has been generated.
    from .version import git_hash
    from .version import version as pycbc_version
except:
    git_hash = 'none'
    pycbc_version = 'none'

__version__ = pycbc_version


def init_logging(verbose=False, format='%(asctime)s %(message)s'):
    """ Common utility for setting up logging in PyCBC.

    Installs a signal handler such that verbosity can be activated at
    run-time by sending a SIGUSR1 to the process.
    """
    def sig_handler(signum, frame):
        logger = logging.getLogger()
        log_level = logger.level
        if log_level == logging.DEBUG:
            log_level = logging.WARN
        else:
            log_level = logging.DEBUG
        logging.warn('Got signal %d, setting log level to %d',
                     signum, log_level)
        logger.setLevel(log_level)

    signal.signal(signal.SIGUSR1, sig_handler)

    if verbose:
        initial_level = logging.DEBUG
    else:
        initial_level = logging.WARN
    logging.getLogger().setLevel(initial_level)
    logging.basicConfig(format=format, level=initial_level)

def makedir(path):
    """
    Make the analysis directory path and any parent directories that don't
    already exist. Will do nothing if path already exists.
    """
    if path is not None and not os.path.exists(path):
        os.makedirs(path)


# PyCBC-Specific Constants

# Set the value we want any aligned memory calls to use
# N.B.: *Not* all pycbc memory will be aligned to multiples
# of this value

PYCBC_ALIGNMENT = 32

# Dynamic range factor: a large constant for rescaling
# GW strains.  This is 2**69 rounded to 17 sig.fig.

DYN_RANGE_FAC =  5.9029581035870565e+20

# String used to separate parameters in configuration file section headers.
# This is used by the distributions and transforms modules
VARARGS_DELIM = '+'

# Check for optional components of the PyCBC Package
try:
    # This is a crude check to make sure that the driver is installed
    try:
        loaded_modules = subprocess.Popen(['lsmod'], stdout=subprocess.PIPE).communicate()[0]
        loaded_modules = loaded_modules.decode()
        if 'nvidia' not in loaded_modules:
            raise ImportError("nvidia driver may not be installed correctly")
    except OSError:
        pass

    # Check that pycuda is installed and can talk to the driver
    import pycuda.driver as _pycudadrv

    HAVE_CUDA=True
except ImportError:
    HAVE_CUDA=False

# Check for MKL capability
try:
    import pycbc.fft.mkl
    HAVE_MKL=True
except ImportError:
    HAVE_MKL=False

# Check for openmp suppport, currently we pressume it exists, unless on
# platforms (mac) that are silly and don't use the standard gcc.
if sys.platform == 'darwin':
    HAVE_OMP = False
else:
    HAVE_OMP = True

# https://pynative.com/python-generate-random-string/
def random_string(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))
