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
import subprocess, os, sys, signal, warnings

# Filter annoying Cython warnings that serve no good purpose.
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import logging
import random
import string
from datetime import datetime as dt

try:
    # This will fail when pycbc is imported during the build process,
    # before version.py has been generated.
    from .version import git_hash
    from .version import version as pycbc_version
except:
    git_hash = 'none'
    pycbc_version = 'none'

__version__ = pycbc_version


class LogFormatter(logging.Formatter):
    """
    Format the logging appropriately
    This will return the log time in the ISO 6801 standard,
    but with millisecond precision
    https://en.wikipedia.org/wiki/ISO_8601
    e.g. 2022-11-18T09:53:01.554+00:00
    """
    converter = dt.fromtimestamp

    def formatTime(self, record, datefmt=None):
        ct = self.converter(record.created).astimezone()
        t = ct.strftime("%Y-%m-%dT%H:%M:%S")
        s = f"{t}.{int(record.msecs):03d}"
        timezone = ct.strftime('%z')
        timezone_colon = f"{timezone[:-2]}:{timezone[-2:]}"
        s += timezone_colon
        return s


def init_logging(verbose=False, format='%(asctime)s %(message)s'):
    """Common utility for setting up logging in PyCBC.

    Installs a signal handler such that verbosity can be activated at
    run-time by sending a SIGUSR1 to the process.

    Parameters
    ----------
    verbose : bool or int, optional
        What level to set the verbosity level to. Accepts either a boolean
        or an integer representing the level to set. If True/False will set to
        ``logging.INFO``/``logging.WARN``. For higher logging levels, pass
        an integer representing the level to set (see the ``logging`` module
        for details). Default is ``False`` (``logging.WARN``).
    format : str, optional
        The format to use for logging messages.
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

    if not verbose:
        initial_level = logging.WARN
    elif int(verbose) == 1:
        initial_level = logging.INFO
    else:
        initial_level = int(verbose)

    logger = logging.getLogger()
    logger.setLevel(initial_level)
    sh = logging.StreamHandler()
    logger.addHandler(sh)
    sh.setFormatter(LogFormatter(fmt=format))


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

    # MacosX after python3.7 switched to 'spawn', however, this does not
    # preserve common state information which we have relied on when using
    # multiprocessing based pools.
    import multiprocessing
    if hasattr(multiprocessing, 'set_start_method'):
        multiprocessing.set_start_method('fork')
else:
    HAVE_OMP = True

# https://pynative.com/python-generate-random-string/
def random_string(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def gps_now():
    """Return the current GPS time as a float using Astropy.
    """
    from astropy.time import Time

    return float(Time.now().gps)
