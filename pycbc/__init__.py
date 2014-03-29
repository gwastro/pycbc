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
import subprocess, os

# Check for optional components of the PyCBC Package
try:
    # This is a crude check to make sure that the driver is installed
    try:
        err = subprocess.call(["nvidia-smi"], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        if err != 0:
            raise ImportError("Cannot access 'nvidia-smi', driver may not be installed correctly")
    except OSError:
        pass

    # Check that pycuda is installed and can talk to the driver
    import pycuda.driver as _pycudadrv

    HAVE_CUDA=True 
except ImportError:
    HAVE_CUDA=False
    
try:
    # This is a crude check to make sure that the driver is installed
    try:
        err = subprocess.call(["nvidia-smi"], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        if err != 0:
            raise ImportError("Cannot access 'nvidia-smi', driver may not be installed correctly")
    except OSError:
        pass

    import pyopencl as _pyopencl
    import pyfft.cl as _pyfftcl
    HAVE_OPENCL=True
except ImportError:
    HAVE_OPENCL=False
    
# Set the value we want any aligned memory calls to use
# N.B.: *Not* all pycbc memory will be aligned to multiples
# of this value

PYCBC_ALIGNMENT = 32

# PYCBC Specfic Constants

DYN_RANGE_FAC =  5.9029581035870565e+20


