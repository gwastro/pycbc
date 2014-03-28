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
import ctypes, ctypes.util

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
    
# Determine whether we can use aligned memory, and if so define
# an aligned_malloc function

HAVE_ALIGNED_MALLOC = False
try:
    libc = ctypes.cdll.LoadLibrary(ctypes.util.find_library("c"))
    if libc is not None:
        if hasattr(libc,'posix_memalign'):
            HAVE_ALIGNED_MALLOC = True
except:
    pass
    
if HAVE_ALIGNED_MALLOC:
    MEM_ALIGNMENT = 32
    def amalloc(n):
        am_func = libc.posix_memalign
        am_func.argtypes = [ctypes.POINTER(ctypes.c_void_p),
                            ctypes.c_ulonglong,ctypes.c_ulonglong]
        newmem = ctypes.c_void_p()
        am_func(ctypes.byref(newmem),MEM_ALIGNMENT,n)
        return newmem

# PYCBC Specfic Constants

DYN_RANGE_FAC =  5.9029581035870565e+20


