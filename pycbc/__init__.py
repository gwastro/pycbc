# Copyright (C) 2012  Alex Nitz
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

# Check for optional components of the PyCBC Package

try:
    import pycuda.driver as _pycudadrv
    HAVE_CUDA=True
except ImportError:
    HAVE_CUDA=False
    
try:
    import pyopencl as _pyopencl
    HAVE_OPENCL=True
except ImportError:
    HAVE_OPENCL=False
    
    
# PYCBC Specfic Constants

DYN_RANGE_FAC =  5.9029581035870565e+20


