#!/usr/bin/python
# Copyright (C) 2011 Gergely Debrenczeni, Karsten Wiesner
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
setup.py file for swig wrap pycbcopencl into pycbc
"""

from distutils.core import setup, Extension
import os; 

# Setting the ocmpier to use g++

os.environ['CC'] = 'g++'; 

vector_module = Extension('_pycbcopencl',
                          sources=['pycbcopencl_wrap.c','pycbcopencl.c'],
                          )


# base source files that do not require special libraries
pycbcopencl_swig_sources = ['pycbcopencl.i']
pycbcopencl_c_sources = ['pycbcopencl.c']
gpu_inspiral_gpuutils_c_sources = ['gpu_inspiral_gpuutils.c']


# define the extension module
pycbcopencl_ext = Extension( '_pycbcopencl', 
  sources = pycbcopencl_swig_sources + pycbcopencl_c_sources + gpu_inspiral_gpuutils_c_sources,
  depends = ['pycbcopencl.h','gpu_inspiral_gpuutils.h'],
  swig_opts = [],
  include_dirs = ['/usr/local/nvidia/sdk-3.2/OpenCL/common/inc/'],
  extra_compile_args = ['-Wall','-fPIC'],
  library_dirs = ['/usr/lib'],
  libraries = ['OpenCL','stdc++'])

setup (name = 'pycbcopencl',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap pycbcopencl""",
       ext_modules = [pycbcopencl_ext],
       py_modules = ["pycbcopencl"],
       )
