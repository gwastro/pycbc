#!/usr/bin/python
# Copyright (C) 2011 Karsten Wiesner
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

"""
setup.py file for swig wrap datavectoropencl into pycbc
"""

from distutils.core import setup, Extension
import os;

# Setting the ocmpier to use g++

os.environ['CC'] = 'g++';
#os.environ['CXX'] = 'g++';
#os.environ['CPP'] = 'g++';
#os.environ['LDSHARED'] = 'g++';


vector_module = Extension('_datavectoropencl',
                          sources=['datavectoropencl_wrap.c','datavectoropencl.c'],
                          )


# base source files that do not require special libraries
datavectoropencl_swig_sources = ['datavectoropencl.i']
datavectoropencl_c_sources = ['datavectoropencl.c']


# define the extension module
datavectoropencl_ext = Extension( '_datavectoropencl', 
  sources = datavectoropencl_swig_sources + datavectoropencl_c_sources,
  depends = ['datavectoropencl.h'],
  swig_opts = [],
  include_dirs = ['/usr/local/nvidia/sdk-3.2/OpenCL/common/inc/','../../../clayer/opencl'],
  extra_compile_args = ['-Wall','-fPIC'],
  library_dirs = ['/usr/lib','/usr/lib','../../../'],
  libraries = ['pycbcopencl','OpenCL','stdc++'])

setup (name = 'datavectoropencl',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap datavectoropencl""",
       ext_modules = [datavectoropencl_ext],
       py_modules = ["datavectoropencl"],
       )
