#!/usr/bin/python
# Copyright (C) 2011 Josh Willis
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
setup.py file for swig wrap fftcpu into pycbc
"""

from distutils.core import setup, Extension

fftcpu_module = Extension('_fftcpu',
                          sources=['fftcpu_wrap.c','fftcpu.c'],
                          )


# base source files that do not require special libraries
fftcpu_swig_sources = ['fftcpu.i']
fftcpu_c_sources = ['fftcpu.c']


# define the extension module
fftcpu_ext = Extension( '_fftcpu', 
  sources = fftcpu_swig_sources + fftcpu_c_sources,
  depends = ['fftcpu.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'fftcpu',
       version = '0.1',
       author = "Josh Willis",
       description = """Swig wrap fftcpu""",
       ext_modules = [fftcpu_ext],
       py_modules = ["fftcpu"],
       )
