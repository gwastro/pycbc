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
setup.py file for swig wrap fft_fftw into pycbc
"""

from distutils.core import setup, Extension

fft_fftw_module = Extension('_fft_fftw',
                          sources=['fft_fftw_wrap.c','fft_fftw.c'],
                          )


# base source files that do not require special libraries
fft_fftw_swig_sources = ['fft_fftw.i']
fft_fftw_c_sources = ['fft_fftw.c']


# define the extension module
fft_fftw_ext = Extension( '_fft_fftw', 
  sources = fft_fftw_swig_sources + fft_fftw_c_sources,
  depends = ['fft_fftw.h','fft_fftw_private.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = ['fftw3','fftw3f'])

setup (name = 'fft_fftw',
       version = '0.1',
       author = "Josh Willis",
       description = """Swig wrap fft_fftw""",
       ext_modules = [fft_fftw_ext],
       py_modules = ["fft_fftw"],
       )
