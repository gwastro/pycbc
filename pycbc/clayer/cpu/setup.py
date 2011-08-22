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
#

"""
setup.py file for swig wrap pycbccpu into pycbc
"""

from distutils.core import setup, Extension

vector_module = Extension('_pycbccpu',
                          sources=['pycbccpu_wrap.c','pycbccpu.c'],
                          )


# base source files that do not require special libraries
pycbccpu_swig_sources = ['pycbccpu.i']
pycbccpu_c_sources = ['pycbccpu.c']


# define the extension module
pycbccpu_ext = Extension( '_pycbccpu', 
  sources = pycbccpu_swig_sources + pycbccpu_c_sources,
  depends = ['pycbccpu.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'pycbccpu',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap pycbccpu""",
       ext_modules = [pycbccpu_ext],
       py_modules = ["pycbccpu"],
       )
