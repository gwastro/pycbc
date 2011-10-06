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
setup.py file for swig wrap datavectorcpu into pycbc
"""

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext
from distutils import file_util


class pycbc_build_ext(build_ext):
    
    def run(self):
        build_ext.run(self)
        # copy the python stub for the _xyz.so ext library to the 
        # build/lib location. this won't be done by the build_ext cmd.
        file_util.copy_file('datavectorcpu.py', self.build_lib)


datavectorcpu_ext = Extension(
    name= '_datavectorcpu', 
    sources = ['datavectorcpu.i',
               'datavectorcpu.c'],
    #not documented option: swig_opts = [],
    include_dirs = ['.','../'],
    #define_macros=[('TESTMACRO', '1')],
    #undef_macros=['TESTMACRO'],
    library_dirs = [],
    libraries = [],
    runtime_library_dirs = [],
    extra_objects = [],
    extra_compile_args = ['-Wall'],
    extra_link_args = [],
    export_symbols = [],
    depends = ['datavectorcpu.c',
               'datavectorcpu.i',
               'datavectorcpu_types.h',
               'datavectorcpu_prototypes.h'],
    language = []
)


setup (name = 'datavectorcpu',
       version = '0.1',
       author = "Karsten Wiesner",
       description = 'swig wrapped python extension datavectorcpu',
       ext_modules = [datavectorcpu_ext],
       cmdclass = {
            'build_ext' : pycbc_build_ext }  # extend build_ext for pycbc

       #py_modules = ["datavectorcpu"],   don't list the swig generated 
       #python extension stub here
       )                                   
