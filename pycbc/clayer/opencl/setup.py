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

import os

from distutils import log
from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext
from distutils.command.clean import clean
from distutils import file_util

# use g++ for c++ code
os.environ['CC'] = 'g++';

class pycbc_build_ext(build_ext):
    
    def run(self):
        build_ext.run(self)
        # copy the python stub for the _xyz.so ext library to the 
        # build/lib location. this won't be done by the build_ext cmd.
        file_util.copy_file('pycbcopencl.py', self.build_lib)
        
        # copy files updir to ensure test run from package dir        
        file_util.copy_file('pycbcopencl.py', '../../')
        shared_lib_ext= os.path.join(self.build_lib, '_pycbcopencl.so')
        file_util.copy_file(shared_lib_ext, '../../' )


class pycbc_clean(clean):
    
    def run(self):
        
        self.all = 1       # clean everything in build/
        clean.run(self)
                           # clean specific files
        self.clean_files = ['pycbcopencl.py','pycbcopencl_wrap.c',
                            '../../pycbcopencl.py', '../../_pycbcopencl.so']
                            
        for f in self.clean_files:
            print "removing {0}".format(f)
            try:
                os.unlink(f)
            except:
                log.warn("'%s' does not exist -- can't clean it" % f)


pycbcopencl_ext = Extension(
    name= '_pycbcopencl', 
    sources = ['pycbcopencl.i',
               'pycbcopencl.c',
               'gpu_inspiral_gpuutils.c'],
    swig_opts = [],                     # add include dirs for swig here
    include_dirs = ['/usr/local/nvidia/sdk-3.2/OpenCL/common/inc/',
                    '/Developer/GPU\ Computing/OpenCL/common/inc/'], # add include dirs for compilers here
    #define_macros=[('TESTMACRO', '1')],
    #undef_macros=['TESTMACRO'],
    library_dirs = ['/usr/lib','/usr/lib'],
    libraries = ['OpenCL','stdc++'],
    runtime_library_dirs = [],
    extra_objects = [],
    extra_compile_args = ['-Wall','-fPIC'],
    extra_link_args = [],
    export_symbols = [],
    depends = ['pycbcopencl.c',
               'pycbcopencl.i',
               'pycbcopencl.h',
               'pycbcopencl_private.h',
               'gpu_inspiral_gpuutils.c',
               'gpu_inspiral_gpuutils.h'],
    language = []
)


setup (name = 'pycbcopencl',
       version = '0.1',
       author = "Karsten Wiesner",
       description = 'swig wrapped python extension pycbcopencl',
       ext_modules = [pycbcopencl_ext],
       cmdclass = {
            'build_ext' : pycbc_build_ext,
            'clean' : pycbc_clean 
       }

       #py_modules = ["pycbcopencl"],   don't list the swig generated 
       #python extension stub here
       )                                   
