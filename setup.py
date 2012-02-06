#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz
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


"""
setup.py file for PyCBC package
"""
import os, fnmatch, sys,commands
from trace import fullmodname
import unittest
from distutils import sysconfig,file_util
from distutils.core import setup,Command,Extension
from distutils.command.clean import clean as _clean
from distutils.command.build import build as _build
from distutils.command.install import install as _install
from distutils.command.config import config as _config
from distutils.command.install_lib import install_lib as _install_lib
import logging

# make sure we use c99. unfortunately there is not better way to do this in
# distutils than setting the environment variable
try:
  os.environ['CFLAGS'] += '--std=c99'
except KeyError:
  os.environ['CFLAGS'] = '--std=c99'


# Get build information from pkg-config
def pkg_config(libraries=[],library_dirs=[],include_dirs=[],pkg_libraries=[]):
    # Check that we have the packages
    for pkg in pkg_libraries:
        if os.system('pkg-config --exists %s 2>/dev/null' % pkg) == 0:
            logging.good("found %s library", pkg)
        else:
            logging.fatal("could not find %s library",pkg)
            sys.exit(1)
        
    # Get the pck-config flags
    if len(pkg_libraries)>0 : 
        for token in commands.getoutput("pkg-config --libs --cflags %s" % ' '.join(pkg_libraries)).split():
            if token[:2]== "-l":
                libraries.append(token[2:])
            if token[:2]== "-L":
                library_dirs.append(token[2:])
            if token[:2] == "-I":
                include_dirs.append(token[2:])
            
    return libraries,library_dirs,include_dirs


pycbc_extensions=[]
pycbc_clean_files = []
# ======== CPU Extension Modules ====================

pycbc_extensions.append(Extension('pycbc.example.ext.cpu._hello_cpu',
    sources = ['pycbc/example/ext/cpu/hello.c','pycbc/example/ext/cpu/hello.i'],
    depends = ['pycbc/example/ext/cpu/hello.h'],
    include_dirs = ['include'],
    libraries = [],
    swig_opts = ['-outdir','pycbc/example/ext/'],
    extra_compile_args = ['-Wall','-fPIC',]
    ))
    
pycbc_clean_files+=['pycbc/example/ext/cpu/hello_wrap.c','pycbc/example/ext/cpu/']
    
# ========= OpenCL Extension Modules ================

# ========= Cuda Extension Modules ==================

# ========= Scripts and Programs ================
pycbc_scripts=[]


# ======== DISTUTILS CONFIGURATION AND BUILD SCRIPTS ==========================

# Configure the pycbc package #FIXME, this is temporary
class config(_config):
    user_options = [('with-opencl=', None,"specify the opencl toolkit location")]
    def initialize_options(self):
        _config.initialize_options(self)
        self.with_opencl=None
    def finalize_options(self):
        _config.finalize_options(self)
        
        #Check for OPENCL headers
        if (self.with_opencl!=None):
            log.warn("looking for opencl headers")
            if not self.check_header("CL/opencl.h",include_dirs=[self.with_opencl+"/include"]):
                log.warn("could not find opencl.h")
                self.with_opencl=None
            else:
                log.good("found opencl headers")
                #Everthing looks ok, build the opencl modules
                self.distribution.ext_modules+=pycbc_opencl_extensions
        else:
            log.warn("will not build opencl modules")
            self.with_opencl=None
    
        #Check for CUDA
        
# Run all of the testing scripts
class test(Command):
    user_options = []
    test_modules = []
    def initialize_options(self):
        self.build_dir = None
    def finalize_options(self):
        #Populate the needed variables
        self.set_undefined_options('build',('build_lib', 'build_dir'))
        
    def find_test_modules(self,pattern):
       # Find all the unittests that match a given string pattern
        modules= []
        for path, dirs, files in os.walk("test"):
            for filename in fnmatch.filter(files, pattern):
                #add the test directories to the path
                sys.path.append(os.path.join(path))
                #save the module name for importing
                modules.append(fullmodname(filename))
        return modules
        
    def run(self):
        # Ensure that we have everything built first
        self.run_command('build')

        # Get the list of cpu test modules
        self.test_modules+= self.find_test_modules("test*.py")
        
        # Run from the build directory
        sys.path.insert(0,self.build_dir)

        # load all of the tests into a suite
        suite = unittest.TestLoader().loadTestsFromNames(self.test_modules)

        # run test suite
        results=unittest.TextTestRunner(verbosity=2).run(suite)

# class to clean up the files automatically generated by swig
class clean(_clean):
    def finalize_options (self):
        _clean.finalize_options(self)
        self.clean_files = pycbc_clean_files

    def run(self):
        _clean.run(self)
        for f in self.clean_files:
            try:
                os.unlink(f)
                print('removed {0}'.format(f))
            except:
                pass

class build(_build):
    # override order of build sub-commands to work around
    # <http://bugs.python.org/issue7562>
    sub_commands = [ ('build_clib', _build.has_c_libraries),
                     ('build_ext', _build.has_ext_modules),
                     ('build_py', _build.has_pure_modules),
                     ('build_scripts', _build.has_scripts) ]

    def run(self):
        _build.run(self)

# do the actual work of building the package
setup (
    name = 'PyCBC',
    version = '0.1',
    description = 'Gravitational wave CBC analysis toolkit',
    author = 'Ligo Virgo Collaboration - PyCBC team',
    url = 'https://sugwg-git.phy.syr.edu/dokuwiki/doku.php?id=pycbc:home',
    cmdclass = { 'config' : config,
                 'clean' : clean,
                 'build' : build,
		         'test'  : test ,},
    ext_modules = [],
    packages = ['pycbc'],
    scripts = pycbc_scripts,
)

