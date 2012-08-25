#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz, Andrew Miller, Josh Willis
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
import os
import fnmatch
import sys
import subprocess
import commands
from trace import fullmodname
import unittest
from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup,Command,Extension,find_packages
from distutils.command.clean import clean as _clean
from distutils.command.build import build as _build
from numpy import get_include as np_get_include


# Below function copied from PyCBC v1 setup.py file.  Not sure whom to credit--Duncan?
def pkg_config(libraries=[],library_dirs=[],include_dirs=[],pkg_libraries=[]):
    # Check that we have the packages
    for pkg in pkg_libraries:
        if os.system('pkg-config --exists %s 2>/dev/null' % pkg) == 0:
            pass
        else:
            print "Could not find library {0}".format(pkg)
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

# Now use the above function to set up our extension library's needs:

ext_libraries, ext_library_dirs, ext_include_dirs = pkg_config(pkg_libraries=["lal","lalsimulation","lalinspiral"])

# Setup our swig options. We need the first two to match with how the swiglal
# wrappings are compiled, so that our module can talk to theirs.  Then we must
# reassemble the include_dirs to add the "-I" so swig can find things
ext_swig_opts = ['-O','-keyword','-builtin','-outdir','pycbc']
for libpath in ext_include_dirs:
    ext_swig_opts.append('-I'+str(libpath))

# Numpy include must be added separately:

ext_include_dirs += [np_get_include()]

# Define our extension module

lalwrap_module = Extension('_lalwrap',
                           sources=['include/lalwrap.i'],
                           depends=['include/pycbc_laltypes.i'],
                           swig_opts=ext_swig_opts,
                           include_dirs=ext_include_dirs,
                           library_dirs=ext_library_dirs,
                           libraries=ext_libraries,
                           extra_compile_args=['-std=c99']
                           )

testlalwrap_module = Extension('_testlalwrap',
                               sources=['include/testlalwrap.i'],
                               depends=['include/pycbc_laltypes.i'],
                               swig_opts=ext_swig_opts,
                               include_dirs=ext_include_dirs,
                               library_dirs=ext_library_dirs,
                               libraries=ext_libraries,
                               extra_compile_args=['-std=c99']
                               )

# Add swig-generated files to the list of things to clean, so they
# get regenerated each time.
class clean(_clean):
    def finalize_options (self):
        _clean.finalize_options(self)
        self.clean_files = ['pycbc/lalwrap.py','pycbc/lalwrap.pyc','include/lalwrap_wrap.c',
                            'pycbc/testlalwrap.py','pycbc/testlalwrap.pyc','include/testlalwrap_wrap.c']

    def run(self):
        _clean.run(self)
        for f in self.clean_files:
            try:
                os.unlink(f)
                print 'removed {0}'.format(f)
            except:
                pass

# Override build order, so swig is handled first.
class build(_build):
    # override order of build sub-commands to work around
    # <http://bugs.python.org/issue7562>
    sub_commands = [ ('build_clib', _build.has_c_libraries),
                     ('build_ext', _build.has_ext_modules),
                     ('build_py', _build.has_pure_modules),
                     ('build_scripts', _build.has_scripts) ]

    def run(self):
        _build.run(self)

test_results = []
# Run all of the testing scripts
class TestBase(Command):
    user_options = []
    test_modules = []
    def initialize_options(self):
        self.scheme = None
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
        self.run_command('build')
        # Get the list of cpu test modules
        self.test_modules = self.find_test_modules("test*.py")
        # Run from the build directory
        os.environ['PYTHONPATH'] = self.build_dir + ":" + os.environ['PYTHONPATH']

        test_results.append("\n" + (self.scheme + " tests ").rjust(30))
        for test in self.test_modules:
            test_command = 'python ' + 'test/' + test + '.py -s ' + self.scheme
            a = subprocess.call(test_command,env=os.environ,shell=True)
            if a != 0:
                result_str = str(test).ljust(30) + ": Fail : " + str(a)
            else:
                result_str = str(test).ljust(30) + ": Pass"
            test_results.append(result_str)

        for test in test_results:
            print test


class test(Command):
    def has_cuda(self):
        try:
            import pycuda
            return True
        except:
            return False

    def has_opencl(self):
        try:
            import pyopencl
            return True
        except:
            return False

    sub_commands = [('test_cpu',None),('test_cuda',has_cuda),('test_opencl',has_opencl)]
    user_options = []
    description = "run the available tests for all compute schemes (cpu,cuda,opencl)"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)

class test_cpu(TestBase):
    description = "run all CPU tests"
    def initialize_options(self):
        TestBase.initialize_options(self)
        self.scheme = 'cpu'

class test_cuda(TestBase):
    description = "run CUDA tests"
    def initialize_options(self):
        TestBase.initialize_options(self)
        self.scheme = 'cuda'

class test_opencl(TestBase):
    description = "run OpenCL tests"
    def initialize_options(self):
        TestBase.initialize_options(self)
        self.scheme = 'opencl'

# do the actual work of building the package
setup (
    name = 'PyCBC',
    version = '0.0.1',
    py_modules = ['distribute_setup'],
    description = 'Gravitational wave CBC analysis toolkit',
    author = 'Ligo Virgo Collaboration - PyCBC team',
    url = 'https://sugwg-git.phy.syr.edu/dokuwiki/doku.php?id=pycbc:home',
    cmdclass = { 'test'  : test,
                 'test_cpu':test_cpu,
                 'test_cuda':test_cuda,
                 'test_opencl':test_opencl,
                 'clean' : clean,
                 'build' : build},
    ext_modules = [lalwrap_module,testlalwrap_module],
    install_requires = [''],
    packages = find_packages(),
    scripts = [],
)

