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
import os
from trace import fullmodname
import unittest
from distutils import sysconfig,file_util
from distutils.core import setup,Command,Extension


# Below function copied from PyCBC v1 setup.py file.  Not sure whom to credit--Duncan?
def pkg_config(libraries=[],library_dirs=[],include_dirs=[],pkg_libraries=[]):
    # Check that we have the packages
    for pkg in pkg_libraries:
        if os.system('pkg-config --exists %s 2>/dev/null' % pkg) == 0:
            pass
        else:
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
    version = '0.1',
    description = 'Gravitational wave CBC analysis toolkit',
    author = 'Ligo Virgo Collaboration - PyCBC team',
    url = 'https://sugwg-git.phy.syr.edu/dokuwiki/doku.php?id=pycbc:home',
    cmdclass = { 'test'  : test , 'test_cpu':test_cpu,'test_cuda':test_cuda,
                 'test_opencl':test_opencl},
    ext_modules = [],
    requires = ['lal'],
    packages = ['pycbc','pycbc.fft','pycbc.types','pycbc.filter'],
    scripts = [],
)

