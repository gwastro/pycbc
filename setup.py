#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz , Andrew Miller
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
from distutils import sysconfig,file_util
from distutils.core import setup,Command,Extension


# ======== DISTUTILS CONFIGURATION AND BUILD SCRIPTS ==========================   
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
        # Get the list of cpu test modules
        self.test_modules = self.find_test_modules("test*.py")     
        # Run from the build directory
        sys.path.insert(0,self.build_dir)

        test_results.append("\n" + (self.scheme + " tests ").rjust(30))
        for test in self.test_modules:
            a = subprocess.call(['python','test/'+test+'.py','-s', self.scheme])
            if a != 0:
                result_str = str(test).ljust(30) + ": Fail : " + str(a)
            else:
                result_str = str(test).ljust(30) + ": Pass" 
            test_results.append(result_str)

class test(Command):
    sub_commands = [('test_cpu',None),('test_cuda',None),('test_opencl',None)]
    user_options = []
    description = "run the available tests for all compute schemes (cpu,cuda,opencl)"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)
        for test in test_results:
            print test
    
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
	requires = ['swiglal'],
    packages = ['pycbc','pycbc.fft','pycbc.types'],
    scripts = [],
)

