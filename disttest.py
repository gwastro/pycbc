#!/usr/bin/env python
import os, fnmatch,sys
from trace import fullmodname
import unittest
from distutils.core import Command

# Run all of the testing scripts
class disttest(Command):
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
        for path, dirs, files in os.walk(os.getcwd()):
            for filename in fnmatch.filter(files, pattern):
                #add the test directories to the path
                sys.path.append(os.path.join(path))
                #save the module name for importing
                modules.append(fullmodname(filename))
        return modules
    def test_cpu(self):
        return self.find_test_modules("test_*cpu*.py")
    def test_opencl(self):
        return self.find_test_modules("test_*opencl*.py")
    def test_cuda(self):
        return self.find_test_modules("test_*cuda*.py")
    def run(self):
        # Ensure that we have everything built first
        self.run_command('build')

        # Get the list of test modules
        self.test_modules+= self.test_cpu()
        
        #FIXME WHEN WE HAVE THE CONFIGURATION DONE 
        # self.test_modules+= self.test_cuda()
        # self.test_modules+= self.test_opencl()
        
        # Run from the build directory
        sys.path.insert(0,self.build_dir)

        # load all of the tests into a suite
        suite = unittest.TestLoader().loadTestsFromNames(self.test_modules)

        # run test suite
        unittest.TextTestRunner().run(suite)
