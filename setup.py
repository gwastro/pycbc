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
setup.py file for PyCBC package
"""
import os, fnmatch, sys
from trace import fullmodname
import unittest
from distutils.core import setup,Command
from distutils.core import Extension
from distutils.command.clean import clean as _clean
from distutils.command.build import build as _build
from distutils.command.install import install as _install

ver = '0.1'
pycbc_extensions = []
pycbc_clean_files = []

platform_sources = ['pycbc/clayer/except.c','pycbc/clayer/log.c']
platform_depends = ['pycbc/clayer/except.h','pycbc/clayer/log.h']
platform_include_dirs = ['./']

# extension modules for the top-level package
pycbc_extensions.append(Extension( 'pycbc.clayer._cpu', 
    sources = ['pycbc/clayer/cpu/pycbccpu.i',
               'pycbc/clayer/cpu/pycbccpu.c'],
    depends = ['pycbc/clayer/cpu/pycbccpu.h',
               'pycbc/clayer/cpu/pycbccpu_private.h'],
    swig_opts = ['-outdir','pycbc/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    ))

pycbc_clean_files.append('pycbc/clayer/cpu.py')
pycbc_clean_files.append('pycbc/clayer/cpu/pycbccpu_wrap.c')


pycbc_extensions.append(Extension( 'pycbc.datavector.clayer._cpu', 
    sources = platform_sources + ['pycbc/datavector/clayer/cpu/datavectorcpu.i',
               'pycbc/datavector/clayer/cpu/datavectorcpu.c'],
    depends = platform_depends + ['pycbc/datavector/clayer/cpu/datavectorcpu.h',
               'pycbc/datavector/clayer/cpu/datavectorcpu_private.h',
               'pycbc/datavector/clayer/datavector.h',
               'pycbc/clayer/cpu/pycbccpu.h'],
    include_dirs = platform_include_dirs + ['pycbc/clayer/cpu',
                    'pycbc/datavector/clayer'],
    swig_opts = ['-outdir','pycbc/datavector/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    ))

pycbc_clean_files.append('pycbc/datavector/clayer/cpu.py')
pycbc_clean_files.append('pycbc/datavector/clayer/cpu/datavectorcpu_wrap.c')


pycbc_extensions.append(Extension( 'pycbc.straindata.clayer._cpu', 
    sources = platform_sources + ['pycbc/straindata/clayer/cpu/straindatacpu.i',
               'pycbc/straindata/clayer/cpu/straindatacpu.c'],
    depends = platform_depends + ['pycbc/straindata/clayer/cpu/straindatacpu.h',
               'pycbc/straindata/clayer/cpu/straindatacpu_private.h'],
    include_dirs = platform_include_dirs + ['pycbc/clayer/cpu',
                    'pycbc/datavector/clayer/cpu'],
    swig_opts = ['-outdir','pycbc/straindata/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    ))

pycbc_clean_files.append('pycbc/straindata/clayer/cpu.py')
pycbc_clean_files.append('pycbc/straindata/clayer/cpu/straindatacpu_wrap.c')


pycbc_extensions.append(Extension( 'pycbc.templatebank.clayer._cpu', 
    sources = ['pycbc/templatebank/clayer/cpu/templatebankcpu.i',
               'pycbc/templatebank/clayer/cpu/templatebankcpu.c'],
    depends = ['pycbc/templatebank/clayer/cpu/templatebankcpu.h',
               'pycbc/templatebank/clayer/cpu/templatebankcpu_private.h'],
    include_dirs = ['pycbc/clayer/cpu',
               'pycbc/datavector/clayer/cpu'],
    swig_opts = ['-outdir','pycbc/templatebank/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    ))

pycbc_clean_files.append('pycbc/templatebank/clayer/cpu.py')
pycbc_clean_files.append('pycbc/templatebank/clayer/cpu/templatebankcpu_wrap.c')


pycbc_extensions.append(Extension( 'pycbc.matchedfilter.clayer._cpu', 
    sources = ['pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.i',
               'pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.c'],
    depends = ['pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.h',
               'pycbc/matchedfilter/clayer/cpu/matchedfiltercpu_private.h'],
    include_dirs = ['pycbc/clayer/cpu',
               'pycbc/datavector/clayer/cpu'],
    swig_opts = ['-outdir','pycbc/matchedfilter/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    ))

pycbc_clean_files.append('pycbc/matchedfilter/clayer/cpu.py')
pycbc_clean_files.append('pycbc/matchedfilter/clayer/cpu/matchedfiltercpu_wrap.c')

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

        # Get the list of cpu test modules
        self.test_modules+= self.test_cpu()
        
        #FIXME WHEN WE HAVE THE CONFIGURATION DONE 
        # self.test_modules+= self.test_cuda()
        #self.test_modules+= self.test_opencl()
        
        # Run from the build directory
        sys.path.insert(0,self.build_dir)

        # load all of the tests into a suite
        suite = unittest.TestLoader().loadTestsFromNames(self.test_modules)

        # run test suite
        results=unittest.TextTestRunner().run(suite)
        if( not results.wasSuccessful()):
	    raise SystemExit(1)    
	

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
                print 'removed {0}'.format(f)
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
        
# Override the standard install to ensure that all the unittests pass
class install(_install):
    def run(self):      
        if (not self.force):  
            self.run_command('test')
        _install.run(self)

# do the actual work of building the package
setup (
    name = 'pycbc',
    version = ver,
    description = 'Gravitational wave CBC analysis toolkit',
    author = 'Ligo Virgo Collaboration - pyCBC team',
    author_email = 'https://sugwg-git.phy.syr.edu/dokuwiki/doku.php?id=pycbc:home',
    cmdclass = { 'clean' : clean,
                 'build' : build,
                 'install':install,
		         'test'  : test },
    ext_modules = pycbc_extensions,
    packages = ['pycbc','pycbc.clayer',
                'pycbc.datavector','pycbc.datavector.clayer',
                'pycbc.bandpass',
                'pycbc.chisqveto',
                'pycbc.injection',
                'pycbc.overwhiteningfilter',
                'pycbc.resample',
                'pycbc.singledetectorevent',
                'pycbc.straindata','pycbc.straindata.clayer',
                'pycbc.templatebank','pycbc.templatebank.clayer',
                'pycbc.matchedfilter','pycbc.matchedfilter.clayer',
                'pycbc.waveform'],
    scripts = ['bin/pycbc_min_cpu_pipeline']
)

