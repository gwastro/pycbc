#!/usr/bin/python
# Copyright (C) 2011 Karsten Wiesner, Duncan Brown, Alex Nitz
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
import os, fnmatch, sys,commands
from trace import fullmodname
import unittest
from distutils import sysconfig,file_util
from distutils.core import setup,Command,Extension
from distutils.command.clean import clean as _clean
from distutils.command.build import build as _build
from distutils.command.install import install as _install
from distutils.command.build_ext import build_ext as _build_ext
from distutils.command.config import config as _config
from distutils.command.install_lib import install_lib as _install_lib
import log

# make sure we use c99. unfortunately there is not better way to do this in
# distutils than setting the environment variable
try:
  os.environ['CFLAGS'] += '--std=c99'
except KeyError:
  os.environ['CFLAGS'] = '--std=c99'


ver = '0.1'
pycbc_clean_files = []

# Create the swig wrapped extension library along with the unwrapped shared library
def pycbc_add_extensions(path,name,sources,
                        swig_source,
                        libraries=[],
                        depends=[], 
                        include_dirs=[],
                        extra_compile_args=[],
                        swig_opts=[],
                        library_dirs=[],
                        pkg_libraries=[]):
    # Check that we have the packages
    for pkg in pkg_libraries:
        if os.system('pkg-config --exists %s 2>/dev/null' % pkg) == 0:
            log.good("found %s library", pkg)
        else:
            log.fatal("could not find %s library",pkg)
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

    #Create the shared library           
    pycbc_lib=Extension("lib"+name,sources=sources,
                        libraries=libraries,
                        depends=depends,
                        library_dirs=library_dirs,
                        include_dirs=include_dirs,
                        extra_compile_args=extra_compile_args,
                        swig_opts=swig_opts)
    
    #Create the swig wrapped extension
    pycbc_ext=Extension(path+"._"+name,sources=sources + swig_source,
                    libraries=libraries,
                    depends=depends,
                    library_dirs=library_dirs,
                    include_dirs=include_dirs,
                    extra_compile_args=extra_compile_args,
                    swig_opts=swig_opts)
                    
    
                    
    return [pycbc_lib,pycbc_ext]


# ======== CPU extension modules for the top-level package ====================
pycbc_extensions = []

pycbc_extensions+=pycbc_add_extensions( 'pycbc.clayer','pycbccpu', 
        sources = ['pycbc/clayer/cpu/pycbccpu.c','pycbc/clayer/except.c'], 
        swig_source = ['pycbc/clayer/cpu/pycbccpu.i'],
        include_dirs= ['pycbc/clayer','./'],
        depends = ['pycbc/clayer/cpu/pycbccpu.h','pycbc/clayer/cpu/pycbccpu_private.h'],
        extra_compile_args = ['-Wall','-fPIC'],
        swig_opts = ['-outdir','pycbc/clayer'],
    )

pycbc_clean_files.append('pycbc/clayer/cpu.py')
pycbc_clean_files.append('pycbc/clayer/cpu/pycbccpu_wrap.c')


pycbc_extensions+=pycbc_add_extensions( 'pycbc.datavector.clayer','datavectorcpu',
    sources = ['pycbc/datavector/clayer/cpu/datavectorcpu.c'],
    swig_source = ['pycbc/datavector/clayer/cpu/datavectorcpu.i'],
    depends = ['pycbc/datavector/clayer/cpu/datavectorcpu.h',
               'pycbc/datavector/clayer/cpu/datavectorcpu_private.h',
               'pycbc/datavector/clayer/datavector.h',
               'pycbc/clayer/cpu/pycbccpu.h'],
    include_dirs = ['pycbc/clayer/cpu', './',
                    'pycbc/datavector/clayer'],
    libraries = ['pycbccpu'],
    swig_opts = ['-outdir','pycbc/datavector/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    )
    

pycbc_clean_files.append('pycbc/datavector/clayer/cpu.py')
pycbc_clean_files.append('pycbc/datavector/clayer/cpu/datavectorcpu_wrap.c')


pycbc_extensions+=pycbc_add_extensions( 'pycbc.straindata.clayer','straindatacpu', 
    sources = ['pycbc/straindata/clayer/cpu/straindatacpu.c'],
    swig_source = ['pycbc/straindata/clayer/cpu/straindatacpu.i'],
    depends = ['pycbc/straindata/clayer/cpu/straindatacpu.h',
               'pycbc/straindata/clayer/cpu/straindatacpu_private.h'],
    include_dirs = ['pycbc/clayer/cpu',
                    'pycbc/datavector/clayer/cpu','pycbc/datavector/clayer'],
    swig_opts = ['-outdir','pycbc/straindata/clayer'],
    libraries = ['pycbccpu'],
    extra_compile_args = ['-Wall','-fPIC']
    )

pycbc_clean_files.append('pycbc/straindata/clayer/cpu.py')
pycbc_clean_files.append('pycbc/straindata/clayer/cpu/straindatacpu_wrap.c')


pycbc_extensions+=pycbc_add_extensions( 'pycbc.templatebank.clayer','templatebankcpu',
    sources = ['pycbc/templatebank/clayer/cpu/templatebankcpu.c'],
    swig_source = ['pycbc/templatebank/clayer/cpu/templatebankcpu.i'],
    depends = ['pycbc/templatebank/clayer/cpu/templatebankcpu.h',
               'pycbc/templatebank/clayer/cpu/templatebankcpu_private.h'],
    include_dirs = ['pycbc/clayer/cpu',
               'pycbc/datavector/clayer/cpu','pycbc/datavector/clayer'],
    swig_opts = ['-outdir','pycbc/templatebank/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    )

pycbc_clean_files.append('pycbc/templatebank/clayer/cpu.py')
pycbc_clean_files.append('pycbc/templatebank/clayer/cpu/templatebankcpu_wrap.c')


pycbc_extensions+=pycbc_add_extensions( 'pycbc.matchedfilter.clayer','matchedfiltercpu', 
    sources = ['pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.c'],
    swig_source = ['pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.i'],
    depends = ['pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.h',
               'pycbc/matchedfilter/clayer/cpu/matchedfiltercpu_private.h'],
    include_dirs = ['pycbc/clayer/cpu',
               'pycbc/datavector/clayer/cpu','pycbc/datavector/clayer'],
    swig_opts = ['-outdir','pycbc/matchedfilter/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    )

pycbc_clean_files.append('pycbc/matchedfilter/clayer/cpu.py')
pycbc_clean_files.append('pycbc/matchedfilter/clayer/cpu/matchedfiltercpu_wrap.c')

# ======== (END) CPU extension modules for the top-level package ==============

# ======== OPENCL extension modules ===========================================

# THIS IS JUST AN EXAMPLE RIGHT NOW

pycbc_opencl_extensions = []

pycbc_opencl_extensions+=pycbc_add_extensions( 'pycbc.clayer','pycbcopencl', 
        sources = ['pycbc/clayer/opencl/pycbcopencl.c',
                    'pycbc/clayer/opencl/openclexcept.c'], 
        swig_source = ['pycbc/clayer/opencl/pycbcopencl.i'],
        include_dirs= ['pycbc/clayer','./',
                     'pycbc/clayer/opencl'],
        depends = ['pycbc/clayer/opencl/pycbcopencl.h',
                    'pycbc/clayer/opencl/pycbcopencl_private.h',
                    'pycbc/clayer/opencl/gpu_inspiral_utils.h'],
        extra_compile_args = ['-Wall','-fPIC'],
        libraries=['pycbccpu','OpenCL'],
        pkg_libraries=[],
        swig_opts = ['-outdir','pycbc/clayer'],
    )

pycbc_clean_files.append('pycbc/opencl/opencl.py')
pycbc_clean_files.append('pycbc/clayer/cpu/pycbcopencl_wrap.c')



# ======== (END) OPENCL extension modules =====================================

# ======== CUDA   extension modules ===========================================

pycbc_cuda_extensions = []

# ======== (END) CUDA extension modules =======================================



# ======== DISTUTILS CONFIGURATION AND BUILD SCRIPTS ==========================

# Configure the pycbc package #FIXME, this is temporary
class config(_config):
    user_options = [('with-opencl=', None,"specify the opencl toolkit location")]
    def initialize_options(self):
        _config.initialize_options(self)
        self.with_opencl=None
    def finalize_options(self):
        _config.finalize_options(self)
        #Check for opencl headers
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
        # self.test_modules+= self.test_opencl()
        
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
          
# Set the libpycbc location in the extensions modules' rpath    
class build_ext(_build_ext):
    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.install_dir = None
        self.opencl_dir = None 
    def finalize_options(self):
        _build_ext.finalize_options(self)
        self.set_undefined_options('install',('install_lib','install_dir')) 
        self.set_undefined_options('config',('with_opencl','opencl_dir'))   
    def run(self):
        self.rpath.append(self.build_lib)
        self.rpath.append(self.install_dir)
        self.library_dirs.append(self.build_lib)
        
        #Add opencl header to include_dir if we are building those modules
        if (self.opencl_dir != None):
            self.include_dirs+=[self.opencl_dir+"/include"]
            
        _build_ext.run(self)

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
    url = 'https://sugwg-git.phy.syr.edu/dokuwiki/doku.php?id=pycbc:home',
    cmdclass = { 'config' : config,
                 'clean' : clean,
                 'build' : build,
                 'install':install,
		         'test'  : test ,
		         'build_ext': build_ext},
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
    scripts = ['bin/pycbc_min_cpu_pipeline'],
)

