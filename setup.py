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

def pkgconfig(*packages, **kw):
    flag_map = {'-I': 'include_dirs','-L':'library_dirs','-l':'libraries'}
    for token in commands.getoutput("pkg-config --libs --clflags %s" % ' '.join(packages)).split():
        kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
    return kw

ver = '0.1'
pycbc_clean_files = []
pycbc_static_libraries = []
pycbc_headers = []

# ======== CPU extension modules for the top-level package ====================
pycbc_extensions = []

pycbc_extensions.append(Extension( 'pycbc.clayer._cpu', 
    sources = ['pycbc/clayer/cpu/pycbccpu.i',
               'pycbc/clayer/cpu/pycbccpu.c'],
    libraries = ['pycbc'],
    depends = ['pycbc/clayer/cpu/pycbccpu.h',
               'pycbc/clayer/cpu/pycbccpu_private.h'],
    swig_opts = ['-outdir','pycbc/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    ))

pycbc_clean_files.append('pycbc/clayer/cpu.py')
pycbc_clean_files.append('pycbc/clayer/cpu/pycbccpu_wrap.c')


pycbc_extensions.append(Extension( 'pycbc.datavector.clayer._cpu',
    sources = ['pycbc/datavector/clayer/cpu/datavectorcpu.i',
               'pycbc/datavector/clayer/cpu/datavectorcpu.c'],
    depends = ['pycbc/datavector/clayer/cpu/datavectorcpu.h',
               'pycbc/datavector/clayer/cpu/datavectorcpu_private.h',
               'pycbc/datavector/clayer/datavector.h',
               'pycbc/clayer/cpu/pycbccpu.h'],
    include_dirs = ['pycbc/clayer/cpu', './',
                    'pycbc/datavector/clayer'],
    libraries = ['pycbc'],
    swig_opts = ['-outdir','pycbc/datavector/clayer'],
    extra_compile_args = ['-Wall','-fPIC']
    ))

pycbc_clean_files.append('pycbc/datavector/clayer/cpu.py')
pycbc_clean_files.append('pycbc/datavector/clayer/cpu/datavectorcpu_wrap.c')


pycbc_extensions.append(Extension( 'pycbc.straindata.clayer._cpu', 
    sources = ['pycbc/straindata/clayer/cpu/straindatacpu.i',
               'pycbc/straindata/clayer/cpu/straindatacpu.c'],
    depends = ['pycbc/straindata/clayer/cpu/straindatacpu.h',
               'pycbc/straindata/clayer/cpu/straindatacpu_private.h'],
    include_dirs = ['pycbc/clayer/cpu',
                    'pycbc/datavector/clayer/cpu'],
    swig_opts = ['-outdir','pycbc/straindata/clayer'],
    libraries = ['pycbc'],
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

# Find all of the pycbc c source files and build into libpycbc.so
pycbc_sources=['pycbc/clayer/except.c','pycbc/clayer/log.c']
pycbc_include_dirs=['./']
pycbc_depends=['pycbc/clayer/except.h','pycbc/clayer/log.h']
pycbc_libraries=[]
pycbc_library_dirs=[]
for ext in pycbc_extensions:
    for source in ext.sources:
        if ".c" in source:
            pycbc_sources.append(source)
    pycbc_include_dirs += ext.include_dirs
    pycbc_depends += ext.depends
    pycbc_libraries += ext.libraries
    pycbc_library_dirs += ext.library_dirs

pycbc_headers+=pycbc_depends

pycbc_extensions.insert(0,Extension( 'pycbc.libpycbc', 
    sources = pycbc_sources,
    include_dirs = pycbc_include_dirs,
    depends = pycbc_depends,
    libraries = pycbc_libraries,
    library_dirs = pycbc_library_dirs,
    extra_compile_args = ['-Wall','-fPIC']
    ))

# create libpycbc.a

pycbc_include_dirs.append(sysconfig.get_python_inc())

pycbc_static_libraries += [ 'pycbc', { 
    'sources' : pycbc_sources,
    'include_dirs' : pycbc_include_dirs,
    'macros' : [] }],


# ======== (END) CPU extension modules for the top-level package ==============

# ======== OPENCL extension modules ===========================================

pycbc_opencl_extensions = []

# ======== (END) OPENCL extension modules =====================================

# ======== CUDA   extension modules ===========================================

pycbc_cuda_extensions = []

# ======== (END) CUDA extension modules =======================================



# ======== DISTUTILS CONFIGURATION AND BUILD SCRIPTS ==========================

# Configure the pycbc package, locate libraries 
class config(_config):
  def run(self):
    # REQUIRED PACKAGES
    # Check For Lal
    # Check For FFTW3
    print "CONFIG"
    #OPTIONAL PACKAGES
    # Check for CUDA
    # Check for OPENCL
    
        
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
          
# Set the libpycbc location in the extensions modules' rpath    
class build_ext(_build_ext):
    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.install_dir = None
    def finalize_options(self):
        _build_ext.finalize_options(self)
        self.set_undefined_options('install',('install_lib','install_dir'))   
    def run(self):
        self.rpath.append(self.build_lib)
        self.rpath.append(self.install_dir)
        self.library_dirs.append(self.build_lib)
        _build_ext.run(self)

class build(_build):
    # override order of build sub-commands to work around
    # <http://bugs.python.org/issue7562>
    sub_commands = [ ('build_clib', _build.has_c_libraries),
                     ('build_ext', _build.has_ext_modules),
                     ('build_py', _build.has_pure_modules),
                     ('build_scripts', _build.has_scripts) ]

    def run(self):
        self.run_command('config')
        _build.run(self)

        
# Override the standard install to ensure that all the unittests pass
class install(_install):
    def run(self):      
        if (not self.force):  
            self.run_command('test')
        _install.run(self)
        
# install the static libraries

class install_lib(_install_lib):
    def run(self):
        _install_lib.run(self)
        build_clib = self.get_finalized_command('build_clib')
        libs = build_clib.get_library_names()
        clib_dir = build_clib.build_clib
        for lib in libs:
            clib = 'lib' + lib + '.a'
            src_file = os.path.join(clib_dir, clib)
            dest_file = os.path.join(self.install_dir, 'pycbc', clib)
            file_util.copy_file(src_file, dest_file)


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
                 'install_lib':install_lib,
		         'test'  : test ,
		         'build_ext': build_ext},
    ext_modules = pycbc_extensions,
    libraries = pycbc_static_libraries,
    headers = pycbc_headers,
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

