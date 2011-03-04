import os
import sys
import glob
from types import *
from distutils import log
from distutils.core import setup, Extension
from distutils.command import config
from distutils.command import build
from distutils.command import build_clib
from distutils.command import install_lib
from distutils.command.clean import clean
from distutils import file_util
from distutils import spawn

set_runtime_library_dirs = False
argv_replace = []

for arg in sys.argv:
    if arg.startswith('--rpath'):
        set_runtime_library_dirs = True
    else:
        argv_replace.append(arg)

sys.argv = argv_replace

# package version number
ver = '0.1'

# classes specific to the pycbc build
class pycbc_clean(clean):
  def finalize_options (self):
    clean.finalize_options(self)
    self.clean_files = [ 'include/pycbc_wrap.c', 'include/pycbc.py' ]

  def run(self):
    clean.run(self)
    for f in self.clean_files:
      self.announce('removing ' + f)
      try:
        os.unlink(f)
      except:
        log.warn("'%s' does not exist -- can't clean it" % f)

class pycbc_config(config.config):
  def run(self):

    # check for required python modules
    try:
      log.info('checking for glue')
      import glue as glue_test
    except ImportError:
      log.error('glue module is missing')
      sys.exit(1)

    # variables to store the output of pkgconfig
    cflags = ''
    libs = ''

    # check for pkgconfig files for the three libraries that we need
    for pkg in ['fftw3f', 'gsl', 'lal', 'lalinspiral', 'cuda']:
      log.info('checking for %s package config file' % pkg )
      if os.system('pkg-config --exists %s 2>/dev/null' % pkg) == 0:
        log.info('found package config file for %s' % pkg )
        pkgcfg = os.popen('pkg-config --cflags ' + pkg)
        cflags += ' ' + pkgcfg.readline().strip()
        pkgcfg.close()
        pkgcfg = os.popen('pkg-config --libs ' + pkg)
        libs += ' ' + pkgcfg.readline().strip()
        pkgcfg.close()
      else:
        # assume that the user has fftw, gsl and cuda in a standard place
        if pkg is 'fftw3f':
          xlibs = ' -lfftw3f -lm'
          libs += xlibs
        if pkg is 'gsl':
          xlibs = ' -lgsl -lgslcblas -lm'
          libs += xlibs
        if pkg is 'cuda':
          xlibs = ' -L/usr/local/cuda/lib -lcublas -lcufft -lcudart -lm'
          xcflags = ' -I/usr/local/cuda/include'
          libs += xlibs
          cflags += xcflags
        if pkg is 'fftw3f' or pkg is 'gsl' or pkg is 'cuda':
          log.warn('no pkgconfig file found for %s, using %s' % (pkg, xlibs))
        else:
          log.warn('no pkgconfig file found for %s' % pkg)

    # parse the output of pkgconfig into python lists
    extra_cflags = [x for x in cflags.split() if x[0:2] != '-I']
    include_dirs = [x[2:] for x in cflags.split() if x[0:2] == '-I']
    library_dirs = [x[2:] for x in libs.split() if x[0:2] == '-L']
    libraries = [x[2:] for x in libs.split() if x[0:2] == '-l']

    # check that the required header files are present
    required_headers = [
      'fftw3.h', 
      'gsl/gsl_errno.h',
      'gsl/gsl_odeiv.h',
      'lal/LALVersion.h',
      'lal/FindChirp.h',
      'cuda.h',
      'cufft.h'
    ]
    for header in required_headers:
      log.info("checking for header %s" % header)
      if not self.check_header(header,include_dirs):
        log.error( 'header %s is missing' % header )
        sys.exit(1)

    # check that we can link against the required libraries
    required_functions = [
      'fftwf_execute_dft()',
      'gsl_strerror()',
      'LALVersion()',
      'XLALSumStrain()', # from lalinspiral/src/NRWaveInject.c
      'cufftPlan1d()'
    ]
    for func in required_functions:
      log.info('checking for function %s' % func)
      if not self.check_func(func, 
      libraries=libraries, library_dirs=library_dirs):
        log.error( 'function %s is missing' % func )
        sys.exit(1)

    # add everything to the build
    self.distribution.libraries[0][1]['include_dirs'] += include_dirs
    self.distribution.ext_modules[0].include_dirs += include_dirs
    self.distribution.ext_modules[0].extra_compile_args += extra_cflags
    self.distribution.ext_modules[0].library_dirs += library_dirs
    self.distribution.ext_modules[0].libraries += libraries

    # tell linker about nonstandard location of runtime libraries
    if set_runtime_library_dirs:
        self.distribution.ext_modules[0].runtime_library_dirs += library_dirs

class pycbc_build(build.build):
  def run(self):
    self.run_command('config')
    build.build.run(self)

  sub_commands = [
    ('build_clib',    build.build.has_c_libraries),
    ('build_ext',     build.build.has_ext_modules),
    ('build_py',      build.build.has_pure_modules),
    ('build_scripts', build.build.has_scripts),
    ]

class pycbc_build_clib(build_clib.build_clib):
  def run(self):
    # if we are using mac os x, remove any existing static libraries
    # so that otool does not complain
    if sys.platform[:6] == "darwin":
      libs = self.get_library_names()
      clib_dir = self.build_clib
      for lib in libs:
        clib = 'lib' + lib + '.a'
        src_file = os.path.join(clib_dir, clib)
        log.info("removing " + src_file)
        try:
          os.unlink(src_file)
        except:
          pass

    # now build the static libraries
    build_clib.build_clib.run(self)

class pycbc_install_lib(install_lib.install_lib):
  def run (self):
    install_lib.install_lib.run(self)

    # install any static c libraries
    if self.distribution.has_c_libraries():
      build_clib = self.get_finalized_command('build_clib')
      libs = build_clib.get_library_names()
      clib_dir = build_clib.build_clib
      for lib in libs:
        clib = 'lib' + lib + '.a'
        src_file = os.path.join(clib_dir, clib)
        dest_file = os.path.join(self.install_dir, clib)
        file_util.copy_file(src_file, dest_file)
        if sys.platform[:6] == "darwin":
          spawn.spawn(['ranlib'] + [dest_file])

# base source files that do not require special libraries
pycbc_swig_sources = [ 'include/pycbc.i' ]
pycbc_c_sources = [ 
  'src/types.c'
  ]

# define the extension module
pycbc_ext = Extension( '_pycbc', 
  sources = pycbc_swig_sources + pycbc_c_sources,
  depends = ['include/pycbc.h'],
  swig_opts = [],
  include_dirs = ['include'],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

# do the actual work of building the package
setup (name = 'pycbc',
  version = ver,
  description = 'inspiral analysis toolkit',
  author = 'Duncan Brown',
  author_email = 'dabrown@physics.syr.edu',
  cmdclass = {
    'config' : pycbc_config, 
    'build' : pycbc_build, 
    'build_clib' : pycbc_build_clib,
    'install_lib' : pycbc_install_lib, 
    'clean' : pycbc_clean},
  ext_modules = [pycbc_ext],
  libraries = [[ 'pycbc', { 
    'sources' : pycbc_c_sources,
    'include_dirs' : ['include'],
    'macros' : [] }]],
  headers = ['include/pycbc.h'],
  package_dir = {'' : 'include'},
  py_modules = ['pycbc'] )
