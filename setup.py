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

ver = 0.1

# classes specific to the pycbc build
class pycbc_clean(clean):
    def finalize_options (self):
        clean.finalize_options(self)
        self.clean_files = [ 'build/lib*', 'build/temp*' ]
    def run(self):
        clean.run(self)
        for f in self.clean_files:
            self.announce('removing ' + f)
        try:
            os.unlink(f)
        except:
            log.warn("'%s' does not exist -- can't clean it" % f)

datavectorcpu_ext = Extension(
    name= '_datavectorcpu', 
    sources = ['pycbc/datavector/clayer/cpu/datavectorcpu.i','pycbc/datavector/clayer/cpu/datavectorcpu.c'],
    #  not in doc swig_opts = [],
    include_dirs = ['pycbc/datavector/clayer/cpu/','pycbc/datavector/clayer/'],
    define_macros=[('TESTMACRO', '1')],
    undef_macros=['TESTMACRO'],
    library_dirs = [],
    libraries = [],
    runtime_library_dirs = [],
    extra_objects = [],
    extra_compile_args = ['-Wall'],
    extra_link_args = [],
    export_symbols = [],
    # later ... depends = ['pycbc/datavector/clayer/cpu/'],
    language = []
)

matchedfiltercpu_ext = Extension( 
    name = '_matchedfiltercpu', 
    sources = ['pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.i','pycbc/matchedfilter/clayer/cpu/matchedfiltercpu.c'],
    # not in doc swig_opts = [],
    include_dirs = ['pycbc/matchedfilter/clayer/cpu/','pycbc/matchedfilter/clayer/'],
    define_macros=[('TESTMACRO', '1')],
    undef_macros=['TESTMACRO'],
    library_dirs = [],
    libraries = [],
    runtime_library_dirs = [],
    extra_objects = [],
    extra_compile_args = ['-Wall'],
    extra_link_args = [],
    export_symbols = [],
    #later ... depends = ['pycbc/matchedfilter/clayer/cpu/'],
    language = []
)


# do the actual work of building the package
setup (
    name = 'pycbc',
    version = ver,
    description = 'inspiral analysis toolkit',
    author = 'Ligo Virgo Collaboration - pyCBC team',
    author_email = 'https://sugwg-git.phy.syr.edu/dokuwiki/doku.php?id=pycbc:home',
    ext_modules = [datavectorcpu_ext, matchedfiltercpu_ext],
    py_modules = ['datavectorcpu','matchedfiltercpu'],
    # library_dirs = [],
    # libraries = [[ 'pycbc', {        #### what is this?:
    #   'sources' : pycbc_c_sources,
    #   'include_dirs' : ['include'],
    #   'macros' : [] }]],
    # headers = ['include/pycbc.h'],
    package_dir = {'pycbc' : 'pycbc'},
    packages = ['pycbc',
                'pycbc.datavector',
                'pycbc.matchedfilter']
)
