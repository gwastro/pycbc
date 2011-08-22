#!/usr/bin/python

"""
setup.py file for SWIG straindataopencl
"""

from distutils.core import setup, Extension

vector_module = Extension('_straindataopencl',
                          sources=['straindataopencl_wrap.c','straindataopencl.c'],
                          )


# base source files that do not require special libraries
straindataopencl_swig_sources = ['straindataopencl.i']
straindataopencl_c_sources = ['straindataopencl.c']


# define the extension module
straindataopencl_ext = Extension( '_straindataopencl', 
  sources = straindataopencl_swig_sources + straindataopencl_c_sources,
  depends = ['straindataopencl.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'straindataopencl',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap straindataopencl""",
       ext_modules = [straindataopencl_ext],
       py_modules = ["straindataopencl"],
       )
