#!/usr/bin/python

"""
setup.py file for SWIG matchedfilteropencl
"""

from distutils.core import setup, Extension

vector_module = Extension('_matchedfilteropencl',
                          sources=['matchedfilteropencl_wrap.c','matchedfilteropencl.c'],
                          )


# base source files that do not require special libraries
matchedfilteropencl_swig_sources = ['matchedfilteropencl.i']
matchedfilteropencl_c_sources = ['matchedfilteropencl.c']


# define the extension module
matchedfilteropencl_ext = Extension( '_matchedfilteropencl', 
  sources = matchedfilteropencl_swig_sources + matchedfilteropencl_c_sources,
  depends = ['matchedfilteropencl.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'matchedfilteropencl',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap matchedfilteropencl""",
       ext_modules = [matchedfilteropencl_ext],
       py_modules = ["matchedfilteropencl"],
       )
