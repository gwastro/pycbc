#!/usr/bin/python

"""
setup.py file for SWIG datavectorcpu
"""

from distutils.core import setup, Extension

vector_module = Extension('_datavectorcpu',
                          sources=['datavectorcpu_wrap.c','datavectorcpu.c'],
                          )


# base source files that do not require special libraries
datavectorcpu_swig_sources = ['datavectorcpu.i']
datavectorcpu_c_sources = ['datavectorcpu.c']


# define the extension module
datavectorcpu_ext = Extension( '_datavectorcpu', 
  sources = datavectorcpu_swig_sources + datavectorcpu_c_sources,
  depends = ['datavectorcpu.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'datavectorcpu',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap datavectorcpu""",
       ext_modules = [datavectorcpu_ext],
       py_modules = ["datavectorcpu"],
       )
