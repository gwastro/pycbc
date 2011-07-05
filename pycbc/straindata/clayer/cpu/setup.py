#!/usr/bin/python

"""
setup.py file for SWIG straindatacpu
"""

from distutils.core import setup, Extension

vector_module = Extension('_straindatacpu',
                          sources=['straindatacpu_wrap.c','straindatacpu.c'],
                          )


# base source files that do not require special libraries
straindatacpu_swig_sources = ['straindatacpu.i']
straindatacpu_c_sources = ['straindatacpu.c']


# define the extension module
straindatacpu_ext = Extension( '_straindatacpu', 
  sources = straindatacpu_swig_sources + straindatacpu_c_sources,
  depends = ['straindatacpu.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'straindatacpu',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap straindatacpu""",
       ext_modules = [straindatacpu_ext],
       py_modules = ["straindatacpu"],
       )
