#!/usr/bin/python

"""
setup.py file for SWIG templatebankcpu
"""

from distutils.core import setup, Extension

vector_module = Extension('_templatebankcpu',
                          sources=['templatebankcpu_wrap.c','templatebankcpu.c'],
                          )


# base source files that do not require special libraries
templatebankcpu_swig_sources = ['templatebankcpu.i']
templatebankcpu_c_sources = ['templatebankcpu.c']


# define the extension module
templatebankcpu_ext = Extension( '_templatebankcpu', 
  sources = templatebankcpu_swig_sources + templatebankcpu_c_sources,
  depends = ['templatebankcpu.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'templatebankcpu',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap templatebankcpu""",
       ext_modules = [templatebankcpu_ext],
       py_modules = ["templatebankcpu"],
       )
