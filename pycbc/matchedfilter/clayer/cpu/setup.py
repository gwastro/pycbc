#!/usr/bin/python

"""
setup.py file for SWIG matchedfiltercpu
"""

from distutils.core import setup, Extension

vector_module = Extension('_matchedfiltercpu',
                          sources=['matchedfiltercpu_wrap.c','matchedfiltercpu.c'],
                          )


# base source files that do not require special libraries
matchedfiltercpu_swig_sources = ['matchedfiltercpu.i']
matchedfiltercpu_c_sources = ['matchedfiltercpu.c']


# define the extension module
matchedfiltercpu_ext = Extension( '_matchedfiltercpu', 
  sources = matchedfiltercpu_swig_sources + matchedfiltercpu_c_sources,
  depends = ['matchedfiltercpu.h'],
  swig_opts = [],
  include_dirs = [],
  extra_compile_args = ['-Wall'],
  library_dirs = [],
  libraries = [])

setup (name = 'matchedfiltercpu',
       version = '0.1',
       author = "Karsten Wiesner",
       description = """Swig wrap matchedfiltercpu""",
       ext_modules = [matchedfiltercpu_ext],
       py_modules = ["matchedfiltercpu"],
       )
