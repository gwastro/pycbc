#!/usr/bin/python

"""
setup.py file for SWIG mymodule types example
"""

from distutils.core import setup, Extension
from numpy import get_include

testlal_module = Extension('_testlal',
                           sources=['testlal_wrap.c'],
                           include_dirs = ['/home/jwillis/dev/root/lalsuite/branches/master/include',get_include()],
                           library_dirs = ['/home/jwillis/dev/root/lalsuite/branches/master/lib'],
                           define_macros=[('SWIG_TYPE_TABLE','swiglal')],
                           extra_compile_args=['-std=c99'],
                           libraries=['lal','fftw3','fftw3f','gsl']
                          )

setup (name = 'testlal',
       version = '0.0',
       author = "Josh Willis",
       description = """Test SWIG wrapping with numpy interface to LAL""",
       ext_modules = [testlal_module],
       py_modules = ["testlal"],
       )
