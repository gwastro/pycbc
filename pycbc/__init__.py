"""PyCBC contains a toolkit for CBC gravitational wave analysis

"""
# Check for optional components of the PyCBC Package

try:
    import pycuda as _pycuda
    HAVE_CUDA=True
except ImportError:
    HAVE_CUDA=False
    
try:
    import pyopencl as _pyopencl
    HAVE_OPENCL=True
except ImportError:
    HAVE_OPENCL=False
