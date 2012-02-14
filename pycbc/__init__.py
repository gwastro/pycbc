"""

"""
# Check for optional components of the PyCBC Package
HAVE_CUDA = False
HAVE_OPENCL = False

try:
    import pycuda as _pycuda
    HAVE_CUDA=True
except ImportError:
    pass
    
try:
    import pyopencl as _pyopencl
    HAVE_OPENCL=True
except ImportError:
    pass
