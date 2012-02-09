
# Check for optional components of the PyCBC Package
have_cuda = False
have_opencl = False

try:
    import pycuda as _pycuda
    have_cuda=True
except ImportError:
    pass
    
try:
    import pyopencl as _pyopencl
    have_opencl=True
except ImportError:
    pass
