
# Check for optional components of the PyCBC Package
try:
    import pycuda as _pycuda
    have_cuda=True
except ImportError:
    have_cuda=False
    
try:
    import pyopencl as _pyopencl
    have_opencl=True
except ImportError:
    have_opencl=False
