
# Check for optional components of the PyCBC Package
try:
    import pycuda
    import pycuda.gpuarray
    have_cuda=True
except ImportError:
    have_cuda=False
    
try:
    import pyopencl
    import pyopencl.array
    have_opencl=True
except ImportError:
    have_opencl=False
