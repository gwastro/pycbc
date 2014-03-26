from pycbc.types import zeros
import numpy
import ctypes
import functools

def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer

#FFTW constants, these are pulled from fftw3.h
FFTW_FORWARD = -1
FFTW_BACKWARD = 1
FFTW_ESTIMATE =  1 << 0
FFTW_MEASURE   =  1 << 1
FFTW_EXHASTIVE = 1 << 2
FFTW_MEASRE = 0
FFTW_DESTROY_INPUT = 1 << 0
FFTW_NALIGNED = 1 << 1
FFTW_CONSERVE_MEMORY = 1 << 2
FFTW_EXHASTIVE = 1 << 3
FFTW_PRESERVE_INPUT = 1 << 4
FFTW_PATIENT = 1 << 5
FFTW_ESTIMATE = 1 << 6
FFTW_WISDOM_ONLY = 1 << 21

# load the single and double precision libraries
double_lib_name = 'libfftw3.so'
double_lib = ctypes.cdll.LoadLibrary(double_lib_name)
float_lib_name = 'libfftw3f.so'
float_lib = ctypes.cdll.LoadLibrary(float_lib_name)

# Create function maps for the dtypes
plan_function = {'float32': {'complex64': float_lib.fftwf_plan_dft_r2c_1d},
                 'float64': {'complex128': double_lib.fftw_plan_dft_r2c_1d},
                 'complex64': {'float32': float_lib.fftwf_plan_dft_c2r_1d,
                             'complex64': float_lib.fftwf_plan_dft_1d,
                              },
                 'complex128': {'float64': double_lib.fftw_plan_dft_c2r_1d,
                                'complex128': double_lib.fftw_plan_dft_1d,
                               }
                }

execute_function = {'float32': {'complex64': float_lib.fftwf_execute_dft_r2c},
                    'float64': {'complex128': double_lib.fftw_execute_dft_r2c},
                    'complex64': {'float32': float_lib.fftwf_execute_dft_c2r,
                                  'complex64': float_lib.fftwf_execute_dft,
                            },
                    'complex128': {'float64': double_lib.fftw_execute_dft_c2r,
                                   'complex128': double_lib.fftw_execute_dft,
                                  }
                   }

def alignment_of(vec):
    """ Retrn the byte alignment of this array
    """
    pointer = vec.data.ctypes.data
    return double_lib.fftw_alignment_of(pointer)

@memoize
def plan(size, idtype, odtype, direction, flags):
    if idtype.kind == 'c' and odtype.kind == 'f':
        size = 2*(size - 1)

    # make some representative arrays
    ip = zeros(size, dtype=idtype)
    op = zeros(size, dtype=odtype)

    # Get the plan function
    idtype = numpy.dtype(idtype)
    odtype = numpy.dtype(odtype)
    f = plan_function[str(idtype)][str(odtype)]
    
    # handle the C2C cases (forward and reverse)
    if idtype.kind == odtype.kind:
        f.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, 
                      ctypes.c_int, ctypes.c_int]
        theplan = f(size, ip.ptr, op.ptr, direction, flags)
    # handle the R2C and C2R case
    else:
        f.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, 
                      ctypes.c_int]
        theplan = f(size, ip.ptr, op.ptr, flags)  
         
    return theplan
    
def execute(plan, invec, outvec):
    f = execute_function[str(invec.dtype)][str(outvec.dtype)]
    f.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    f(plan, invec.ptr, outvec.ptr)
    
def fft(invec, outvec, prec, itype, otype):
    theplan = plan(len(invec), invec.dtype, outvec.dtype, FFTW_FORWARD, FFTW_MEASURE)
    execute(theplan, invec, outvec)
    
def ifft(invec, outvec, prec, itype, otype):
    theplan = plan(len(invec), invec.dtype, outvec.dtype, FFTW_BACKWARD, FFTW_MEASURE)
    execute(theplan, invec, outvec)

    
    
