import ctypes, functools, pycbc.libutils
from pycbc.types import zeros
from .core import _BaseFFT, _BaseIFFT
import pycbc.scheme as _scheme

def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer

lib = pycbc.libutils.get_ctypes_library('mkl_rt', [])
if lib is None:
    raise ImportError

#MKL constants  taken from mkl_df_defines.h
DFTI_FORWARD_DOMAIN = 0
DFTI_DIMENSION = 1
DFTI_LENGTHS = 2
DFTI_PRECISION = 3
DFTI_FORWARD_SCALE  = 4
DFTI_BACKWARD_SCALE = 5
DFTI_NUMBER_OF_TRANSFORMS = 7
DFTI_COMPLEX_STORAGE = 8
DFTI_REAL_STORAGE = 9
DFTI_CONJUGATE_EVEN_STORAGE = 10
DFTI_PLACEMENT = 11
DFTI_INPUT_STRIDES = 12
DFTI_OUTPUT_STRIDES = 13
DFTI_INPUT_DISTANCE = 14
DFTI_OUTPUT_DISTANCE = 15
DFTI_WORKSPACE = 17
DFTI_ORDERING = 18
DFTI_TRANSPOSE = 19
DFTI_DESCRIPTOR_NAME = 20
DFTI_PACKED_FORMAT = 21
DFTI_COMMIT_STATUS = 22
DFTI_VERSION = 23
DFTI_NUMBER_OF_USER_THREADS = 26
DFTI_THREAD_LIMIT = 27
DFTI_COMMITTED = 30
DFTI_UNCOMMITTED = 31
DFTI_COMPLEX = 32
DFTI_REAL = 33
DFTI_SINGLE = 35
DFTI_DOUBLE = 36
DFTI_COMPLEX_COMPLEX = 39
DFTI_COMPLEX_REAL = 40
DFTI_REAL_COMPLEX = 41
DFTI_REAL_REAL = 42
DFTI_INPLACE = 43         
DFTI_NOT_INPLACE = 44      
DFTI_ORDERED = 48
DFTI_BACKWARD_SCRAMBLED = 49
DFTI_ALLOW = 51            
DFTI_AVOID = 52
DFTI_NONE = 53
DFTI_CCS_FORMAT = 54       
DFTI_PACK_FORMAT = 55    
DFTI_PERM_FORMAT = 56      
DFTI_CCE_FORMAT = 57      

mkl_prec = {'single': DFTI_SINGLE,
            'double': DFTI_DOUBLE,
           }
            
mkl_domain = {'real': {'complex': DFTI_REAL},
              'complex': {'real': DFTI_REAL,
                          'complex':DFTI_COMPLEX,
                         }
             }

def check_status(status):
    """ Check the status of a mkl functions and raise a python exeption if 
    there is an error.
    """
    if status:
        msg = lib.DftiErrorMessage(status)
        msg = ctypes.c_char_p(msg).value
        raise RuntimeError(msg)
   
@memoize     
def create_descriptor(size, idtype, odtype, inplace):
    invec = zeros(1, dtype=idtype)
    outvec = zeros(1, dtype=odtype)
    
    desc = ctypes.c_void_p(1)
    f = lib.DftiCreateDescriptor
    f.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int]
    
    prec = mkl_prec[invec.precision]
    domain = mkl_domain[str(invec.kind)][str(outvec.kind)]
    
    status = f(ctypes.byref(desc), prec, domain, 1, size)
    if inplace:
        lib.DftiSetValue(desc, DFTI_PLACEMENT, DFTI_INPLACE)
    else:
        lib.DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    lib.DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_CCS_FORMAT)
    lib.DftiCommitDescriptor(desc)
    check_status(status)
    return desc
    
def fft(invec, outvec, prec, itype, otype):
    descr = create_descriptor(max(len(invec), len(outvec)), invec.dtype,
                              outvec.dtype, (invec.ptr == outvec.ptr))
    f = lib.DftiComputeForward
    f.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    status = f(descr, invec.ptr, outvec.ptr)
    check_status(status)

def ifft(invec, outvec, prec, itype, otype):
    descr = create_descriptor(max(len(invec), len(outvec)), invec.dtype,
                              outvec.dtype, (invec.ptr == outvec.ptr))
    f = lib.DftiComputeBackward
    f.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    status = f(descr, invec.ptr, outvec.ptr)
    check_status(status)
 
# Class based API

_create_descr = lib.DftiCreateDescriptor
_create_descr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int]

def _get_desc(fftobj):
    desc = ctypes.c_void_p(1)
    prec = mkl_prec[fftobj.invec.precision]
    domain = mkl_domain[str(fftobj.invec.kind)][str(fftobj.outvec.kind)]
    status = _create_descr(ctypes.byref(desc), prec, domain, 1, fftobj.size)
    check_status(status)
    # Now we set various things depending on exactly what kind of transform we're
    # performing.

    lib.DftiSetValue.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]

    # The following only matters if the transform is C2R or R2C
    status = lib.DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE,
                              DFTI_COMPLEX_COMPLEX)
    check_status(status)

    # In-place or out-of-place:
    if fftobj.inplace:
        status = lib.DftiSetValue(desc, DFTI_PLACEMENT, DFTI_INPLACE)
    else:
        status = lib.DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    check_status(status)

    # If we are performing a batched transform:
    if fftobj.nbatch > 1:
        status = lib.DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, fftobj.nbatch)
        check_status(status)
        status = lib.DftiSetValue(desc, DFTI_INPUT_DISTANCE, fftobj.idist)
        check_status(status)
        status = lib.DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, fftobj.odist)
        check_status(status)

    # Knowing how many threads will be allowed may help select a better transform
    nthreads = _scheme.mgr.state.num_threads
    status = lib.DftiSetValue(desc, DFTI_THREAD_LIMIT, nthreads)
    check_status(status)

    # Now everything's ready, so commit
    status = lib.DftiCommitDescriptor(desc)
    check_status(status)

    return desc

class FFT(_BaseFFT):
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(FFT, self).__init__(invec, outvec, nbatch, size)
        self.iptr = self.invec.ptr
        self.optr = self.outvec.ptr
        self._efunc = lib.DftiComputeForward
        self._efunc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
        self.desc = _get_desc(self)

    def execute(self):
        self._efunc(self.desc, self.iptr, self.optr)

class IFFT(_BaseIFFT):
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(IFFT, self).__init__(invec, outvec, nbatch, size)
        self.iptr = self.invec.ptr
        self.optr = self.outvec.ptr
        self._efunc = lib.DftiComputeBackward
        self._efunc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
        self.desc = _get_desc(self)

    def execute(self):
        self._efunc(self.desc, self.iptr, self.optr)
