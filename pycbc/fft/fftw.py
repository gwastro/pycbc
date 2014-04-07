from pycbc.types import zeros
import numpy
import ctypes
import functools

# IMPORTANT NOTE TO PYCBC DEVELOPERS:
# Because this module is loaded automatically when present, and because
# no FFTW function should be called until the user has had the chance
# to set the threading backend, it is ESSENTIAL that simply loading this
# module should not actually *call* ANY functions.

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

FFTW_MEASURE = 0
FFTW_DESTROY_INPUT = 1 << 0
FFTW_UNALIGNED = 1 << 1
FFTW_CONSERVE_MEMORY = 1 << 2
FFTW_EXHAUSTIVE = 1 << 3
FFTW_PRESERVE_INPUT = 1 << 4
FFTW_PATIENT = 1 << 5
FFTW_ESTIMATE = 1 << 6
FFTW_WISDOM_ONLY = 1 << 21

# Load the single and double precision libraries
# We need to construct them directly with CDLL so
# we can give the RTLD_GLOBAL mode, which we must do
# in order to use the threaded libraries as well.
double_lib_name = 'libfftw3.so'
double_lib = ctypes.CDLL(double_lib_name,mode=ctypes.RTLD_GLOBAL)
float_lib_name = 'libfftw3f.so'
float_lib = ctypes.CDLL(float_lib_name,mode=ctypes.RTLD_GLOBAL)

# Even if we don't have threaded support, we define a number
# of threads.  But it's 0 without threaded support (which is 
# an invalid argument) which allows our planning wrapper
# to distinguish threaded from nonthreaded plans
_fftw_nthreads = 0
def get_nthreads():
    return _fftw_nthreads

_fftw_threaded_lib = None
_fftw_threaded_set = False

HAVE_FFTW_THREADED = False
def set_threads_backend(backend):
    global _fftw_threaded_set
    global _fftw_threaded_lib
    if _fftw_threaded_set:
        raise RuntimeError(
            "Threading backend for FFTW already set to {0}; cannot be changed".format(_fftw_threaded_lib))
    if backend == 'pthreads':
        try:
            double_threaded_lib = ctypes.CDLL('libfftw3_threads.so',mode=ctypes.RTLD_GLOBAL)
            float_threaded_lib = ctypes.CDLL('libfftw3f_threads.so',mode=ctypes.RTLD_GLOBAL)
            HAVE_FFTW_THREADED = True
            _fftw_threaded_lib = 'pthreads'
            _fftw_threaded_set = True
        except:
            raise RuntimeError("Could not load 'pthreads' backend")
    elif backend == 'openmp':
        try:
            double_threaded_lib = ctypes.CDLL('libfftw3_omp.so',mode=ctypes.RTLD_GLOBAL)
            float_threaded_lib = ctypes.CDLL('libfftw3f_omp.so',mode=ctypes.RTLD_GLOBAL)
            HAVE_FFTW_THREADED = True
            _fftw_threaded_lib = 'openmp'
            _fftw_threaded_set = True
        except:
            raise RuntimeError("Could not load 'openmp' backend")
    # We only use the following internally; users should just not call this func when
    # not using threads
    elif backend == 'unthreaded':
        _fftw_threaded_lib = 'unthreaded'
        _fftw_threaded_set = True
    else:
        raise ValueError("Invalid input to set_threads_backend(): must be 'pthreads' or 'openmp'")
    # Call each of these exactly once, before anything else...
    if backend != 'unthreaded':
        dret = double_threaded_lib.fftw_init_threads()
        fret = float_threaded_lib.fftwf_init_threads()
        # FFTW for some reason uses *0* to indicate failure.  In C.
        if (dret == 0) or (fret == 0):
            raise RuntimeError("Threaded FFTW found and loaded, but could not initialize")

if HAVE_FFTW_THREADED:
    # Now a function to use a given number of threads
    def use_nthreads(nthreads):
        """
        Set the current number of threads used in FFTW planning/
        execution.  Must be an non-negative integer.
        """
        global _fftw_nthreads
        if not (isinstance(nthreads,int) and (nthreads>0)):
            raise ValueError("nthreads must be nonegative integer")
        _fftw_nthreads = nthreads

        dplanwthr = double_threaded_lib.fftw_plan_with_nthreads
        fplanwthr = float_threaded_lib.fftwf_plan_with_nthreads
        dplanwthr.restype = None
        fplanwthr.restype = None
        dplanwthr(nthreads)
        fplanwthr(nthreads)

    # We need to initialize to something valid
    use_nthreads(1)

# Function to import system-wide wisdom files.

def import_sys_wisdom():
    if not _fftw_threaded_set:
        set_threads_backend('unthreaded')
    double_lib.fftw_import_system_wisdom()
    float_lib.fftwf_import_system_wisdom()
    
# We provide an interface for changing the "measure level"
# By default 1, which does some but not much planning,
# but we provide functions to read and set it

_default_measurelvl = 1
def get_measure_level():
    """
    Get the current 'measure level' used in deciding how much effort to put into
    creating FFTW plans.  From least effort (and shortest planning time) to most
    they are 0 to 3.  No arguments.
    """
    return _default_measurelvl

def set_measure_level(mlvl):
    """
    Set the current 'measure level' used in deciding how much effort to expend
    creating FFTW plans.  Must be an integer from 0 (least effort, shortest time)
    to 3 (most effort and time).
    """

    global _default_measurelvl
    if mlvl not in (0,1,2,3):
        raise ValueError("Measure level can only be one of 0, 1, 2, or 3")
    _default_measurelvl = mlvl

_flag_dict = {0: FFTW_ESTIMATE,
              1: FFTW_MEASURE,
              2: FFTW_MEASURE|FFTW_PATIENT,
              3: FFTW_MEASURE|FFTW_PATIENT|FFTW_EXHAUSTIVE}
def get_flag(mlvl,aligned):
    if aligned:
        return _flag_dict[mlvl]
    else:
        return (_flag_dict[mlvl]|FFTW_UNALIGNED)

# Add the ability to read/store wisdom to filenames

def import_single_wisdom_from_filename(filename):
    if not _fftw_threaded_set:
        set_threads_backend('unthreaded')
    f = float_lib.fftwf_import_wisdom_from_filename
    f.argtypes = [ctypes.c_char_p]
    retval = f(filename)
    if retval == 0:
        raise RuntimeError("Could not import wisdom from file {0}".format(filename))

def import_double_wisdom_from_filename(filename):
    if not _fftw_threaded_set:
        set_threads_backend('unthreaded')
    f = double_lib.fftw_import_wisdom_from_filename
    f.argtypes = [ctypes.c_char_p]
    retval = f(filename)
    if retval == 0:
        raise RuntimeError("Could not import wisdom from file {0}".format(filename))

def export_single_wisdom_to_filename(filename):
    if not _fftw_threaded_set:
        set_threads_backend('unthreaded')
    f = float_lib.fftwf_export_wisdom_to_filename
    f.argtypes = [ctypes.c_char_p]
    retval = f(filename)
    if retval == 0:
        raise RuntimeError("Could not export wisdom to file {0}".format(filename))

def export_double_wisdom_to_filename(filename):
    if not _fftw_threaded_set:
        set_threads_backend('unthreaded')
    f = double_lib.fftw_export_wisdom_to_filename
    f.argtypes = [ctypes.c_char_p]
    retval = f(filename)
    if retval == 0:
        raise RuntimeError("Could not export wisdom to file {0}".format(filename))

# Create function maps for the dtypes
plan_function = {'float32': {'complex64': float_lib.fftwf_plan_dft_r2c_1d},
                 'float64': {'complex128': double_lib.fftw_plan_dft_r2c_1d},
                 'complex64': {'float32': float_lib.fftwf_plan_dft_c2r_1d,
                               'complex64': float_lib.fftwf_plan_dft_1d},
                 'complex128': {'float64': double_lib.fftw_plan_dft_c2r_1d,
                                'complex128': double_lib.fftw_plan_dft_1d}
                }

execute_function = {'float32': {'complex64': float_lib.fftwf_execute_dft_r2c},
                    'float64': {'complex128': double_lib.fftw_execute_dft_r2c},
                    'complex64': {'float32': float_lib.fftwf_execute_dft_c2r,
                                  'complex64': float_lib.fftwf_execute_dft},
                    'complex128': {'float64': double_lib.fftw_execute_dft_c2r,
                                   'complex128': double_lib.fftw_execute_dft}
                   }

@memoize
def plan(size, idtype, odtype, direction, mlvl, nthreads, aligned):
    if not _fftw_threaded_set:
        set_threads_backend('unthreaded')
    # Convert a measure-level to flags
    flags = get_flag(mlvl,aligned)

    if (idtype == odtype):
        # We're in the complex-to-complex case, so lengths are the same
        isize = size
        osize = size
    elif (idtype.kind == 'c') and (odtype.kind == 'f'):
        # Complex-to-real (reverse), so size is length of real array
        isize = size/2+1
        osize = size
    else:
        # Real-to-complex (forward), and size is still that of real
        isize = size
        osize = size/2+1

    # make some representative arrays
    ip = zeros(isize, dtype=idtype)
    op = zeros(osize, dtype=odtype)

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


# Note that we don't need to check whether we've set the threading backend
# in the following functions, since execute is not called directly and
# the fft and ifft will call plan first.    
def execute(plan, invec, outvec):
    f = execute_function[str(invec.dtype)][str(outvec.dtype)]
    f.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    f(plan, invec.ptr, outvec.ptr)
    
def fft(invec, outvec, prec, itype, otype):
    theplan = plan(len(invec), invec.dtype, outvec.dtype, FFTW_FORWARD,
                   get_measure_level(),get_nthreads(),
                   (invec._data.isaligned and outvec._data.isaligned))
    execute(theplan, invec, outvec)
    
def ifft(invec, outvec, prec, itype, otype):
    theplan = plan(len(outvec), invec.dtype, outvec.dtype, FFTW_BACKWARD,
                   get_measure_level(),get_nthreads(),
                   (invec._data.isaligned and outvec._data.isaligned))
    execute(theplan, invec, outvec)

    
def insert_fft_options(optgroup):
    """
    Inserts the options that affect the behavior of this backend

    Parameters
    ----------
    optgroup: fft_option
       OptionParser argument group whose options are extended
    """
    optgroup.add_argument("--fftw-measure-level", 
                      help="Determines the measure level used in planning "
                           "FFTW FFTs; allowed values are: " + str([0,1,2,3]), 
                      type=int, default=_default_measurelvl)
    optgroup.add_argument("--fftw-threads-backend", 
                      help="Give 'pthreads' or 'openmp' to specify which threaded FFTW to use",
                      default=None)
    optgroup.add_argument("--fftw-use-nthreads", 
                      help="Number of threads to use in FFT planning/execution",
                      type=int, default=0) # 0 is sentinel for no thread support
    optgroup.add_argument("--fftw-input-float-wisdom-file", 
                      help="Filename from which to read single-precision wisdom",
                      default=None)
    optgroup.add_argument("--fftw-input-double-wisdom-file", 
                      help="Filename from which to read double-precision wisdom",
                      default=None)
    optgroup.add_argument("--fftw-output-float-wisdom-file", 
                      help="Filename to which to write single-precision wisdom",
                      default=None)
    optgroup.add_argument("--fftw-output-double-wisdom-file", 
                      help="Filename to which to write double-precision wisdom",
                      default=None)

def verify_fft_options(opt,parser):
    """Parses the FFT options and verifies that they are 
       reasonable. 
         
    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes.
    parser : object
        OptionParser instance.
    """
    if opt.fftw_measure_level not in [0,1,2,3]:
        parser.error("{0} is not a valid FFTW measure level.".format(opt.fftw_measure_level))
    if (opt.fftw_use_nthreads > 1) and (opt.fftw_threads_backend is None):
        parser.error("You specified a number of threads, but not a threads backend (pthreads or openmp)")
    if (opt.fftw_threads_backend is not None):
        if (opt.fftw_threads_backend not in ['openmp','pthreads']):
            parser.error("Invalid threads backend; must be 'pthreads' or 'openmp'")

def from_cli(opt):
    # If we specified a threading backend on the command line, we should set
    # that first.
    set_threads_backend(opt.fftw_threads_backend)
    set_measure_level(opt.fftw_measure_level)
    if (opt.fftw_use_nthreads > 1):
        use_nthreads(opt.fftw_use_nthreads)
    # We don't go ahead and import/export wisdom, because that depends on 
    # plan creation.  Instead just return those in the kwd returns
    kwdrets = {}
    if opt.fftw_input_float_wisdom_file is not None:
        kwdrets.update({"input_float_wisdom_file":opt.fftw_input_float_wisdom_file})
    if opt.fftw_input_double_wisdom_file is not None:
        kwdrets.update({"input_double_wisdom_file":opt.fftw_input_double_wisdom_file})
    if opt.fftw_output_float_wisdom_file is not None:
        kwdrets.update({"output_float_wisdom_file":opt.fftw_output_float_wisdom_file})
    if opt.fftw_output_double_wisdom_file is not None:
        kwdrets.update({"output_double_wisdom_file":opt.fftw_output_double_wisdom_file})
    return kwdrets
