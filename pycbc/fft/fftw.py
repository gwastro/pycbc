from pycbc.types import zeros, complex64, complex128
import numpy as _np
import ctypes
import functools
import pycbc.scheme as _scheme
from pycbc.libutils import get_ctypes_library
from .core import _BaseFFT, _BaseIFFT

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
double_lib = get_ctypes_library('fftw3',['fftw3'],mode=ctypes.RTLD_GLOBAL)
float_lib = get_ctypes_library('fftw3f',['fftw3f'],mode=ctypes.RTLD_GLOBAL)
if (double_lib is None) or (float_lib is None):
    raise ImportError("Unable to find FFTW libraries")

# Support for FFTW's two different threading backends
_fftw_threaded_lib = None
_fftw_threaded_set = False
_double_threaded_lib = None
_float_threaded_lib = None

HAVE_FFTW_THREADED = False

# Although we set the number of threads based on the scheme,
# we need a private variable that records the last value used so
# we know whether we need to call plan_with_nthreads() again.
_fftw_current_nthreads = 0

# This function sets the number of threads used internally by FFTW
# in planning. It just takes a number of threads, rather than itself
# looking at scheme.mgr.num_threads, because it should not be called
# directly, but only by functions that get the value they use from
# scheme.mgr.num_threads

def _fftw_plan_with_nthreads(nthreads):
    global _fftw_current_nthreads
    if not HAVE_FFTW_THREADED:
        if (nthreads > 1):
                raise ValueError("Threading is NOT enabled, but {0} > 1 threads specified".format(nthreads))
        else:
            _pycbc_current_threads = nthreads
    else:
        dplanwthr = _double_threaded_lib.fftw_plan_with_nthreads
        fplanwthr = _float_threaded_lib.fftwf_plan_with_nthreads
        dplanwthr.restype = None
        fplanwthr.restype = None
        dplanwthr(nthreads)
        fplanwthr(nthreads)
        _fftw_current_nthreads = nthreads

# This is a global dict-of-dicts used when initializing threads and
# setting the threading library

_fftw_threading_libnames = { 'unthreaded' : {'double' : None, 'float' : None},
                             'openmp' : {'double' : 'fftw3_omp', 'float' : 'fftw3f_omp'},
                             'pthreads' : {'double' : 'fftw3_threads', 'float' : 'fftw3f_threads'}}

def _init_threads(backend):
    # This function actually sets the backend and initializes. It returns zero on
    # success and 1 if given a valid backend but that cannot be loaded.  It raises
    # an exception if called after the threading backend has already been set, or
    # if given an invalid backend.
    global _fftw_threaded_set
    global _fftw_threaded_lib
    global HAVE_FFTW_THREADED
    global _double_threaded_lib
    global _float_threaded_lib
    if _fftw_threaded_set:
        raise RuntimeError(
            "Threading backend for FFTW already set to {0}; cannot be changed".format(_fftw_threaded_lib))
    try:
        double_threaded_libname = _fftw_threading_libnames[backend]['double']
        float_threaded_libname =  _fftw_threading_libnames[backend]['float']
    except KeyError:
        raise ValueError("Backend {0} for FFTW threading does not exist!".format(backend))
    if double_threaded_libname is not None:
        try:
            # Note that the threaded libraries don't have their own pkg-config files;
            # we must look for them wherever we look for double or single FFTW itself
            _double_threaded_lib = get_ctypes_library(double_threaded_libname,['fftw3'],mode=ctypes.RTLD_GLOBAL)
            _float_threaded_lib =  get_ctypes_library(float_threaded_libname,['fftw3f'],mode=ctypes.RTLD_GLOBAL)
            if (_double_threaded_lib is None) or (_float_threaded_lib is None):
                raise RuntimeError("Unable to load threaded libraries {0} or {1}".format(double_threaded_libname,
                                                                                         float_threaded_libname))
            dret = _double_threaded_lib.fftw_init_threads()
            fret = _float_threaded_lib.fftwf_init_threads()
            # FFTW for some reason uses *0* to indicate failure.  In C.
            if (dret == 0) or (fret == 0):
                return 1
            HAVE_FFTW_THREADED = True
            _fftw_threaded_set = True
            _fftw_threaded_lib = backend
            return 0
        except:
            return 1
    else:
        # We get here when we were given the 'unthreaded' backend
        HAVE_FFTW_THREADED = False
        _fftw_threaded_set = True
        _fftw_threaded_lib = backend
        return 0

def set_threads_backend(backend=None):
    # This is the user facing function.  If given a backend it just
    # calls _init_threads and lets it do the work.  If not (the default)
    # then it cycles in order through threaded backends,
    global _fftw_threaded_set
    if backend is not None:
        retval = _init_threads(backend)
        # Since the user specified this backend raise an exception if the above failed
        if retval != 0:
            raise RuntimeError("Could not initialize FFTW threading backend {0}".format(backend))
    else:
        # Note that we pop() from the end, so 'openmp' is the first thing tried
        _backend_list = ['unthreaded','pthreads','openmp']
        while not _fftw_threaded_set:
            _next_backend = _backend_list.pop()
            retval = _init_threads(_next_backend)

# Function to import system-wide wisdom files.

def import_sys_wisdom():
    if not _fftw_threaded_set:
        set_threads_backend()
    double_lib.fftw_import_system_wisdom()
    float_lib.fftwf_import_system_wisdom()

# We provide an interface for changing the "measure level"
# By default this is 0, which does no planning,
# but we provide functions to read and set it

_default_measurelvl = 0
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
        set_threads_backend()
    f = float_lib.fftwf_import_wisdom_from_filename
    f.argtypes = [ctypes.c_char_p]
    retval = f(filename)
    if retval == 0:
        raise RuntimeError("Could not import wisdom from file {0}".format(filename))

def import_double_wisdom_from_filename(filename):
    if not _fftw_threaded_set:
        set_threads_backend()
    f = double_lib.fftw_import_wisdom_from_filename
    f.argtypes = [ctypes.c_char_p]
    retval = f(filename)
    if retval == 0:
        raise RuntimeError("Could not import wisdom from file {0}".format(filename))

def export_single_wisdom_to_filename(filename):
    if not _fftw_threaded_set:
        set_threads_backend()
    f = float_lib.fftwf_export_wisdom_to_filename
    f.argtypes = [ctypes.c_char_p]
    retval = f(filename)
    if retval == 0:
        raise RuntimeError("Could not export wisdom to file {0}".format(filename))

def export_double_wisdom_to_filename(filename):
    if not _fftw_threaded_set:
        set_threads_backend()
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
def plan(size, idtype, odtype, direction, mlvl, aligned, nthreads, inplace):
    if not _fftw_threaded_set:
        set_threads_backend()
    if nthreads != _fftw_current_nthreads:
        _fftw_plan_with_nthreads(nthreads)
    # Convert a measure-level to flags
    flags = get_flag(mlvl,aligned)
    
    # We make our arrays of the necessary type and size.  Things can be
    # tricky, especially for in-place transforms with one of input or
    # output real.
    if (idtype == odtype):
        # We're in the complex-to-complex case, so lengths are the same
        ip = zeros(size, dtype=idtype)
        if inplace:
            op = ip
        else:
            op = zeros(size, dtype=odtype)
    elif (idtype.kind == 'c') and (odtype.kind == 'f'):
        # Complex-to-real (reverse), so size is length of real array.
        # However the complex array may be larger (in bytes) and
        # should therefore be allocated first and reused for an in-place
        # transform
        ip = zeros(size/2+1, dtype=idtype)
        if inplace:
            op = ip.view(dtype=odtype)[0:size]
        else:
            op = zeros(size, dtype=odtype)
    else:
        # Real-to-complex (forward), and size is still that of real.
        # However it is still true that the complex array may be larger
        # (in bytes) and should therefore be allocated first and reused
        # for an in-place transform
        op = zeros(size/2+1, dtype=odtype)
        if inplace:
            ip = op.view(dtype=idtype)[0:size]
        else:
            ip = zeros(size, dtype=idtype)

    # Get the plan function
    idtype = _np.dtype(idtype)
    odtype = _np.dtype(odtype)
    f = plan_function[str(idtype)][str(odtype)]
    f.restype = ctypes.c_void_p

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

    # We don't need ip or op anymore
    del ip, op

    # And done...
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
                   get_measure_level(),(invec._data.isaligned and outvec._data.isaligned),
                   _scheme.mgr.state.num_threads, (invec.ptr == outvec.ptr))
    execute(theplan, invec, outvec)

def ifft(invec, outvec, prec, itype, otype):
    theplan = plan(len(outvec), invec.dtype, outvec.dtype, FFTW_BACKWARD,
                   get_measure_level(),(invec._data.isaligned and outvec._data.isaligned),
                   _scheme.mgr.state.num_threads, (invec.ptr == outvec.ptr))
    execute(theplan, invec, outvec)

# Class based API

# First, set up a lot of different ctypes functions:
plan_many_c2c_f = float_lib.fftwf_plan_many_dft
plan_many_c2c_f.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_int, ctypes.c_uint]
plan_many_c2c_f.restype = ctypes.c_void_p

plan_many_c2c_d = double_lib.fftw_plan_many_dft
plan_many_c2c_d.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_int, ctypes.c_uint]
plan_many_c2c_d.restype = ctypes.c_void_p

plan_many_c2r_f = float_lib.fftwf_plan_many_dft_c2r
plan_many_c2r_f.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_uint]
plan_many_c2r_f.restype = ctypes.c_void_p

plan_many_c2r_d = double_lib.fftw_plan_many_dft_c2r
plan_many_c2r_d.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_uint]
plan_many_c2r_d.restype = ctypes.c_void_p

plan_many_r2c_f = float_lib.fftwf_plan_many_dft_r2c
plan_many_r2c_f.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_uint]
plan_many_r2c_f.restype = ctypes.c_void_p

plan_many_r2c_d = double_lib.fftw_plan_many_dft_r2c
plan_many_r2c_d.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_uint]
plan_many_r2c_d.restype = ctypes.c_void_p

# Now set up a dictionary indexed by (str(input_dtype), str(output_dtype)) to
# translate input and output dtypes into the correct planning function.

_plan_funcs_dict = { ('complex64', 'complex64') : plan_many_c2c_f,
                     ('complex64', 'float32') : plan_many_r2c_f,
                     ('float32', 'complex64') : plan_many_c2r_f,
                     ('complex128', 'complex128') : plan_many_c2c_d,
                     ('complex128', 'float64') : plan_many_r2c_d,
                     ('float64', 'complex128') : plan_many_c2r_d }

# To avoid multiple-inheritance, we set up a function that returns much
# of the initialization that will need to be handled in __init__ of both
# classes.

def _fftw_setup(fftobj):
        n = _np.asarray([fftobj.size], dtype=_np.int32)
        inembed = _np.asarray([len(fftobj.invec)], dtype=_np.int32)
        onembed = _np.asarray([len(fftobj.outvec)], dtype=_np.int32)
        nthreads = _scheme.mgr.state.num_threads
        if not _fftw_threaded_set:
            set_threads_backend()
        if nthreads != _fftw_current_nthreads:
            _fftw_plan_with_nthreads(nthreads)  
        mlvl = get_measure_level()
        aligned = fftobj.invec.data.isaligned and fftobj.outvec.data.isaligned
        flags = get_flag(mlvl, aligned)
        plan_func = _plan_funcs_dict[ (str(fftobj.invec.dtype), str(fftobj.outvec.dtype)) ]
        tmpin = zeros(len(fftobj.invec), dtype = fftobj.invec.dtype)
        tmpout = zeros(len(fftobj.outvec), dtype = fftobj.outvec.dtype)
        # C2C, forward
        if fftobj.forward and (fftobj.outvec.dtype in [complex64, complex128]):
            plan = plan_func(1, n.ctypes.data, fftobj.nbatch,
                             tmpin.ptr, inembed.ctypes.data, 1, fftobj.idist,
                             tmpout.ptr, onembed.ctypes.data, 1, fftobj.odist,
                             FFTW_FORWARD, flags)
        # C2C, backward
        elif not fftobj.forward and (fftobj.invec.dtype in [complex64, complex128]):
            plan = plan_func(1, n.ctypes.data, fftobj.nbatch,
                             tmpin.ptr, inembed.ctypes.data, 1, fftobj.idist,
                             tmpout.ptr, onembed.ctypes.data, 1, fftobj.odist,
                             FFTW_BACKWARD, flags)
        # R2C or C2R (hence no direction argument for plan creation)
        else:
            plan = plan_func(1, n.ctypes.data, fftobj.nbatch,
                             tmpin.ptr, inembed.ctypes.data, 1, fftobj.idist,
                             tmpout.ptr, onembed.ctypes.data, 1, fftobj.odist,
                             flags)
        del tmpin
        del tmpout
        return plan

class FFT(_BaseFFT):
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(FFT, self).__init__(invec, outvec, nbatch, size)
        self.iptr = self.invec.ptr
        self.optr = self.outvec.ptr
        self._efunc = execute_function[str(self.invec.dtype)][str(self.outvec.dtype)]
        self._efunc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
        self.plan = _fftw_setup(self)

    def execute(self):
        self._efunc(self.plan, self.iptr, self.optr)

class IFFT(_BaseIFFT):
    def __init__(self, invec, outvec, nbatch=1, size=None):
        super(IFFT, self).__init__(invec, outvec, nbatch, size)
        self.iptr = self.invec.ptr
        self.optr = self.outvec.ptr
        self._efunc = execute_function[str(self.invec.dtype)][str(self.outvec.dtype)]
        self._efunc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
        self.plan = _fftw_setup(self)

    def execute(self):
        self._efunc(self.plan, self.iptr, self.optr)

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
                      help="Give 'openmp', 'pthreads' or 'unthreaded' to specify which threaded FFTW to use",
                      default=None)
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
    optgroup.add_argument("--fftw-import-system-wisdom",
                          help = "If given, call fftw[f]_import_system_wisdom()",
                          action = "store_true")

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

    if opt.fftw_import_system_wisdom and ((opt.fftw_input_float_wisdom_file is not None) 
                                          or (opt.fftw_input_double_wisdom_file is not None)):
        parser.error("If --fftw-import-system-wisdom is given, then you cannot give"
                     " either of --fftw-input-float-wisdom-file or --fftw-input-double-wisdom-file")

    if opt.fftw_threads_backend is not None:
        if opt.fftw_threads_backend not in ['openmp','pthreads','unthreaded']:
            parser.error("Invalid threads backend; must be 'openmp', 'pthreads' or 'unthreaded'")

def from_cli(opt):
    # Since opt.fftw_threads_backend defaults to None, the following is always
    # appropriate:
    set_threads_backend(opt.fftw_threads_backend)

    # Import system wisdom.  Should really add error checking and logging to that...
    if opt.fftw_import_system_wisdom:
        import_sys_wisdom()

    # Read specified user-provided wisdom files
    if opt.fftw_input_float_wisdom_file is not None:
        import_single_wisdom_from_filename(opt.fftw_input_float_wisdom_file)        

    if opt.fftw_input_double_wisdom_file is not None:
        import_double_wisdom_from_filename(opt.fftw_input_double_wisdom_file)        

    # Set the user-provided measure level
    set_measure_level(opt.fftw_measure_level)
