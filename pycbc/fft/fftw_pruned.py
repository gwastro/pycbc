"""This module provides a functions to perform a pruned FFT based on FFTW

This should be considered a test and example module, as the functionality
can and should be generalized to other FFT backends, and precisions.

These functions largely implemented the generic FFT decomposition as
described rather nicely by wikipedia.

http://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

I use a similar naming convention here, with minor simplifications to the
twiddle factors.
"""
import numpy, ctypes, pycbc.types
from pycbc.libutils import get_ctypes_library
import logging
from .fftw_pruned_cython import second_phase_cython

warn_msg = ("The FFTW_pruned module can be used to speed up computing SNR "
            "timeseries by computing first at a low sample rate and then "
            "computing at full sample rate only at certain samples. This code "
            "has not yet been used in production, and has no test case. "
            "This was also ported to Cython in this state. "
            "This code would need verification before trusting results. "
            "Please do contribute test cases.")

logging.warning(warn_msg)

# FFTW constants
FFTW_FORWARD = -1
FFTW_BACKWARD = 1
FFTW_MEASURE = 0
FFTW_PATIENT = 1 << 5
FFTW_ESTIMATE = 1 << 6
float_lib = get_ctypes_library('fftw3f', ['fftw3f'],mode=ctypes.RTLD_GLOBAL)
fexecute = float_lib.fftwf_execute_dft
fexecute.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]

ftexecute = float_lib.fftwf_execute_dft
ftexecute.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]

def plan_transpose(N1, N2):
    """
    Create a plan for transposing internally to the pruned_FFT calculation.
    (Alex to provide a write up with more details.)

    Parameters
    -----------
    N1 : int
        Number of rows.
    N2 : int
        Number of columns.

    Returns
    --------
    plan : FFTWF plan
        The plan for performing the FFTW transpose.
    """

    rows = N1
    cols = N2

    iodim = numpy.zeros(6, dtype=numpy.int32)
    iodim[0] = rows
    iodim[1] = 1
    iodim[2] = cols
    iodim[3] = cols
    iodim[4] = rows
    iodim[5] = 1

    N = N1*N2
    vin = pycbc.types.zeros(N, dtype=numpy.complex64)
    vout = pycbc.types.zeros(N, dtype=numpy.complex64)

    f = float_lib.fftwf_plan_guru_dft
    f.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                  ctypes.c_void_p, ctypes.c_void_p,
                  ctypes.c_void_p, ctypes.c_void_p,
                  ctypes.c_int]
    f.restype = ctypes.c_void_p
    return f(0, None, 2, iodim.ctypes.data, vin.ptr, vout.ptr, None, FFTW_MEASURE)

def plan_first_phase(N1, N2):
    """
    Create a plan for the first stage of the pruned FFT operation.
    (Alex to provide a write up with more details.)

    Parameters
    -----------
    N1 : int
        Number of rows.
    N2 : int
        Number of columns.

    Returns
    --------
    plan : FFTWF plan
        The plan for performing the first phase FFT.
    """
    N = N1*N2
    vin = pycbc.types.zeros(N, dtype=numpy.complex64)
    vout = pycbc.types.zeros(N, dtype=numpy.complex64)
    f = float_lib.fftwf_plan_many_dft
    f.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                  ctypes.c_void_p, ctypes.c_void_p,
                  ctypes.c_int, ctypes.c_int,
                  ctypes.c_void_p, ctypes.c_void_p,
                  ctypes.c_int, ctypes.c_int,
                  ctypes.c_int, ctypes.c_int]
    f.restype = ctypes.c_void_p
    return f(1, ctypes.byref(ctypes.c_int(N2)), N1,
             vin.ptr, None, 1, N2,
             vout.ptr, None, 1, N2, FFTW_BACKWARD, FFTW_MEASURE)

_theplan = None
def first_phase(invec, outvec, N1, N2):
    """
    This implements the first phase of the FFT decomposition, using
    the standard FFT many plans.

    Parameters
    -----------
    invec : array
        The input array.
    outvec : array
        The output array.
    N1 : int
        Number of rows.
    N2 : int
        Number of columns.
    """
    global _theplan
    if _theplan is None:
        _theplan = plan_first_phase(N1, N2)
    fexecute(_theplan, invec.ptr, outvec.ptr)

def second_phase(invec, indices, N1, N2):
    """
    This is the second phase of the FFT decomposition that actually performs
    the pruning. It is an explicit calculation for the subset of points. Note
    that there seem to be some numerical accumulation issues at various values
    of N1 and N2.

    Parameters
    ----------
    invec :
        The result of the first phase FFT
    indices : array of ints
        The index locations to calculate the FFT
    N1 : int
        The length of the second phase "FFT"
    N2 : int
        The length of the first phase FFT

    Returns
    -------
    out : array of floats
    """
    invec = numpy.array(invec.data, copy=False)
    NI = len(indices) # pylint:disable=unused-variable
    N1 = int(N1)
    N2 = int(N2)
    out = numpy.zeros(len(indices), dtype=numpy.complex64)
    indices = numpy.array(indices, dtype=numpy.uint32)

    # Note, the next step if this needs to be faster is to invert the loops
    second_phase_cython(N1, N2, NI, indices, out, invec)

    return out

_thetransposeplan = None
def fft_transpose_fftw(vec):
    """
    Perform an FFT transpose from vec into outvec.
    (Alex to provide more details in a write-up.)

    Parameters
    -----------
    vec : array
        Input array.

    Returns
    --------
    outvec : array
        Transposed output array.
    """
    global _thetransposeplan
    outvec = pycbc.types.zeros(len(vec), dtype=vec.dtype)
    if _theplan is None:
        N1, N2 = splay(vec)
        _thetransposeplan = plan_transpose(N1, N2)
    ftexecute(_thetransposeplan, vec.ptr, outvec.ptr)
    return  outvec

fft_transpose = fft_transpose_fftw

def splay(vec):
    """ Determine two lengths to split stride the input vector by
    """
    N2 = 2 ** int(numpy.log2( len(vec) ) / 2)
    N1 = len(vec) / N2
    return N1, N2

def pruned_c2cifft(invec, outvec, indices, pretransposed=False):
    """
    Perform a pruned iFFT, only valid for power of 2 iffts as the
    decomposition is easier to choose. This is not a strict requirement of the
    functions, but it is unlikely to the optimal to use anything but power
    of 2. (Alex to provide more details in write up.

    Parameters
    -----------
    invec : array
        The input vector. This should be the correlation between the data and
        the template at full sample rate. Ideally this is pre-transposed, but
        if not this will be transposed in this function.
    outvec : array
        The output of the first phase of the pruned FFT.
    indices : array of ints
        The indexes at which to calculate the full sample-rate SNR.
    pretransposed : boolean, default=False
        Used to indicate whether or not invec is pretransposed.

    Returns
    --------
    SNRs : array
        The complex SNRs at the indexes given by indices.
    """
    N1, N2 = splay(invec)

    if not pretransposed:
        invec = fft_transpose(invec)
    first_phase(invec, outvec, N1=N1, N2=N2)
    out = second_phase(outvec, indices, N1=N1, N2=N2)
    return out
