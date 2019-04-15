%%cython -a

import numpy
cimport numpy
#from libc.stdlib cimport malloc, free
#from libc.complex cimport creal, cimag
from libc.math cimport sin, cos # This imports c's sin and cos function from the math library
from cython import wraparound, boundscheck, cdivision

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

@boundscheck(False)
@wraparound(False)
@cdivision(True)
def second_phase_cython(int N1, int N2, int NI, numpy.ndarray[numpy.uint32_t, ndim=1] indices,
                        numpy.ndarray[numpy.complex64_t, ndim=1] out,
                        numpy.ndarray[COMPLEXTYPE, ndim=1] invec):
    cdef float pi = 3.14159265359
    cdef int N = N1 * N2
    cdef numpy.uint32_t k, n1, k2
    cdef float phase_inc, sp, sc, n1f, curr_angle
    cdef COMPLEXTYPE val
    for i in range(NI):
        val = 0 + 0j
        k = indices[i]
        k2 = k % N2
        phase_inc = (2 * pi * k) / <float> N
        for n1 in range(N1):
            n1f = <float> n1
            # NOTE: Previously this used sincosf, which doesn't seem to be present in Cython
            curr_angle = phase_inc * n1
            sp = sin(curr_angle)
            cp = cos(curr_angle)
            val += (cp + sp * 1j) * invec[k2 + N2*n1]  
        out[i] = val


@boundscheck(False)
@wraparound(False)
@cdivision(True)
def fast_second_phase_cython(int N1, int N2, int NI,
                             numpy.ndarray[numpy.uint32_t, ndim=1] indices,
                             numpy.ndarray[numpy.complex64_t, ndim=1] out,
                             numpy.ndarray[COMPLEXTYPE, ndim=1] invec):
    cdef float pi = 3.14159265359
    cdef int N = N1 * N2
    cdef float sp, cp, phase_inc, n1f
    cdef numpy.uint32_t k, k2, n1
    cdef COMPLEXTYPE val, twiddle_inc, twiddle
    
    for i in range(NI):
        val = 0 + 0j
        k = indices[i]
        k2 = k % N2
        phase_inc = (2 * pi * k) / <float> N
        sp = sin(phase_inc)
        cp = cos(phase_inc)
        twiddle_inc = cp + cp * 1j
        twiddle = 1 + 0j
        
        for n1 in range(N1):
            n1f = <float> n1
            val += twiddle * invec[k2 + N2*n1]
            twiddle *= twiddle_inc
        out[i] = val
