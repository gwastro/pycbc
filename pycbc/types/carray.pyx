cimport cython, numpy

ctypedef fused A:
    float
    double

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_real(numpy.ndarray [A, ndim = 1] a, numpy.ndarray [A, ndim = 1] b):
    cdef double total = 0
    cdef unsigned int xmax = a.shape[0]
    cdef unsigned int i
    
    cdef A* x = &a[0]
    cdef A* y = &a[0]

    for i in range(xmax):
        total += x[i] * y[i]

    return total
    
ctypedef fused B:
    float complex
    double complex

def abs_arg_max_complex(numpy.ndarray [B, ndim=1] a):
    cdef unsigned int xmax = a.shape[0]
    cdef double mag
    cdef double magmax = 0
    cdef unsigned int idx = 0

    for i in range(xmax):
        mag = a[i].real * a[i].real + a[i].imag * a[i].imag
        if mag > magmax:
            magmax = mag
            idx = i
            
    return idx
       

