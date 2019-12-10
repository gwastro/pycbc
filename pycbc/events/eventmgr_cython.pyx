import numpy
cimport numpy
from cython import wraparound, boundscheck, cdivision


ctypedef fused REALTYPE:
    float
    double


@boundscheck(False)
@wraparound(False)
@cdivision(True)
def findchirp_cluster_over_window_cython\
        (numpy.ndarray[numpy.int32_t, ndim=1] times,
         numpy.ndarray[REALTYPE, ndim=1] absvalues, int window_length,
         numpy.ndarray[numpy.int32_t, ndim=1] indices, int tlen):
    cdef int j = 0
    cdef int curr_ind = 0
    cdef int i
    
    for i in range(tlen):
        if ((times[i] - times[curr_ind]) > window_length):
            j += 1
            indices[j] = i
            curr_ind = i
        elif (absvalues[i] > absvalues[curr_ind]):
            indices[j] = i
            curr_ind = i
    return j
