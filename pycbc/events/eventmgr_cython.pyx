import numpy as np
cimport numpy as cnp
from cython import wraparound, boundscheck, cdivision


ctypedef fused REALTYPE:
    float
    double


@boundscheck(False)
@wraparound(False)
@cdivision(True)
def findchirp_cluster_over_window_cython\
        (cnp.ndarray[cnp.int32_t, ndim=1] times,
         cnp.ndarray[REALTYPE, ndim=1] absvalues, int window_length,
         cnp.ndarray[cnp.int32_t, ndim=1] indices, int tlen):
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


@boundscheck(False)
@wraparound(False)
@cdivision(True)
def test_internal_one_three(
    double[:] tdif,
    double[:] pdif,
    double[:] sdif,
    int[:] tbin,
    int[:] pbin,
    int[:] sbin,
    double twidth,
    double pwidth,
    double swidth,
    int length
):
    cdef:
        int idx

    for idx in range(length):
        tbin[idx] = <int>(tdif[idx] / twidth)
        pbin[idx] = <int>(pdif[idx] / pwidth)
        sbin[idx] = <int>(sdif[idx] / swidth)

    return None

