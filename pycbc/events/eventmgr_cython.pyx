import numpy as np
cimport numpy as cnp
from cython import wraparound, boundscheck, cdivision
from libc.math cimport M_PI, sqrt


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
def logsignalrateinternals_computepsignalbins(
    double[:] pdif,
    double[:] tdif,
    double[:] sdif,
    int[:] pbin,
    int[:] tbin,
    int[:] sbin,
    float[:] p,
    double[:] t,
    float[:] s,
    float[:] sig,
    float[:] pref,
    double[:] tref,
    float[:] sref,
    float[:] sigref,
    double[:] shift,
    long int[:] rtype,
    double sense,
    double senseref,
    double twidth,
    double pwidth,
    double swidth,
    int to_shift_ref,
    int to_shift_ifo,
    int length
):
    cdef:
        int idx, ridx

    for idx in range(length):
        ridx = rtype[idx]
        pdif[idx] = (pref[ridx] - p[ridx]) % (M_PI * 2)
        if pdif[idx] < 0:
            # C modulus operator is not same as python's, correct for this
            pdif[idx] += (M_PI * 2)
        tdif[idx] = shift[ridx] * to_shift_ref + tref[ridx] - shift[ridx] * to_shift_ifo - t[ridx]
        sdif[idx] = (s[ridx] * sense * sqrt(sigref[ridx])) / (sref[ridx] * senseref * sqrt(sig[ridx]))

    for idx in range(length):
        tbin[idx] = <int>(tdif[idx] / twidth)
        pbin[idx] = <int>(pdif[idx] / pwidth)
        sbin[idx] = <int>(sdif[idx] / swidth)

@boundscheck(False)
@wraparound(False)
@cdivision(True)
def logsignalrateinternals_compute2detrate(
    int[:] nbinned0,
    int[:] nbinned1,
    int[:] nbinned2,
    long int c0_size,
    long int c1_size,
    long int c2_size,
    float[:] rate,
    long int[:] rtype,
    float[:] sref,
    float[:,:,::1] two_det_weights, # This declares a C-contiguous array
    float max_penalty,
    float ref_snr,
    int length
):
    cdef:
        int idx, ridx, id0, id1, id2
        float rescale_fac

    for idx in range(length):
        ridx = rtype[idx]
        id0 = nbinned0[idx] + c0_size / 2
        id1 = nbinned1[idx] + c1_size / 2
        id2 = nbinned2[idx] + c2_size / 2

        # For bins that exist in the signal pdf histogram, apply that pdf value,
        # otherwise apply the "max penalty" value.
        # The bins are specified by 3 indexes (corresponding to the time
        # difference, phase difference and relative sensitivity dimensions).
        if (id0 > 0) and (id0 < c0_size) and (id1 > 0) and (id1 < c1_size) and (id2 > 0) and (id2 < c2_size):
            rate[ridx] = two_det_weights[id0, id1, id2]
        else:
            rate[ridx] = max_penalty
        # Scale by signal population SNR
        rescale_fac = ref_snr / sref[ridx]
        rate[ridx] *= (rescale_fac*rescale_fac*rescale_fac*rescale_fac)


@boundscheck(False)
@wraparound(False)
@cdivision(True)
def coincbuffer_expireelements(
    float[:] cbuffer,
    int[:] timer1,
    int[:] timer2,
    int time1,
    int time2,
    int expiration,
    int length,
):
    cdef:
        int idx, keep_count

    keep_count = 0
    for idx in range(length):
        if (timer1[idx] >= (time1 - expiration)) and (timer2[idx] >= (time2 - expiration)):
            cbuffer[keep_count] = cbuffer[idx]
            timer1[keep_count] = timer1[idx]
            timer2[keep_count] = timer2[idx]
            keep_count += 1

    return keep_count


@boundscheck(False)
@wraparound(False)
@cdivision(True)
def coincbuffer_numgreater(
    float[:] cbuffer,
    int length,
    float value
):
    cdef:
        int idx, count

    count = 0
    for idx in range(length):
        if cbuffer[idx] > value:
            count += 1
    return count
