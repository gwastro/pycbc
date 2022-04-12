# cython: profile=True
import numpy
cimport numpy
from libc.stdlib cimport malloc, free
from libc.math cimport cos, sin, sqrt # This imports c's sin and cos function from the math library
from cython import wraparound, boundscheck, cdivision
from pycbc.types import real_same_precision_as, complex_same_precision_as

ctypedef fused REALTYPE:
    float
    double

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

@boundscheck(False)
@wraparound(False)
def chisq_accum_bin_cython(numpy.ndarray[REALTYPE, ndim=1] chisq,
                           numpy.ndarray[COMPLEXTYPE, ndim=1] q, int N):
    # This can be made parallelized?!
    for i in range(N):
        chisq[i] += q[i].real*q[i].real+q[i].imag*q[i].imag

@boundscheck(False)
@wraparound(False)
@cdivision(True)
def point_chisq_code(numpy.ndarray[REALTYPE, ndim=1] chisq,
                     numpy.ndarray[COMPLEXTYPE, ndim=1] v1,
                     int n, int slen,
                     numpy.ndarray[REALTYPE, ndim=1] shifts,
                     numpy.ndarray[numpy.uint32_t, ndim=1] bins,
                     int blen):
    # Do I need to declare UINT vs INT??
    cdef int num_parallel_regions, bstart, bend, i, j, k, r, start, end
    cdef REALTYPE *outr
    cdef REALTYPE *outi
    cdef REALTYPE *pr
    cdef REALTYPE *pi
    cdef REALTYPE *vsr
    cdef REALTYPE *vsi
    cdef REALTYPE *outr_tmp
    cdef REALTYPE *outi_tmp
    cdef COMPLEXTYPE v
    cdef REALTYPE vr, vi, t1, t2, k1, k2, k3, vs, va

    num_parallel_regions = 1

    outr = <REALTYPE *> malloc(n * sizeof(REALTYPE))
    outi = <REALTYPE *> malloc(n * sizeof(REALTYPE))
    pr = <REALTYPE *> malloc(n * sizeof(REALTYPE))
    pi = <REALTYPE *> malloc(n * sizeof(REALTYPE))
    vsr = <REALTYPE *> malloc(n * sizeof(REALTYPE))
    vsi = <REALTYPE *> malloc(n * sizeof(REALTYPE))
    outr_tmp = <REALTYPE *> malloc(n * sizeof(REALTYPE))
    outi_tmp = <REALTYPE *> malloc(n * sizeof(REALTYPE))

    for r in range(blen):
        bstart = bins[r] # int
        bend = bins[r+1] # int
        blen = bend - bstart # int

        #outr = numpy.zeros(n, dtype=real_type)
        #outi = numpy.zeros(n, dtype=real_type)
        for i in range(n):
            outr[i] = 0.
            outi[i] = 0.


        # CAN WE MAKE THIS LOOP PARALLEL?!
        for k in range(num_parallel_regions):
            start = blen * k / num_parallel_regions + bstart # uint
            end = blen * (k + 1) / num_parallel_regions + bstart # uint

            # start the cumulative rotations at the offset point
            #pr = numpy.zeros(n, dtype=real_type)
            #pi = numpy.zeros(n, dtype=real_type)
            #vsr = numpy.zeros(n, dtype=real_type)
            #vsi = numpy.zeros(n, dtype=real_type)
            #outr_tmp = numpy.zeros(n, dtype=real_type)
            #outi_tmp = numpy.zeros(n, dtype=real_type)

            for i in range(n):
                pr[i] = cos(2 * 3.141592653 * shifts[i] * (start) / slen)
                pi[i] = sin(2 * 3.141592653 * shifts[i] * (start) / slen)
                vsr[i] = cos(2 * 3.141592653 * shifts[i] / slen)
                vsi[i] = sin(2 * 3.141592653 * shifts[i] / slen)
                outr_tmp[i] = 0
                outi_tmp[i] = 0

            for j in range(start, end): # uint
                v = v1[j]
                vr = v.real
                vi = v.imag
                vs = vr + vi;
                va = vi - vr;

                for i in range(n):
                    t1 = pr[i]
                    t2 = pi[i]

                    # Complex multiply pr[i] * v
                    k1 = vr * (t1 + t2)
                    k2 = t1 * va
                    k3 = t2 * vs

                    outr_tmp[i] += k1 - k3
                    outi_tmp[i] += k1 + k2

                    # phase shift for the next time point
                    pr[i] = t1 * vsr[i] - t2 * vsi[i]
                    pi[i] = t1 * vsi[i] + t2 * vsr[i]

            # WHAT DOES THIS MEAN?? : pragma omp critical
            for i in range(n):
                outr[i] += outr_tmp[i]
                outi[i] += outi_tmp[i]


        for i in range(n):
            chisq[i] += outr[i]*outr[i] + outi[i]*outi[i]

    free(outr)
    free(outi)
    free(pr)
    free(pi)
    free(vsr)
    free(vsi)
    free(outr_tmp)
    free(outi_tmp)

@boundscheck(False)
@wraparound(False)
@cdivision(True)
def point_chisq_multibin_code(numpy.ndarray[COMPLEXTYPE, ndim=2] chisq,
                              numpy.ndarray[COMPLEXTYPE, ndim=1] v1,
                              int n, int slen,
                              numpy.ndarray[REALTYPE, ndim=1] shifts,
                              numpy.ndarray[COMPLEXTYPE, ndim=1] scales,
                              numpy.ndarray[numpy.uint32_t, ndim=2] bins,
                              int bnum):

    cdef int blen, i, j, k, bi
    cdef REALTYPE pr, pi, sr, si, scale_r, scale_i, phase_diff_r, phase_diff_i, vr, vi, bsum_r, bsum_i
    cdef COMPLEXTYPE s, v

    for i in range(n):
        k = 0
        k += bins[i, 0]

        pr = cos(2 * 3.141592653 * shifts[i] * (k) / slen)
        pi = sin(2 * 3.141592653 * shifts[i] * (k) / slen)
        
        s = scales[i]
        sr = s.real
        si = s.imag
        
        scale_r = pr * sr - pi * si
        scale_i = pr * si + pi * sr
        
        phase_diff_r = cos(2 * 3.141592653 * shifts[i] / slen)
        phase_diff_i = sin(2 * 3.141592653 * shifts[i] / slen)

        for j in range(bnum):
            blen = bins[i, j+1] - bins[i, j]
            bsum_r = 0
            bsum_i = 0
            
            for bi in range(blen):
                v = v1[k]
                vr = v.real
                vi = v.imag
                
                sr = scale_r
                si = scale_i

                bsum_r += vr * sr - vi * si
                bsum_i += vr * si + vi * sr

                k += 1
                scale_r = sr * phase_diff_r - si * phase_diff_i
                scale_i = sr * phase_diff_i + si * phase_diff_r

            chisq[i, j].real = bsum_r
            chisq[i, j].imag = bsum_i

def chisq_accum_bin_numpy(chisq, q):
    chisq += q.squared_norm()

def chisq_accum_bin(chisq, q):
    chisq = numpy.array(chisq.data, copy=False)
    q = numpy.array(q.data, copy=False)
    N = len(chisq)
    chisq_accum_bin_cython(chisq, q, N)

def shift_sum(v1, shifts, bins):
    real_type = real_same_precision_as(v1)
    shifts = numpy.array(shifts, dtype=real_type)

    bins = numpy.array(bins, dtype=numpy.uint32)
    blen = len(bins) - 1 #
    v1 = numpy.array(v1.data, copy=False)
    slen = len(v1)
    n = int(len(shifts))

    # Create some output memory
    chisq = numpy.zeros(n, dtype=real_type)

    point_chisq_code(chisq, v1, n, slen, shifts, bins, blen)

    return  chisq

def shift_multibin(v1, shifts, scales, bins):
    real_type = real_same_precision_as(v1)
    complex_type = complex_same_precision_as(v1)
    
    shifts = numpy.array(shifts, dtype=real_type)
    scales = numpy.array(scales, dtype=complex_type)

    bins = numpy.array(bins, dtype=numpy.uint32)
    bnum = bins.shape[1] - 1
    
    v1 = numpy.array(v1.data, copy=False)
    slen = len(v1)
    n = int(len(shifts))

    # Create some output memory
    chisq = numpy.zeros((n, bnum), dtype=complex_type)

    point_chisq_multibin_code(chisq, v1, n, slen, shifts, scales, bins, bnum)

    return chisq
