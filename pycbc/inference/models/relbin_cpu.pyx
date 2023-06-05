""" Optimized inner loop functions for the relative likelihood model
"""
cdef extern from "complex.h":
    double complex exp(double complex)
    double norm(double complex)
    double complex conj(double complex)
    double real(double complex)
    double imag(double complex)

import numpy
cimport cython, numpy

# Used for calculating the cross terms of
# two signals when analyzing multiple signals
cpdef likelihood_parts_multi(double [::1] freqs,
                     double fp,
                     double fc,
                     double dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double fp2,
                     double fc2,
                     double dtc2,
                     double complex[::1] hp2,
                     double complex[::1] hc2,
                     double complex[::1] h002,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1

    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(-2.0j * 3.141592653 * dtc * freqs[i])
               * (fp * hp[i] + fc * hc[i])) / h00[i]
        r0n *= conj((exp(-2.0j * 3.141592653 * dtc2 * freqs[i])
               * (fp2 * hp2[i] + fc2 * hc2[i])) / h002[i])
        r1 = r0n - r0
        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1

        r0 = r0n
    return hd

# Used for calculating the cross terms of
# two signals when analyzing multiple signals
# + allows for frequency-varying antenna response
cpdef likelihood_parts_multi_v(double [::1] freqs,
                     double[::1] fp,
                     double[::1] fc,
                     double[::1] dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double[::1] fp2,
                     double[::1] fc2,
                     double[::1] dtc2,
                     double complex[::1] hp2,
                     double complex[::1] hc2,
                     double complex[::1] h002,
                     double complex[::1] a0,
                     double complex[::1] a1,

                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1

    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(-2.0j * 3.141592653 * dtc[i] * freqs[i])
               * (fp[i] * hp[i] + fc[i] * hc[i])) / h00[i]
        r0n *= conj((exp(-2.0j * 3.141592653 * dtc2[i] * freqs[i])
               * (fp2[i] * hp2[i] + fc2[i] * hc2[i])) / h002[i])
        r1 = r0n - r0

        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1

        r0 = r0n
    return hd

# Standard likelihood
cpdef likelihood_parts(double [::1] freqs,
                     double fp,
                     double fc,
                     double dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1, x0, x1, x0n;
    cdef double hh=0

    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(-2.0j * 3.141592653 * dtc * freqs[i])
               * (fp * hp[i] + fc * hc[i])) / h00[i]
        r1 = r0n - r0

        x0n = norm(r0n)
        x1 = x0n - x0

        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1
            hh += real(b0[i-1] * x0 + b1[i-1] * x1)

        r0 = r0n
        x0 = x0n
    return conj(hd), hh

# Likelihood where no antenna response is applied
cpdef likelihood_parts_det(double [::1] freqs,
                     double dtc,
                     double complex[::1] hp,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1, x0, x1, x0n;
    cdef double hh=0
    cdef int N

    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(-2.0j * 3.141592653 * dtc * freqs[i])
               * (hp[i])) / h00[i]
        r1 = r0n - r0

        x0n = norm(r0n)
        x1 = x0n - x0

        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1
            hh += real(b0[i-1] * x0 + b1[i-1] * x1)

        r0 = r0n
        x0 = x0n
    return conj(hd), hh

# Used where the antenna response may be frequency varying
cpdef likelihood_parts_v(double [::1] freqs,
                     double[::1] fp,
                     double[::1] fc,
                     double[::1] dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1, x0, x0n, x1
    cdef double hh=0

    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(-2.0j * 3.141592653 * dtc[i] * freqs[i])
               * (fp[i] * hp[i] + fc[i] * hc[i])) / h00[i]
        r1 = r0n - r0

        x0n = norm(r0n)
        x1 = x0n - x0

        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1
            hh += real(b0[i-1] * x0 + b1[i-1] * x1)

        r0 = r0n
        x0 = x0n
    return conj(hd), hh

# Used where the antenna response may be frequency varying
# and there is a polarization vector marginalization
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef likelihood_parts_v_pol(double [::1] freqs,
                     double[::1] fp,
                     double[::1] fc,
                     double[::1] dtc,
                     double complex[::1] pol_phase,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1, x0, x0n, x1, fp2, fc2
    cdef double hh=0

    N = freqs.shape[0]
    num_samples = pol_phase.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)

    for j in range(num_samples):
        hh = 0
        hd = 0
        for i in range(N):
        
            f = (fp[i] + 1.0j * fc[i]) * pol_phase[j]
            fp2 = real(f)
            fc2 = imag(f)
        
            r0n = (exp(-2.0j * 3.141592653 * dtc[i] * freqs[i])
                   * (fp2 * hp[i] + fc2 * hc[i])) / h00[i]
            r1 = r0n - r0

            x0n = norm(r0n)
            x1 = x0n - x0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1
                hh += real(b0[i-1] * x0 + b1[i-1] * x1)

            r0 = r0n
            x0 = x0n
            
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv

# Used where the antenna response may be frequency varying
# and there is a polarization vector marginalization
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef likelihood_parts_v_time(double [::1] freqs,
                     double[::1] fp,
                     double[::1] fc,
                     double[::1] times,
                     double[::1] dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1, x0, x0n, x1
    cdef double hh=0, ttime;

    N = freqs.shape[0]
    num_samples = dtc.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)

    for j in range(num_samples):
        hh = 0
        hd = 0
        for i in range(N):   
            # This allows for multiple time offsets
            ttime = times[i] + dtc[j]
            r0n = (exp(-2.0j * 3.141592653 * ttime * freqs[i])
                   * (fp[i] * hp[i] + fc[i] * hc[i])) / h00[i]
            r1 = r0n - r0

            x0n = norm(r0n)
            x1 = x0n - x0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1
                hh += real(b0[i-1] * x0 + b1[i-1] * x1)

            r0 = r0n
            x0 = x0n
            
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv

# Used where the antenna response may be frequency varying
# and there is a polarization vector marginalization
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef likelihood_parts_v_pol_time(double [::1] freqs,
                     double[::1] fp,
                     double[::1] fc,
                     double[::1] times,
                     double[::1] dtc,
                     double complex[::1] pol_phase,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1, x0, x0n, x1, fp2, fc2
    cdef double hh=0, ttime;

    N = freqs.shape[0]
    num_samples = pol_phase.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)

    for j in range(num_samples):
        hh = 0
        hd = 0
        for i in range(N):
        
            f = (fp[i] + 1.0j * fc[i]) * pol_phase[j]
            fp2 = real(f)
            fc2 = imag(f)
        
            # This allows for multiple time offsets
            ttime = times[i] + dtc[j]
            r0n = (exp(-2.0j * 3.141592653 * ttime * freqs[i])
                   * (fp2 * hp[i] + fc2 * hc[i])) / h00[i]
            r1 = r0n - r0

            x0n = norm(r0n)
            x1 = x0n - x0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1
                hh += real(b0[i-1] * x0 + b1[i-1] * x1)

            r0 = r0n
            x0 = x0n
            
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv

# Standard likelihood but simultaneously handling multiple sky or time points
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef likelihood_parts_vector(double [::1] freqs,
                     double[::1] fp,
                     double[::1] fc,
                     double[::1] dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd, r0, r0n, r1, x0, x0n, x1
    cdef double hh
    N = freqs.shape[0]
    num_samples = fp.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)

    for j in range(num_samples):
        hd = 0
        hh = 0
        for i in range(N):
            r0n = (exp(-2.0j * 3.141592653 * dtc[j] * freqs[i])
                   * (fp[j] * hp[i] + fc[j] * hc[i])) / h00[i]
            r1 = r0n - r0

            x0n = norm(r0n)
            x1 = x0n - x0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1
                hh += real(b0[i-1] * x0 + b1[i-1] * x1)

            r0 = r0n
            x0 = x0n
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv
    
# Standard likelihood but simultaneously handling multiple time points
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef likelihood_parts_vectort(double [::1] freqs,
                     double fp,
                     double fc,
                     double[::1] dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd, r0, r0n, r1, x0, x0n, x1
    cdef double hh
    N = freqs.shape[0]
    num_samples = dtc.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)

    for j in range(num_samples):
        hd = 0
        hh = 0
        for i in range(N):
            r0n = (exp(-2.0j * 3.141592653 * dtc[j] * freqs[i])
                   * (fp * hp[i] + fc * hc[i])) / h00[i]
            r1 = r0n - r0

            x0n = norm(r0n)
            x1 = x0n - x0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1
                hh += real(b0[i-1] * x0 + b1[i-1] * x1)

            r0 = r0n
            x0 = x0n
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv

# Like above, but if only polarization is marginalized
# this is a slow implementation and the loop should be inverted /
# refactored to do each pol separately and then combine
# included as is for testing purposes
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef likelihood_parts_vectorp(double [::1] freqs,
                     double[::1] fp,
                     double[::1] fc,
                     double dtc,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ) :
    cdef size_t i
    cdef double complex hd, r0, r0n, r1, x0, x0n, x1
    cdef double hh
    N = freqs.shape[0]
    num_samples = fp.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)

    for j in range(num_samples):
        hd = 0
        hh = 0
        for i in range(N):
            r0n = (exp(-2.0j * 3.141592653 * dtc * freqs[i])
                   * (fp[j] * hp[i] + fc[j] * hc[i])) / h00[i]
            r1 = r0n - r0

            x0n = norm(r0n)
            x1 = x0n - x0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1
                hh += real(b0[i-1] * x0 + b1[i-1] * x1)

            r0 = r0n
            x0 = x0n
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv

# Standard likelihood but simultaneously handling multiple sky or time points
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef snr_predictor(double [::1] freqs,
                     double tstart,
                     double delta_t,
                     int num_samples,
                     double complex[::1] hp,
                     double complex[::1] hc,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ):
    cdef size_t i
    cdef double complex hd, r0, r0n, r1, x0, x0n, x1
    cdef double hh

    cdef double complex chd, cr0, cr0n, cr1, cx0, cx0n, cx1
    cdef double chh

    N = freqs.shape[0]

    cdef numpy.ndarray[numpy.float64_t, ndim=1] snr = numpy.empty(num_samples, dtype=numpy.float64)
    cdef numpy.ndarray[numpy.complex128_t, ndim=1] twiddle = numpy.empty(N, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.complex128_t, ndim=1] rotate = numpy.empty(N, dtype=numpy.complex128)

    hh = 0
    chh = 0
    for i in range(N):
        twiddle[i] = exp(-2.0j * 3.141592653 * tstart * freqs[i]) / h00[i]
        rotate[i] =  exp(-2.0j * 3.141592653 * delta_t * freqs[i])

        r0n =  hp[i] / h00[i]
        cr0n = hc[i] / h00[i]

        x0n = norm(r0n)
        cx0n = norm(cr0n)
        x1 = x0n - x0
        cx1 = cx0n - cx0
        if i > 0:
            hh += real(b0[i-1] * x0 + b1[i-1] * x1)
            chh += real(b0[i-1] * cx0 + b1[i-1] * cx1)

        x0 = x0n
        cx0 = cx0n

    for j in range(num_samples):
        hd = 0
        chd = 0

        # Calculate the average SNR for the hp / hc waveforms
        for i in range(N):
            r0n =  twiddle[i] * hp[i]
            cr0n = twiddle[i] * hc[i]

            r1 = r0n - r0
            cr1 = cr0n - cr0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1
                chd += a0[i-1] * cr0 + a1[i-1] * cr1

            r0 = r0n
            cr0 = cr0n

            # Time shift the twiddle factors to the next time samples
            twiddle[i] *= rotate[i]

        snr[j] = (norm(hd) / hh / 2.0 + norm(chd) / chh / 2.0) ** 0.5
    return snr


# calculate a rough SNR for use in predicting sky location consistency
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cpdef snr_predictor_dom(double [::1] freqs,
                     double tstart,
                     double delta_t,
                     int num_samples,
                     double complex[::1] hp,
                     double complex[::1] h00,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     double [::1] b0,
                     double [::1] b1,
                     ):
    cdef size_t i
    cdef double complex hd, r0, r0n, r1, x0, x0n, x1
    cdef double hh

    N = freqs.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] sh = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.complex128_t, ndim=1] twiddle = numpy.empty(N, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.complex128_t, ndim=1] rotate = numpy.empty(N, dtype=numpy.complex128)

    hh = 0
    for i in range(N):
        twiddle[i] = exp(-2.0j * 3.141592653 * tstart * freqs[i]) / h00[i]
        rotate[i] =  exp(-2.0j * 3.141592653 * delta_t * freqs[i])

        r0n =  hp[i] / h00[i]

        x0n = norm(r0n)
        x1 = x0n - x0
        if i > 0:
            hh += real(b0[i-1] * x0 + b1[i-1] * x1)

        x0 = x0n

    for j in range(num_samples):
        hd = 0
        for i in range(N):
            r0n =  twiddle[i] * hp[i]
            r1 = r0n - r0

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * r1

            r0 = r0n
            # Time shift the twiddle factors to the next time samples
            twiddle[i] *= rotate[i]

        sh[j] = conj(hd)
    return sh, hh
