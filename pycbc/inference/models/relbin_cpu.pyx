""" Optimized inner loop functions for the relative likelihood model
"""
cdef extern from "complex.h":
    double complex exp(double complex)
    double norm(double complex)
    double complex conj(double complex)
    double real(double complex)
    double imag(double complex)

# The phase pass below wants the real sine and cosine on their own, rather
# than the complex exponential declared above, so that it can write them
# straight into a pair of real buffers.
cdef extern from "math.h":
    double cos(double)
    double sin(double)

import numpy
cimport cython, numpy

# The pi here is deliberately the truncated 3.141592653 rather than M_PI.
# It appears only in the 2*pi of a time shift, so its whole effect is an
# exact rescaling dtc -> dtc*(1 - 1.877e-10): <h|h> is untouched, and on a
# real event the peak moves by 1.7e-12 s against a posterior width of
# 6.8e-5.  Replacing it would move published numbers for no correctness
# gain, so it stays, named once instead of spelled out twenty-one times.
cdef double MINUS_TWO_PI = -2.0 * 3.141592653
cdef double complex MINUS_TWO_PI_J = -2.0j * 3.141592653

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
        r0n = (exp(MINUS_TWO_PI_J * dtc * freqs[i])
               * (fp * hp[i] + fc * hc[i])) / h00[i]
        r0n *= conj((exp(MINUS_TWO_PI_J * dtc2 * freqs[i])
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
        r0n = (exp(MINUS_TWO_PI_J * dtc[i] * freqs[i])
               * (fp[i] * hp[i] + fc[i] * hc[i])) / h00[i]
        r0n *= conj((exp(MINUS_TWO_PI_J * dtc2[i] * freqs[i])
               * (fp2[i] * hp2[i] + fc2[i] * hc2[i])) / h002[i])
        r1 = r0n - r0

        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1

        r0 = r0n
    return hd

# Used for calculating the cross terms of
# two signals when analyzing multiple signals
# with no antenna response applied
cpdef likelihood_parts_det_multi(double [::1] freqs,
                     double dtc,
                     double complex[::1] hp,
                     double complex[::1] h00,
                     double dtc2,
                     double complex[::1] hp2,
                     double complex[::1] h002,
                     double complex[::1] a0,
                     double complex[::1] a1,
                     ) :
    cdef size_t i
    cdef double complex hd=0, r0, r0n, r1
    cdef int N

    N = freqs.shape[0]
    for i in range(N):
        r0n = conj((exp(MINUS_TWO_PI_J * dtc * freqs[i])
               * (hp[i])) / h00[i])
        r0n *= ((exp(MINUS_TWO_PI_J * dtc2 * freqs[i])
               * (hp2[i])) / h002[i])
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
        r0n = (exp(MINUS_TWO_PI_J * dtc * freqs[i])
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
    cdef double complex hd=0, r0, r0n, r1, x0=0, x1, x0n;
    cdef double hh=0
    cdef int N

    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(MINUS_TWO_PI_J * dtc * freqs[i])
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
        r0n = (exp(MINUS_TWO_PI_J * dtc[i] * freqs[i])
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
    cdef size_t i, j
    cdef double complex hd, r0, r0n, tw
    cdef double app, acc, apc, cp, sp, ph
    cdef int N = freqs.shape[0]
    cdef int num_samples = pol_phase.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)
    cdef double complex[::1] ar = numpy.empty(N, dtype=numpy.complex128)
    cdef double complex[::1] br = numpy.empty(N, dtype=numpy.complex128)

    for i in range(N):
        ph = MINUS_TWO_PI * dtc[i] * freqs[i]
        tw = (cos(ph) + 1.0j * sin(ph)) / h00[i]
        ar[i] = tw * (fp[i] * hp[i] + fc[i] * hc[i])
        br[i] = tw * (fp[i] * hc[i] - fc[i] * hp[i])
    hh_coefficients(ar, br, b0, b1, N, &app, &acc, &apc)

    for j in range(num_samples):
        cp = real(pol_phase[j])
        sp = imag(pol_phase[j])

        hd = 0
        for i in range(N):
            r0n = cp * ar[i] + sp * br[i]

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * (r0n - r0)

            r0 = r0n
        hdv[j] = conj(hd)
        hhv[j] = cp * cp * app + sp * sp * acc + cp * sp * apc
    return hdv, hhv
# Used where the antenna response may be frequency varying
# and there is a polarization vector marginalization
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
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
    cdef size_t i, j
    cdef double complex hd, r0, r0n
    cdef double hh=0, x0=0, x0n, tw, ph
    cdef int N = freqs.shape[0]
    cdef int num_samples = dtc.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)
    cdef double complex[::1] hr = numpy.empty(N, dtype=numpy.complex128)
    cdef double [::1] sphase = numpy.empty(N, dtype=numpy.float64)
    cdef double [::1] cphase = numpy.empty(N, dtype=numpy.float64)

    for i in range(N):
        ph = MINUS_TWO_PI * times[i] * freqs[i]
        hr[i] = ((cos(ph) + 1.0j * sin(ph))
                 * (fp[i] * hp[i] + fc[i] * hc[i])) / h00[i]

        x0n = norm(hr[i])
        if i > 0:
            hh += b0[i-1] * x0 + b1[i-1] * (x0n - x0)

        x0 = x0n

    for j in range(num_samples):
        tw = MINUS_TWO_PI * dtc[j]
        for i in range(N):
            sphase[i] = sin(tw * freqs[i])
            cphase[i] = cos(tw * freqs[i])

        hd = 0
        for i in range(N):
            r0n = (cphase[i] + 1.0j * sphase[i]) * hr[i]

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * (r0n - r0)

            r0 = r0n
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv
# Used where the antenna response may be frequency varying
# and there is a polarization vector marginalization
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
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
    cdef size_t i, j
    cdef double complex hd, r0, r0n, tw
    cdef double app, acc, apc, cp, sp, ph, twj
    cdef int N = freqs.shape[0]
    cdef int num_samples = pol_phase.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)
    cdef double complex[::1] ar = numpy.empty(N, dtype=numpy.complex128)
    cdef double complex[::1] br = numpy.empty(N, dtype=numpy.complex128)
    cdef double [::1] sphase = numpy.empty(N, dtype=numpy.float64)
    cdef double [::1] cphase = numpy.empty(N, dtype=numpy.float64)

    for i in range(N):
        ph = MINUS_TWO_PI * times[i] * freqs[i]
        tw = (cos(ph) + 1.0j * sin(ph)) / h00[i]
        ar[i] = tw * (fp[i] * hp[i] + fc[i] * hc[i])
        br[i] = tw * (fp[i] * hc[i] - fc[i] * hp[i])
    hh_coefficients(ar, br, b0, b1, N, &app, &acc, &apc)

    for j in range(num_samples):
        cp = real(pol_phase[j])
        sp = imag(pol_phase[j])
        twj = MINUS_TWO_PI * dtc[j]
        for i in range(N):
            sphase[i] = sin(twj * freqs[i])
            cphase[i] = cos(twj * freqs[i])

        hd = 0
        for i in range(N):
            r0n = ((cphase[i] + 1.0j * sphase[i])
                   * (cp * ar[i] + sp * br[i]))

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * (r0n - r0)

            r0 = r0n
        hdv[j] = conj(hd)
        hhv[j] = cp * cp * app + sp * sp * acc + cp * sp * apc
    return hdv, hhv
# The <h|h> of a whole bank of samples, without a sample loop.
#
# Every sample scales the two polarizations by a real pair (fp, fc) and then
# applies a time shift, and the time shift has unit modulus, so it drops out
# of |fp * hp + fc * hc|^2 entirely.  What is left is a quadratic form in
# (fp, fc), and summing its three coefficient series once with the relative
# binning weights gives constants from which any sample's <h|h> follows as
# fp * fp * app + fc * fc * acc + fp * fc * apc.  Both the factor of two in
# the cross term and its being the real part are needed; neither is optional.
# This is the same reasoning snr_predictor already uses to lift its own
# <h|h> out of the sample loop.
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)     # Disable checking for dividing by zero
cdef void hh_coefficients(double complex[::1] hpr,
                     double complex[::1] hcr,
                     double [::1] b0,
                     double [::1] b1,
                     int N,
                     double *app,
                     double *acc,
                     double *apc,
                     ) :
    cdef size_t i
    cdef double p0=0, c0=0, x0=0, p0n, c0n, x0n
    cdef double sp=0, sc=0, sx=0

    for i in range(N):
        p0n = norm(hpr[i])
        c0n = norm(hcr[i])
        x0n = 2.0 * (real(hpr[i]) * real(hcr[i])
                     + imag(hpr[i]) * imag(hcr[i]))

        if i > 0:
            sp += b0[i-1] * p0 + b1[i-1] * (p0n - p0)
            sc += b0[i-1] * c0 + b1[i-1] * (c0n - c0)
            sx += b0[i-1] * x0 + b1[i-1] * (x0n - x0)

        p0 = p0n
        c0 = c0n
        x0 = x0n
    app[0] = sp
    acc[0] = sc
    apc[0] = sx

# Standard likelihood but simultaneously handling multiple sky or time points
#
# Three things are kept out of the inner loop that do not belong there.  The
# division by h00 depends on the frequency bin alone, so it is done once for
# the bank instead of once per (bin, sample).  <h|h> is not accumulated here
# at all, it comes from the constants above.  And the phase is evaluated in
# a pass of its own, into scratch buffers; nothing in that pass depends on
# the previous element, so the compiler is free to vectorize the sine and
# cosine, which it cannot do while the serial r0 carry shares the loop.
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
    cdef size_t i, j
    cdef double complex hd, r0, r0n
    cdef double app, acc, apc, fpj, fcj, tw
    cdef int N = freqs.shape[0]
    cdef int num_samples = fp.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)
    cdef double complex[::1] hpr = numpy.empty(N, dtype=numpy.complex128)
    cdef double complex[::1] hcr = numpy.empty(N, dtype=numpy.complex128)
    cdef double [::1] sphase = numpy.empty(N, dtype=numpy.float64)
    cdef double [::1] cphase = numpy.empty(N, dtype=numpy.float64)

    for i in range(N):
        hpr[i] = hp[i] / h00[i]
        hcr[i] = hc[i] / h00[i]
    hh_coefficients(hpr, hcr, b0, b1, N, &app, &acc, &apc)

    for j in range(num_samples):
        fpj = fp[j]
        fcj = fc[j]
        tw = MINUS_TWO_PI * dtc[j]
        for i in range(N):
            sphase[i] = sin(tw * freqs[i])
            cphase[i] = cos(tw * freqs[i])

        hd = 0
        for i in range(N):
            r0n = ((cphase[i] + 1.0j * sphase[i])
                   * (fpj * hpr[i] + fcj * hcr[i]))

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * (r0n - r0)

            r0 = r0n
        hdv[j] = conj(hd)
        hhv[j] = fpj * fpj * app + fcj * fcj * acc + fpj * fcj * apc
    return hdv, hhv

# Standard likelihood but simultaneously handling multiple time points
#
# Here fp and fc are scalars, so every sample sees the same waveform ratio
# and differs only by a time shift.  A time shift cannot change <h|h>, so
# there is one value of it for the whole bank and it is computed once, with
# the ratio, before the sample loop is entered.  The phase again gets a pass
# of its own so that the sine and cosine can vectorize.
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
    cdef size_t i, j
    cdef double complex hd, r0, r0n
    cdef double hh=0, x0=0, x0n, tw
    cdef int N = freqs.shape[0]
    cdef int num_samples = dtc.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)
    cdef double complex[::1] hr = numpy.empty(N, dtype=numpy.complex128)
    cdef double [::1] sphase = numpy.empty(N, dtype=numpy.float64)
    cdef double [::1] cphase = numpy.empty(N, dtype=numpy.float64)

    for i in range(N):
        hr[i] = (fp * hp[i] + fc * hc[i]) / h00[i]

        x0n = norm(hr[i])
        if i > 0:
            hh += b0[i-1] * x0 + b1[i-1] * (x0n - x0)

        x0 = x0n

    for j in range(num_samples):
        tw = MINUS_TWO_PI * dtc[j]
        for i in range(N):
            sphase[i] = sin(tw * freqs[i])
            cphase[i] = cos(tw * freqs[i])

        hd = 0
        for i in range(N):
            r0n = (cphase[i] + 1.0j * sphase[i]) * hr[i]

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * (r0n - r0)

            r0 = r0n
        hdv[j] = conj(hd)
        hhv[j] = hh
    return hdv, hhv

# Like above, but if only polarization is marginalized
#
# Here dtc is a scalar, so the time shift is the same for every sample and
# there is no transcendental left to evaluate inside the sample loop at all:
# the phase is folded into the per bin waveform ratios once.  Folding it in
# does not disturb the <h|h> constants, since a unit modulus factor common
# to both polarizations cancels out of every term of the quadratic form.
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
    cdef size_t i, j
    cdef double complex hd, r0, r0n, tw
    cdef double app, acc, apc, fpj, fcj
    cdef int N = freqs.shape[0]
    cdef int num_samples = fp.shape[0]

    cdef numpy.ndarray[numpy.complex128_t, ndim=1] hdv = numpy.empty(num_samples, dtype=numpy.complex128)
    cdef numpy.ndarray[numpy.float64_t, ndim=1] hhv = numpy.empty(num_samples, dtype=numpy.float64)
    cdef double complex[::1] hpr = numpy.empty(N, dtype=numpy.complex128)
    cdef double complex[::1] hcr = numpy.empty(N, dtype=numpy.complex128)

    for i in range(N):
        tw = exp(MINUS_TWO_PI_J * dtc * freqs[i]) / h00[i]
        hpr[i] = tw * hp[i]
        hcr[i] = tw * hc[i]
    hh_coefficients(hpr, hcr, b0, b1, N, &app, &acc, &apc)

    for j in range(num_samples):
        fpj = fp[j]
        fcj = fc[j]

        hd = 0
        for i in range(N):
            r0n = fpj * hpr[i] + fcj * hcr[i]

            if i > 0:
                hd += a0[i-1] * r0 + a1[i-1] * (r0n - r0)

            r0 = r0n
        hdv[j] = conj(hd)
        hhv[j] = fpj * fpj * app + fcj * fcj * acc + fpj * fcj * apc
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
        twiddle[i] = exp(MINUS_TWO_PI_J * tstart * freqs[i]) / h00[i]
        rotate[i] =  exp(MINUS_TWO_PI_J * delta_t * freqs[i])

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
        twiddle[i] = exp(MINUS_TWO_PI_J * tstart * freqs[i]) / h00[i]
        rotate[i] =  exp(MINUS_TWO_PI_J * delta_t * freqs[i])

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
