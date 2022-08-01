""" Optimized inner loop functions for the relative likelihood model
"""
cdef extern from "complex.h":
    double complex exp(double complex)
    double norm(double complex)
    double complex conj(double complex)
    double real(double complex)

cimport cython

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
    return hd, hh

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
    return hd, hh
