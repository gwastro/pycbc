""" Optimized inner loop functions for the relative likelihood model
"""
cdef extern from "complex.h":
    double complex exp(double complex)
    double norm(double complex)
    double complex conj(double complex)
    double real(double complex)

cimport cython    

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
    cdef double complex hd=0, r0, r0n, r1
    cdef double hh=0
    
    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(-2.0j * 3.141592653 * dtc * freqs[i]) 
               * (fp * hp[i] + fc * hc[i])) / h00[i]        
        r1 = r0n - r0
        
        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1
            hh += b0[i-1] * norm(r0) + 2.0 * b1[i-1] * real(r1 * conj(r0))
    
        r0 = r0n
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
    cdef double complex hd=0, r0, r0n, r1
    cdef double hh=0
    
    N = freqs.shape[0]
    for i in range(N):
        r0n = (exp(-2.0j * 3.141592653 * dtc[i] * freqs[i])
               * (fp[i] * hp[i] + fc[i] * hc[i])) / h00[i]        
        r1 = r0n - r0
        
        if i > 0:
            hd += a0[i-1] * r0 + a1[i-1] * r1
            hh += b0[i-1] * norm(r0) + 2.0 * b1[i-1] * real(r1 * conj(r0))
    
        r0 = r0n
    return hd, hh
