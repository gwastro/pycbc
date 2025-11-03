/*
 * Copyright (C) 2016  Alex Nitz, Collin Capano
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*
 * Utilities for handling frequency compressed an unequally spaced frequency
 * domain waveforms.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Added for memset
#include <stdint.h> // Added for int64_t
#include <complex>


/* =====================================================================
 *
 * FAST LINEAR INTERPOLATION (HELPER)
 *
 * =====================================================================
 */

template <class T, class T_COMPLEX>
static inline void _decomp_ccode_segment(
    T_COMPLEX* h_seg,
    T f1, T f2, T a1, T a2, T p1, T p2,
    T df, int64_t k, const int64_t kmax,
    const int update_interval)
{
    T* outptr = (T*) h_seg;
    T inv_sdf = 1.0 / (f2 - f1);
    T m_amp = (a2 - a1) * inv_sdf;
    T b_amp = a1 - m_amp * f1;
    T m_phi = (p2 - p1) * inv_sdf;
    T b_phi = p1 - m_phi * f1;

    T h_re, h_im, g_re, g_im, dphi_re, dphi_im;
    T incrh_re, incrh_im, incrg_re, incrg_im, f;
    int64_t findex = k;
    int64_t k_sub_max;

    while (findex < kmax) {
        f = findex * df;
        T interp_amp = m_amp * f + b_amp;
        T interp_phi = m_phi * f + b_phi;
        dphi_re = cos(m_phi * df);
        dphi_im = sin(m_phi * df);
        h_re = interp_amp * cos(interp_phi);
        h_im = interp_amp * sin(interp_phi);
        g_re = m_amp * df * cos(interp_phi);
        g_im = m_amp * df * sin(interp_phi);

        *outptr = h_re;
        *(outptr + 1) = h_im;
        outptr += 2;
        findex++;

        k_sub_max = findex + update_interval;
        if (k_sub_max > kmax) k_sub_max = kmax;

        while (findex < k_sub_max) {
            incrh_re = h_re * dphi_re - h_im * dphi_im;
            incrh_im = h_re * dphi_im + h_im * dphi_re;
            incrg_re = g_re * dphi_re - g_im * dphi_im;
            incrg_im = g_re * dphi_im + g_im * dphi_re;
            h_re = incrh_re + incrg_re;
            h_im = incrh_im + incrg_im;
            g_re = incrg_re;
            g_im = incrg_im;

            *outptr = h_re;
            *(outptr + 1) = h_im;
            outptr += 2;
            findex++;
        }
    }
}


/* =====================================================================
 *
 * FAST QUADRATIC INTERPOLATION (HELPER)
 *
 * =====================================================================
 */

template <class T, class T_COMPLEX>
static inline void _decomp_qcode_segment(
    T_COMPLEX* h_seg,
    T f0, T f1, T f2, T a0, T a1, T a2, T p0, T p1, T p2,
    T df, int64_t k_start, const int64_t kmax, // Renamed k to k_start
    const int update_interval)
{
    T f;
    T h2 = df * df;

    // Denominators for Lagrange basis (constant for segment)
    T denom0 = (f0 - f1) * (f0 - f2);
    T denom1 = (f1 - f0) * (f1 - f2);
    T denom2 = (f2 - f0) * (f2 - f1);

    // --- Get power-basis coefficients: c2*f^2 + c1*f + c0 ---
    T c_a2 = a0/denom0 + a1/denom1 + a2/denom2;
    T c_p2 = p0/denom0 + p1/denom1 + p2/denom2;
    T c_a1 = -a0*(f1+f2)/denom0 - a1*(f0+f2)/denom1 - a2*(f0+f1)/denom2;
    T c_p1 = -p0*(f1+f2)/denom0 - p1*(f0+f2)/denom1 - p2*(f0+f1)/denom2;
    T c_a0 = a0*f1*f2/denom0 + a1*f0*f2/denom1 + a2*f0*f1/denom2;
    T c_p0 = p0*f1*f2/denom0 + p1*f0*f2/denom1 + p2*f0*f1/denom2;

    // --- Finite difference constants (constant for segment) ---
    T d2_a_const = 2 * c_a2 * h2;
    T d2_p_const = 2 * c_p2 * h2;
    T_COMPLEX d2_phase_const = std::polar((T)1.0, d2_p_const);

    // Stepper variables
    T a, p, d1_a, d1_p;
    T_COMPLEX phase, d_phase;
    int64_t k_sub_max;
    int64_t k = k_start; // k is now the local loop variable

    for (; k < kmax; ) {
        // --- Re-calculate steppers to correct FP error ---
        f = k * df;
        k_sub_max = k + update_interval;
        if (k_sub_max > kmax) k_sub_max = kmax;

        a = c_a2*f*f + c_a1*f + c_a0;
        p = c_p2*f*f + c_p1*f + c_p0;
        d1_a = c_a2*(2*f*df + h2) + c_a1*df;
        d1_p = c_p2*(2*f*df + h2) + c_p1*df;
        phase = std::polar((T)1.0, p);
        d_phase = std::polar((T)1.0, d1_p);

        // --- Fast Inner Loop ---
        for (; k < k_sub_max; k++) {
            h_seg[k - k_start] = a * phase; // Corrected indexing
            phase = phase * d_phase;
            d_phase = d_phase * d2_phase_const;
            a = a + d1_a;
            d1_a = d1_a + d2_a_const;
        }
    }
}


/* =====================================================================
 *
 * FAST CUBIC INTERPOLATION (HELPER)
 *
 * =====================================================================
 */

template <class T, class T_COMPLEX>
static inline void _decomp_tcode_segment(
    T_COMPLEX* h_seg,
    T f0, T f1, T f2, T f3, T a0, T a1, T a2, T a3, T p0, T p1, T p2, T p3,
    T df, int64_t k_start, const int64_t kmax, // Renamed k to k_start
    const int update_interval)
{
    T f;
    T h2 = df * df;
    T h3 = h2 * df;

    // --- Denominators for Lagrange basis (constant for segment) ---
    T denom0 = (f0-f1)*(f0-f2)*(f0-f3);
    T denom1 = (f1-f0)*(f1-f2)*(f1-f3);
    T denom2 = (f2-f0)*(f2-f1)*(f2-f3);
    T denom3 = (f3-f0)*(f3-f1)*(f3-f2);

    // --- Get power-basis coefficients: c3*f^3 + c2*f^2 + c1*f + c0 ---
    // --- Amplitude ---
    T c_a3 = a0/denom0 + a1/denom1 + a2/denom2 + a3/denom3;
    T c_a2 = -a0*(f1+f2+f3)/denom0 - a1*(f0+f2+f3)/denom1 - a2*(f0+f1+f3)/denom2 - a3*(f0+f1+f2)/denom3;
    T c_a1 = a0*(f1*f2 + f1*f3 + f2*f3)/denom0 + a1*(f0*f2 + f0*f3 + f2*f3)/denom1 + a2*(f0*f1 + f0*f3 + f1*f3)/denom2 + a3*(f0*f1 + f0*f2 + f1*f2)/denom3;
    T c_a0 = -a0*(f1*f2*f3)/denom0 - a1*(f0*f2*f3)/denom1 - a2*(f0*f1*f3)/denom2 - a3*(f0*f1*f2)/denom3;
    // --- Phase ---
    T c_p3 = p0/denom0 + p1/denom1 + p2/denom2 + p3/denom3;
    T c_p2 = -p0*(f1+f2+f3)/denom0 - p1*(f0+f2+f3)/denom1 - p2*(f0+f1+f3)/denom2 - p3*(f0+f1+f2)/denom3;
    T c_p1 = p0*(f1*f2 + f1*f3 + f2*f3)/denom0 + p1*(f0*f2 + f0*f3 + f2*f3)/denom1 + p2*(f0*f1 + f0*f3 + f1*f3)/denom2 + p3*(f0*f1 + f0*f2 + f1*f2)/denom3;
    T c_p0 = -p0*(f1*f2*f3)/denom0 - p1*(f0*f2*f3)/denom1 - p2*(f0*f1*f3)/denom2 - p3*(f0*f1*f2)/denom3;

    // --- Finite difference constants (constant for segment) ---
    T d3_a_const = 6 * c_a3 * h3;
    T d3_p_const = 6 * c_p3 * h3;
    T_COMPLEX d3_phase_const = std::polar((T)1.0, d3_p_const);

    // Stepper variables
    T a, p, d1_a, d1_p, d2_a, d2_p;
    T_COMPLEX phase, d1_phase, d2_phase;
    int64_t k_sub_max;
    int64_t k = k_start; // k is now the local loop variable

    for (; k < kmax; ) {
        // --- Re-calculate steppers to correct FP error ---
        f = k * df;
        k_sub_max = k + update_interval;
        if (k_sub_max > kmax) k_sub_max = kmax;

        // Initial values for this sub-block
        a = c_a3*f*f*f + c_a2*f*f + c_a1*f + c_a0;
        p = c_p3*f*f*f + c_p2*f*f + c_p1*f + c_p0;

        // First differences at f (Delta_f)
        d1_a = c_a3*(3*f*f*df + 3*f*h2 + h3) + c_a2*(2*f*df + h2) + c_a1*df;
        d1_p = c_p3*(3*f*f*df + 3*f*h2 + h3) + c_p2*(2*f*df + h2) + c_p1*df;

        // Second differences at f (Delta^2_f)
        d2_a = c_a3*(6*f*h2 + 6*h3) + c_a2*(2*h2);
        d2_p = c_p3*(6*f*h2 + 6*h3) + c_p2*(2*h2);

        // Complex phase steppers
        phase = std::polar((T)1.0, p);
        d1_phase = std::polar((T)1.0, d1_p);
        d2_phase = std::polar((T)1.0, d2_p);

        // --- Fast Inner Loop ---
        for (; k < k_sub_max; k++) {
            h_seg[k - k_start] = a * phase; // Corrected indexing

            // Step phase forward
            phase = phase * d1_phase;
            d1_phase = d1_phase * d2_phase;
            d2_phase = d2_phase * d3_phase_const;

            // Step amplitude forward
            a = a + d1_a;
            d1_a = d1_a + d2_a;
            d2_a = d2_a + d3_a_const;
        }
    }
}


/* =====================================================================
 *
 * FAST QUARTIC INTERPOLATION (HELPER) - NEWLY ADDED
 *
 * =====================================================================
 */

template <class T, class T_COMPLEX>
static inline void _decomp_Qcode_segment(
    T_COMPLEX* h_seg,
    T f0, T f1, T f2, T f3, T f4, // 5 freqs
    T a0, T a1, T a2, T a3, T a4, // 5 amps
    T p0, T p1, T p2, T p3, T p4, // 5 phases
    T df, int64_t k_start, const int64_t kmax,
    const int update_interval)
{
    T f;
    T h2 = df * df;
    T h3 = h2 * df;
    T h4 = h3 * df;

    // --- Denominators for Lagrange basis (constant for segment) ---
    T denom0 = (f0-f1)*(f0-f2)*(f0-f3)*(f0-f4);
    T denom1 = (f1-f0)*(f1-f2)*(f1-f3)*(f1-f4);
    T denom2 = (f2-f0)*(f2-f1)*(f2-f3)*(f2-f4);
    T denom3 = (f3-f0)*(f3-f1)*(f3-f2)*(f3-f4);
    T denom4 = (f4-f0)*(f4-f1)*(f4-f2)*(f4-f3);

    // --- Get power-basis coefficients: c4*f^4 + c3*f^3 + c2*f^2 + c1*f + c0 ---
    // --- Amplitude ---
    T c_a4 = a0/denom0 + a1/denom1 + a2/denom2 + a3/denom3 + a4/denom4;
    T c_a3 = -a0*(f1+f2+f3+f4)/denom0 - a1*(f0+f2+f3+f4)/denom1 - a2*(f0+f1+f3+f4)/denom2 - a3*(f0+f1+f2+f4)/denom3 - a4*(f0+f1+f2+f3)/denom4;
    T c_a2 = a0*(f1*f2+f1*f3+f1*f4+f2*f3+f2*f4+f3*f4)/denom0 + a1*(f0*f2+f0*f3+f0*f4+f2*f3+f2*f4+f3*f4)/denom1 + a2*(f0*f1+f0*f3+f0*f4+f1*f3+f1*f4+f3*f4)/denom2 + a3*(f0*f1+f0*f2+f0*f4+f1*f2+f1*f4+f2*f4)/denom3 + a4*(f0*f1+f0*f2+f0*f3+f1*f2+f1*f3+f2*f3)/denom4;
    T c_a1 = -a0*(f1*f2*f3+f1*f2*f4+f1*f3*f4+f2*f3*f4)/denom0 - a1*(f0*f2*f3+f0*f2*f4+f0*f3*f4+f2*f3*f4)/denom1 - a2*(f0*f1*f3+f0*f1*f4+f0*f3*f4+f1*f3*f4)/denom2 - a3*(f0*f1*f2+f0*f1*f4+f0*f2*f4+f1*f2*f4)/denom3 - a4*(f0*f1*f2+f0*f1*f3+f0*f2*f3+f1*f2*f3)/denom4;
    T c_a0 = a0*(f1*f2*f3*f4)/denom0 + a1*(f0*f2*f3*f4)/denom1 + a2*(f0*f1*f3*f4)/denom2 + a3*(f0*f1*f2*f4)/denom3 + a4*(f0*f1*f2*f3)/denom4;
    // --- Phase ---
    T c_p4 = p0/denom0 + p1/denom1 + p2/denom2 + p3/denom3 + p4/denom4;
    T c_p3 = -p0*(f1+f2+f3+f4)/denom0 - p1*(f0+f2+f3+f4)/denom1 - p2*(f0+f1+f3+f4)/denom2 - p3*(f0+f1+f2+f4)/denom3 - p4*(f0+f1+f2+f3)/denom4;
    T c_p2 = p0*(f1*f2+f1*f3+f1*f4+f2*f3+f2*f4+f3*f4)/denom0 + p1*(f0*f2+f0*f3+f0*f4+f2*f3+f2*f4+f3*f4)/denom1 + p2*(f0*f1+f0*f3+f0*f4+f1*f3+f1*f4+f3*f4)/denom2 + p3*(f0*f1+f0*f2+f0*f4+f1*f2+f1*f4+f2*f4)/denom3 + p4*(f0*f1+f0*f2+f0*f3+f1*f2+f1*f3+f2*f3)/denom4;
    T c_p1 = -p0*(f1*f2*f3+f1*f2*f4+f1*f3*f4+f2*f3*f4)/denom0 - p1*(f0*f2*f3+f0*f2*f4+f0*f3*f4+f2*f3*f4)/denom1 - p2*(f0*f1*f3+f0*f1*f4+f0*f3*f4+f1*f3*f4)/denom2 - p3*(f0*f1*f2+f0*f1*f4+f0*f2*f4+f1*f2*f4)/denom3 - p4*(f0*f1*f2+f0*f1*f3+f0*f2*f3+f1*f2*f3)/denom4;
    T c_p0 = p0*(f1*f2*f3*f4)/denom0 + p1*(f0*f2*f3*f4)/denom1 + p2*(f0*f1*f3*f4)/denom2 + p3*(f0*f1*f2*f4)/denom3 + p4*(f0*f1*f2*f3)/denom4;

    // --- Finite difference constants (constant for segment) ---
    T d4_a_const = 24 * c_a4 * h4;
    T d4_p_const = 24 * c_p4 * h4;
    T_COMPLEX d4_phase_const = std::polar((T)1.0, d4_p_const);

    // Stepper variables
    T a, p, d1_a, d1_p, d2_a, d2_p, d3_a, d3_p;
    T_COMPLEX phase, d1_phase, d2_phase, d3_phase;
    int64_t k_sub_max;
    int64_t k = k_start; // k is now the local loop variable

    for (; k < kmax; ) {
        // --- Re-calculate steppers to correct FP error ---
        f = k * df;
        k_sub_max = k + update_interval;
        if (k_sub_max > kmax) k_sub_max = kmax;

        // Initial values for this sub-block
        a = c_a4*f*f*f*f + c_a3*f*f*f + c_a2*f*f + c_a1*f + c_a0;
        p = c_p4*f*f*f*f + c_p3*f*f*f + c_p2*f*f + c_p1*f + c_p0;

        // First differences at f (Delta_f)
        d1_a = c_a4*(4*f*f*f*df + 6*f*f*h2 + 4*f*h3 + h4) + c_a3*(3*f*f*df + 3*f*h2 + h3) + c_a2*(2*f*df + h2) + c_a1*df;
        d1_p = c_p4*(4*f*f*f*df + 6*f*f*h2 + 4*f*h3 + h4) + c_p3*(3*f*f*df + 3*f*h2 + h3) + c_p2*(2*f*df + h2) + c_p1*df;

        // Second differences at f (Delta^2_f)
        d2_a = c_a4*(12*f*f*h2 + 24*f*h3 + 14*h4) + c_a3*(6*f*h2 + 6*h3) + c_a2*(2*h2);
        d2_p = c_p4*(12*f*f*h2 + 24*f*h3 + 14*h4) + c_p3*(6*f*h2 + 6*h3) + c_p2*(2*h2);
        
        // Third differences at f (Delta^3_f)
        d3_a = c_a4*(24*f*h3 + 36*h4) + c_a3*(6*h3);
        d3_p = c_p4*(24*f*h3 + 36*h4) + c_p3*(6*h3);

        // Complex phase steppers
        phase = std::polar((T)1.0, p);
        d1_phase = std::polar((T)1.0, d1_p);
        d2_phase = std::polar((T)1.0, d2_p);
        d3_phase = std::polar((T)1.0, d3_p);

        // --- Fast Inner Loop ---
        for (; k < k_sub_max; k++) {
            h_seg[k - k_start] = a * phase; // Corrected indexing

            // Step phase forward
            phase = phase * d1_phase;
            d1_phase = d1_phase * d2_phase;
            d2_phase = d2_phase * d3_phase;
            d3_phase = d3_phase * d4_phase_const;

            // Step amplitude forward
            a = a + d1_a;
            d1_a = d1_a + d2_a;
            d2_a = d2_a + d3_a;
            d3_a = d3_a + d4_a_const;
        }
    }
}


/* =====================================================================
 *
 * ORIGINAL LINEAR INTERPOLATOR (MAIN)
 *
 * =====================================================================
 */

void _decomp_ccode_double(std::complex<double> * h,
                          double delta_f,
                          const int64_t hlen,
                          const int64_t start_index,
                          double * sample_frequencies,
                          double * amp,
                          double * phase,
                          const int64_t sflen,
                          const int64_t imin)
{
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    double f1, f2, a1, a2, p1, p2;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<double>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // Call the fast linear helper
        _decomp_ccode_segment<double, std::complex<double>>(
            &h[k], f1, f2, a1, a2, p1, p2,
            delta_f, k, kmax, update_interval);

        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<double>)*(hlen-last_findex));
    }
}

void _decomp_ccode_float(std::complex<float> * h,
                         float delta_f,
                         const int64_t hlen,
                         const int64_t start_index,
                         float * sample_frequencies,
                         float * amp,
                         float * phase,
                         const int64_t sflen,
                         const int64_t imin)
{
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    float f1, f2, a1, a2, p1, p2;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<float>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // Call the fast linear helper
        _decomp_ccode_segment<float, std::complex<float>>(
            &h[k], f1, f2, a1, a2, p1, p2,
            delta_f, k, kmax, update_interval);

        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<float>)*(hlen-last_findex));
    }
}


/* =====================================================================
 *
 * FAST QUADRATIC INTERPOLATOR (MAIN)
 *
 * =====================================================================
 */

void _decomp_qcode_double(std::complex<double> * h,
                          double delta_f,
                          const int64_t hlen,
                          const int64_t start_index,
                          double * sample_frequencies,
                          double * amp,
                          double * phase,
                          const int64_t sflen,
                          const int64_t imin)
{
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    double f0, f1, f2, a0, a1, a2, p0, p1, p2;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<double>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // === Use robust linear interpolation for the first segment ===
        if (i == imin) {
            _decomp_ccode_segment<double, std::complex<double>>(
                &h[k], f1, f2, a1, a2, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Use fast quadratic interpolation for all subsequent segments ===
        else {
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];

            _decomp_qcode_segment<double, std::complex<double>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<double>)*(hlen-last_findex));
    }
}

void _decomp_qcode_float(std::complex<float> * h,
                         float delta_f,
                         const int64_t hlen,
                         const int64_t start_index,
                         float * sample_frequencies,
                         float * amp,
                         float * phase,
                         const int64_t sflen,
                         const int64_t imin)
{
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    float f0, f1, f2, a0, a1, a2, p0, p1, p2;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<float>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // === Use robust linear interpolation for the first segment ===
        if (i == imin) {
            _decomp_ccode_segment<float, std::complex<float>>(
                &h[k], f1, f2, a1, a2, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Use fast quadratic interpolation for all subsequent segments ===
        else {
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];

            _decomp_qcode_segment<float, std::complex<float>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<float>)*(hlen-last_findex));
    }
}


/* =====================================================================
 *
 * FAST CUBIC INTERPOLATOR (MAIN)
 *
 * =====================================================================
 */

void _decomp_tcode_double(std::complex<double> * h,
                          double delta_f,
                          const int64_t hlen,
                          const int64_t start_index,
                          double * sample_frequencies,
                          double * amp,
                          double * phase,
                          const int64_t sflen,
                          const int64_t imin)
{
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    double f0, f1, f2, f3, a0, a1, a2, a3, p0, p1, p2, p3;

    // A cubic interpolator requires at least 3 segments
    // (4 points) to run. If we have less, we must degrade.
    // sflen=2 (1 seg) -> linear only
    // sflen=3 (2 segs) -> linear, quadratic
    // sflen=4 (3 segs) -> linear, quadratic, quadratic
    // sflen>=5 (>=4 segs) -> linear, quadratic, cubic..., quadratic

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<double>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // === Boundary 1: First segment (i == imin) ===
        // Must use linear (points i, i+1)
        if (i == imin) {
            _decomp_ccode_segment<double, std::complex<double>>(
                &h[k], f1, f2, a1, a2, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 2: Last segment (i == sflen - 2) or not enough points ===
        // Must use quadratic (points i-1, i, i+1)
        else if (i == sflen - 2 || sflen < 5) { // Degrade if not enough points for cubic
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];
            _decomp_qcode_segment<double, std::complex<double>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === All other segments: Use fast cubic (points i-1, i, i+1, i+2) ===
        else {
            f0 = sample_frequencies[i-1];
            f3 = sample_frequencies[i+2];
            a0 = amp[i-1];
            a3 = amp[i+2];
            p0 = phase[i-1];
            p3 = phase[i+2];

            _decomp_tcode_segment<double, std::complex<double>>(
                &h[k], f0, f1, f2, f3, a0, a1, a2, a3, p0, p1, p2, p3,
                delta_f, k, kmax, update_interval);
        }
        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<double>)*(hlen-last_findex));
    }
}


void _decomp_tcode_float(std::complex<float> * h,
                         float delta_f,
                         const int64_t hlen,
                         const int64_t start_index,
                         float * sample_frequencies,
                         float * amp,
                         float * phase,
                         const int64_t sflen,
                         const int64_t imin)
{
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    float f0, f1, f2, f3, a0, a1, a2, a3, p0, p1, p2, p3;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<float>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // === Boundary 1: First segment (i == imin) ===
        if (i == imin) {
            _decomp_ccode_segment<float, std::complex<float>>(
                &h[k], f1, f2, a1, a2, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 2: Last segment (i == sflen - 2) or not enough points ===
        else if (i == sflen - 2 || sflen < 5) { // Degrade if not enough points for cubic
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];
            _decomp_qcode_segment<float, std::complex<float>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === All other segments: Use fast cubic (points i-1, i, i+1, i+2) ===
        else {
            f0 = sample_frequencies[i-1];
            f3 = sample_frequencies[i+2];
            a0 = amp[i-1];
            a3 = amp[i+2];
            p0 = phase[i-1];
            p3 = phase[i+2];

            _decomp_tcode_segment<float, std::complex<float>>(
                &h[k], f0, f1, f2, f3, a0, a1, a2, a3, p0, p1, p2, p3,
                delta_f, k, kmax, update_interval);
        }
        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<float>)*(hlen-last_findex));
    }
}


/* =====================================================================
 *
 * FAST QUARTIC INTERPOLATOR (MAIN) - NEWLY ADDED
 *
 * =====================================================================
 */

// Main function for quartic interpolation
void _decomp_Qcode_double(std::complex<double> * h,
                          double delta_f,
                          const int64_t hlen,
                          const int64_t start_index,
                          double * sample_frequencies,
                          double * amp,
                          double * phase,
                          const int64_t sflen,
                          const int64_t imin)
{
    // Quartic requires 5 points (4 segments).
    // If fewer than 6 points (5 segments) are available,
    // degrade to the cubic function which can handle all smaller cases.
    if (sflen < 6) {
        _decomp_tcode_double(h, delta_f, hlen, start_index,
                             sample_frequencies, amp, phase, sflen, imin);
        return;
    }
    
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    double f0, f1, f2, f3, f4, a0, a1, a2, a3, a4, p0, p1, p2, p3, p4;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<double>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        // Get segment boundaries
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        // Calculate start and end indices for this segment
        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // === Boundary 1: First segment (i == imin) ===
        if (i == imin) {
            _decomp_ccode_segment<double, std::complex<double>>(
                &h[k], f1, f2, a1, a2, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 2: Second segment (i == imin + 1) ===
        else if (i == imin + 1) {
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];
            _decomp_qcode_segment<double, std::complex<double>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 3: Last segment (i == sflen - 2) ===
        else if (i == sflen - 2) {
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];
            _decomp_qcode_segment<double, std::complex<double>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 4: Second-to-last segment (i == sflen - 3) ===
        else if (i == sflen - 3) {
            f0 = sample_frequencies[i-1];
            f3 = sample_frequencies[i+2];
            a0 = amp[i-1];
            a3 = amp[i+2];
            p0 = phase[i-1];
            p3 = phase[i+2];
            _decomp_tcode_segment<double, std::complex<double>>(
                &h[k], f0, f1, f2, f3, a0, a1, a2, a3, p0, p1, p2, p3,
                delta_f, k, kmax, update_interval);
        }
        // === All other segments: Use fast quartic (points i-1, i, i+1, i+2, i+3) ===
        else {
            f0 = sample_frequencies[i-1];
            f3 = sample_frequencies[i+2];
            f4 = sample_frequencies[i+3];
            a0 = amp[i-1];
            a3 = amp[i+2];
            a4 = amp[i+3];
            p0 = phase[i-1];
            p3 = phase[i+2];
            p4 = phase[i+3];

            _decomp_Qcode_segment<double, std::complex<double>>(
                &h[k], f0, f1, f2, f3, f4, a0, a1, a2, a3, a4, p0, p1, p2, p3, p4,
                delta_f, k, kmax, update_interval);
        }
        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<double>)*(hlen-last_findex));
    }
}

// Main function for quartic interpolation (float)
void _decomp_Qcode_float(std::complex<float> * h,
                         float delta_f,
                         const int64_t hlen,
                         const int64_t start_index,
                         float * sample_frequencies,
                         float * amp,
                         float * phase,
                         const int64_t sflen,
                         const int64_t imin)
{
    // Degrade to cubic if not enough points
    if (sflen < 6) {
        _decomp_tcode_float(h, delta_f, hlen, start_index,
                            sample_frequencies, amp, phase, sflen, imin);
        return;
    }
    
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    float f0, f1, f2, f3, f4, a0, a1, a2, a3, a4, p0, p1, p2, p3, p4;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<float>)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        f1 = sample_frequencies[i];
        f2 = sample_frequencies[i+1];
        a1 = amp[i];
        a2 = amp[i+1];
        p1 = phase[i];
        p2 = phase[i+1];

        if (i == imin) {
            k = start_index;
        } else {
            k = (int64_t)ceil(f1 / delta_f);
        }

        if (i == sflen - 2) {
           kmax = (int64_t)(f2 / delta_f) + 1;
        } else {
           kmax = (int64_t)(f2 / delta_f);
        }
        if (kmax > hlen) kmax = hlen;

        // === Boundary 1: First segment (i == imin) ===
        if (i == imin) {
            _decomp_ccode_segment<float, std::complex<float>>(
                &h[k], f1, f2, a1, a2, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 2: Second segment (i == imin + 1) ===
        else if (i == imin + 1) {
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];
            _decomp_qcode_segment<float, std::complex<float>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 3: Last segment (i == sflen - 2) ===
        else if (i == sflen - 2) {
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];
            _decomp_qcode_segment<float, std::complex<float>>(
                &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                delta_f, k, kmax, update_interval);
        }
        // === Boundary 4: Second-to-last segment (i == sflen - 3) ===
        else if (i == sflen - 3) {
            f0 = sample_frequencies[i-1];
            f3 = sample_frequencies[i+2];
            a0 = amp[i-1];
            a3 = amp[i+2];
            p0 = phase[i-1];
            p3 = phase[i+2];
            _decomp_tcode_segment<float, std::complex<float>>(
                &h[k], f0, f1, f2, f3, a0, a1, a2, a3, p0, p1, p2, p3,
                delta_f, k, kmax, update_interval);
        }
        // === All other segments: Use fast quartic (points i-1, i, i+1, i+2, i+3) ===
        else {
            f0 = sample_frequencies[i-1];
            f3 = sample_frequencies[i+2];
            f4 = sample_frequencies[i+3];
            a0 = amp[i-1];
            a3 = amp[i+2];
            a4 = amp[i+3];
            p0 = phase[i-1];
            p3 = phase[i+2];
            p4 = phase[i+3];

            _decomp_Qcode_segment<float, std::complex<float>>(
                &h[k], f0, f1, f2, f3, f4, a0, a1, a2, a3, a4, p0, p1, p2, p3, p4,
                delta_f, k, kmax, update_interval);
        }
        last_findex = kmax;
    }
    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<float>)*(hlen-last_findex));
    }
}
