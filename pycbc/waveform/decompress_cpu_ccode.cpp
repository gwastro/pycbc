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
    double f1, double f2, double a1, double a2, double p1, double p2,
    double df, int64_t k, const int64_t kmax,
    const int update_interval)
{
    T* outptr = (T*) h_seg;
    double inv_sdf = 1.0 / (f2 - f1);
    double m_amp = (a2 - a1) * inv_sdf;
    double b_amp = a1 - m_amp * f1;
    double m_phi = (p2 - p1) * inv_sdf;
    double b_phi = p1 - m_phi * f1;

    double h_re, h_im, g_re, g_im, dphi_re, dphi_im;
    double incrh_re, incrh_im, incrg_re, incrg_im, f;
    int64_t findex = k;
    int64_t k_sub_max;

    while (findex < kmax) {
        f = findex * df;
        double interp_amp = m_amp * f + b_amp;
        double interp_phi = m_phi * f + b_phi;
        dphi_re = cos(m_phi * df);
        dphi_im = sin(m_phi * df);
        h_re = interp_amp * cos(interp_phi);
        h_im = interp_amp * sin(interp_phi);
        g_re = m_amp * df * cos(interp_phi);
        g_im = m_amp * df * sin(interp_phi);

        *outptr = (T)h_re;
        *(outptr + 1) = (T)h_im;
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

            *outptr = (T)h_re;
            *(outptr + 1) = (T)h_im;
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
    double f0, double f1, double f2, double a0, double a1, double a2, double p0, double p1, double p2,
    double df, int64_t k_start, const int64_t kmax, // Renamed k to k_start
    const int update_interval)
{
    double f;
    double h2 = df * df;

    // Denominators for Lagrange basis (constant for segment)
    double denom0 = (f0 - f1) * (f0 - f2);
    double denom1 = (f1 - f0) * (f1 - f2);
    double denom2 = (f2 - f0) * (f2 - f1);

    // --- Get power-basis coefficients: c2*f^2 + c1*f + c0 ---
    double c_a2 = a0/denom0 + a1/denom1 + a2/denom2;
    double c_p2 = p0/denom0 + p1/denom1 + p2/denom2;
    double c_a1 = -a0*(f1+f2)/denom0 - a1*(f0+f2)/denom1 - a2*(f0+f1)/denom2;
    double c_p1 = -p0*(f1+f2)/denom0 - p1*(f0+f2)/denom1 - p2*(f0+f1)/denom2;
    double c_a0 = a0*f1*f2/denom0 + a1*f0*f2/denom1 + a2*f0*f1/denom2;
    double c_p0 = p0*f1*f2/denom0 + p1*f0*f2/denom1 + p2*f0*f1/denom2;

    // --- Finite difference constants (constant for segment) ---
    double d2_a_const = 2 * c_a2 * h2;
    double d2_p_const = 2 * c_p2 * h2;
    std::complex<double> d2_phase_const = std::polar(1.0, d2_p_const);

    // Stepper variables
    double a, p, d1_a, d1_p;
    std::complex<double> phase, d_phase;
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
        phase = std::polar(1.0, p);
        d_phase = std::polar(1.0, d1_p);

        // --- Fast Inner Loop ---
        for (; k < k_sub_max; k++) {
            h_seg[k - k_start] = (T_COMPLEX)(a * phase); // Corrected indexing and cast
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
    double f0, double f1, double f2, double f3, 
    double a0, double a1, double a2, double a3, 
    double p0, double p1, double p2, double p3,
    double df, int64_t k_start, const int64_t kmax, // Renamed k to k_start
    const int update_interval)
{
    double f;
    double h2 = df * df;
    double h3 = h2 * df;

    // --- Denominators for Lagrange basis (constant for segment) ---
    double denom0 = (f0-f1)*(f0-f2)*(f0-f3);
    double denom1 = (f1-f0)*(f1-f2)*(f1-f3);
    double denom2 = (f2-f0)*(f2-f1)*(f2-f3);
    double denom3 = (f3-f0)*(f3-f1)*(f3-f2);

    // --- Get power-basis coefficients: c3*f^3 + c2*f^2 + c1*f + c0 ---
    // --- Amplitude ---
    double c_a3 = a0/denom0 + a1/denom1 + a2/denom2 + a3/denom3;
    double c_a2 = -a0*(f1+f2+f3)/denom0 - a1*(f0+f2+f3)/denom1 - a2*(f0+f1+f3)/denom2 - a3*(f0+f1+f2)/denom3;
    double c_a1 = a0*(f1*f2 + f1*f3 + f2*f3)/denom0 + a1*(f0*f2 + f0*f3 + f2*f3)/denom1 + a2*(f0*f1 + f0*f3 + f1*f3)/denom2 + a3*(f0*f1 + f0*f2 + f1*f2)/denom3;
    double c_a0 = -a0*(f1*f2*f3)/denom0 - a1*(f0*f2*f3)/denom1 - a2*(f0*f1*f3)/denom2 - a3*(f0*f1*f2)/denom3;
    // --- Phase ---
    double c_p3 = p0/denom0 + p1/denom1 + p2/denom2 + p3/denom3;
    double c_p2 = -p0*(f1+f2+f3)/denom0 - p1*(f0+f2+f3)/denom1 - p2*(f0+f1+f3)/denom2 - p3*(f0+f1+f2)/denom3;
    double c_p1 = p0*(f1*f2 + f1*f3 + f2*f3)/denom0 + p1*(f0*f2 + f0*f3 + f2*f3)/denom1 + p2*(f0*f1 + f0*f3 + f1*f3)/denom2 + p3*(f0*f1 + f0*f2 + f1*f2)/denom3;
    double c_p0 = -p0*(f1*f2*f3)/denom0 - p1*(f0*f2*f3)/denom1 - p2*(f0*f1*f3)/denom2 - p3*(f0*f1*f2)/denom3;

    // --- Finite difference constants (constant for segment) ---
    double d3_a_const = 6 * c_a3 * h3;
    double d3_p_const = 6 * c_p3 * h3;
    std::complex<double> d3_phase_const = std::polar(1.0, d3_p_const);

    // Stepper variables
    double a, p, d1_a, d1_p, d2_a, d2_p;
    std::complex<double> phase, d1_phase, d2_phase;
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
        phase = std::polar(1.0, p);
        d1_phase = std::polar(1.0, d1_p);
        d2_phase = std::polar(1.0, d2_p);

        // --- Fast Inner Loop ---
        for (; k < k_sub_max; k++) {
            h_seg[k - k_start] = (T_COMPLEX)(a * phase); // Corrected indexing and cast

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
 * FAST QUARTIC INTERPOLATION (HELPER)
 *
 * =====================================================================
 */

template <class T, class T_COMPLEX>
static inline void _decomp_Qcode_segment(
    T_COMPLEX* h_seg,
    double f0, double f1, double f2, double f3, double f4, // 5 freqs
    double a0, double a1, double a2, double a3, double a4, // 5 amps
    double p0, double p1, double p2, double p3, double p4, // 5 phases
    double df, int64_t k_start, const int64_t kmax,
    const int update_interval)
{
    double f;
    double h2 = df * df;
    double h3 = h2 * df;
    double h4 = h3 * df;

    // --- Denominators for Lagrange basis (constant for segment) ---
    double denom0 = (f0-f1)*(f0-f2)*(f0-f3)*(f0-f4);
    double denom1 = (f1-f0)*(f1-f2)*(f1-f3)*(f1-f4);
    double denom2 = (f2-f0)*(f2-f1)*(f2-f3)*(f2-f4);
    double denom3 = (f3-f0)*(f3-f1)*(f3-f2)*(f3-f4);
    double denom4 = (f4-f0)*(f4-f1)*(f4-f2)*(f4-f3);

    // --- Get power-basis coefficients: c4*f^4 + c3*f^3 + c2*f^2 + c1*f + c0 ---
    // --- Amplitude ---
    double c_a4 = a0/denom0 + a1/denom1 + a2/denom2 + a3/denom3 + a4/denom4;
    double c_a3 = -a0*(f1+f2+f3+f4)/denom0 - a1*(f0+f2+f3+f4)/denom1 - a2*(f0+f1+f3+f4)/denom2 - a3*(f0+f1+f2+f4)/denom3 - a4*(f0+f1+f2+f3)/denom4;
    double c_a2 = a0*(f1*f2+f1*f3+f1*f4+f2*f3+f2*f4+f3*f4)/denom0 + a1*(f0*f2+f0*f3+f0*f4+f2*f3+f2*f4+f3*f4)/denom1 + a2*(f0*f1+f0*f3+f0*f4+f1*f3+f1*f4+f3*f4)/denom2 + a3*(f0*f1+f0*f2+f0*f4+f1*f2+f1*f4+f2*f4)/denom3 + a4*(f0*f1+f0*f2+f0*f3+f1*f2+f1*f3+f2*f3)/denom4;
    double c_a1 = -a0*(f1*f2*f3+f1*f2*f4+f1*f3*f4+f2*f3*f4)/denom0 - a1*(f0*f2*f3+f0*f2*f4+f0*f3*f4+f2*f3*f4)/denom1 - a2*(f0*f1*f3+f0*f1*f4+f0*f3*f4+f1*f3*f4)/denom2 - a3*(f0*f1*f2+f0*f1*f4+f0*f2*f4+f1*f2*f4)/denom3 - a4*(f0*f1*f2+f0*f1*f3+f0*f2*f3+f1*f2*f3)/denom4;
    double c_a0 = a0*(f1*f2*f3*f4)/denom0 + a1*(f0*f2*f3*f4)/denom1 + a2*(f0*f1*f3*f4)/denom2 + a3*(f0*f1*f2*f4)/denom3 + a4*(f0*f1*f2*f3)/denom4;
    // --- Phase ---
    double c_p4 = p0/denom0 + p1/denom1 + p2/denom2 + p3/denom3 + p4/denom4;
    double c_p3 = -p0*(f1+f2+f3+f4)/denom0 - p1*(f0+f2+f3+f4)/denom1 - p2*(f0+f1+f3+f4)/denom2 - p3*(f0+f1+f2+f4)/denom3 - p4*(f0+f1+f2+f3)/denom4;
    double c_p2 = p0*(f1*f2+f1*f3+f1*f4+f2*f3+f2*f4+f3*f4)/denom0 + p1*(f0*f2+f0*f3+f0*f4+f2*f3+f2*f4+f3*f4)/denom1 + p2*(f0*f1+f0*f3+f0*f4+f1*f3+f1*f4+f3*f4)/denom2 + p3*(f0*f1+f0*f2+f0*f4+f1*f2+f1*f4+f2*f4)/denom3 + p4*(f0*f1+f0*f2+f0*f3+f1*f2+f1*f3+f2*f3)/denom4;
    double c_p1 = -p0*(f1*f2*f3+f1*f2*f4+f1*f3*f4+f2*f3*f4)/denom0 - p1*(f0*f2*f3+f0*f2*f4+f0*f3*f4+f2*f3*f4)/denom1 - p2*(f0*f1*f3+f0*f1*f4+f0*f3*f4+f1*f3*f4)/denom2 - p3*(f0*f1*f2+f0*f1*f4+f0*f2*f4+f1*f2*f4)/denom3 - p4*(f0*f1*f2+f0*f1*f3+f0*f2*f3+f1*f2*f3)/denom4;
    double c_p0 = p0*(f1*f2*f3*f4)/denom0 + p1*(f0*f2*f3*f4)/denom1 + p2*(f0*f1*f3*f4)/denom2 + p3*(f0*f1*f2*f4)/denom3 + p4*(f0*f1*f2*f3)/denom4;

    // --- Finite difference constants (constant for segment) ---
    double d4_a_const = 24 * c_a4 * h4;
    double d4_p_const = 24 * c_p4 * h4;
    std::complex<double> d4_phase_const = std::polar(1.0, d4_p_const);

    // Stepper variables
    double a, p, d1_a, d1_p, d2_a, d2_p, d3_a, d3_p;
    std::complex<double> phase, d1_phase, d2_phase, d3_phase;
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
        phase = std::polar(1.0, p);
        d1_phase = std::polar(1.0, d1_p);
        d2_phase = std::polar(1.0, d2_p);
        d3_phase = std::polar(1.0, d3_p);

        // --- Fast Inner Loop ---
        for (; k < k_sub_max; k++) {
            h_seg[k - k_start] = (T_COMPLEX)(a * phase); // Corrected indexing and cast

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
 * MASTER TEMPLATED INTERPOLATOR (MAIN) - NEWLY ADDED
 *
 * =====================================================================
 */
template <class T, class T_COMPLEX>
static inline void _decomp_main_loop(
    int degree, // The degree of interpolation to aim for
    T_COMPLEX * h,
    double delta_f, // Force double here
    const int64_t hlen,
    const int64_t start_index,
    T * sample_frequencies,
    T * amp,
    T * phase,
    const int64_t sflen,
    const int64_t imin)
{
    int64_t k, kmax;
    int64_t last_findex = start_index;
    const int update_interval = 128;
    // Use double for all local scalar logic
    double f0, f1, f2, f3, f4, a0, a1, a2, a3, a4, p0, p1, p2, p3, p4;

    // Determine the maximum possible degree based on number of points
    int max_possible_degree;
    if (sflen < 3) max_possible_degree = 1; // Not enough for quadratic
    else if (sflen < 4) max_possible_degree = 2; // Not enough for cubic
    else if (sflen < 5) max_possible_degree = 3; // Not enough for quartic
    else max_possible_degree = 4; // We only implement up to quartic

    // Use the lower of the requested degree or the max possible
    if (degree > max_possible_degree) {
        degree = max_possible_degree;
    }

    // zero out the beginning
    memset(h, 0, sizeof(T_COMPLEX)*start_index);

    for (int64_t i = imin; i < sflen-1; i++) {
        // Get segment boundaries (always needed) - cast to double
        f1 = (double)sample_frequencies[i];
        f2 = (double)sample_frequencies[i+1];
        a1 = (double)amp[i];
        a2 = (double)amp[i+1];
        p1 = (double)phase[i];
        p2 = (double)phase[i+1];

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

        // --- Gracefully degrade at the boundaries ---
        // Start with the highest requested degree
        int current_degree = degree; 

        // Degrade based on start position
        if (i == imin) current_degree = 1; // First segment must be linear
        else if (i == imin + 1 && current_degree > 2) current_degree = 2;
        else if (i == imin + 2 && current_degree > 3) current_degree = 3;

        // Degrade based on end position (Quartic needs i+3, Cubic needs i+2)
        if (current_degree > 3 && i >= sflen - 3) current_degree = 3;
        if (current_degree > 2 && i >= sflen - 2) current_degree = 2;

        // --- Call the correct helper based on current_degree ---
        // Note: Helpers now accept explicit doubles for scalars
        switch (current_degree) {
            case 1:
                _decomp_ccode_segment<T, T_COMPLEX>(
                    &h[k], f1, f2, a1, a2, p1, p2,
                    delta_f, k, kmax, update_interval);
                break;

            case 2:
                f0 = (double)sample_frequencies[i-1];
                a0 = (double)amp[i-1]; p0 = (double)phase[i-1];
                _decomp_qcode_segment<T, T_COMPLEX>(
                    &h[k], f0, f1, f2, a0, a1, a2, p0, p1, p2,
                    delta_f, k, kmax, update_interval);
                break;

            case 3:
                f0 = (double)sample_frequencies[i-1];
                f3 = (double)sample_frequencies[i+2];
                a0 = (double)amp[i-1]; a3 = (double)amp[i+2];
                p0 = (double)phase[i-1]; p3 = (double)phase[i+2];
                _decomp_tcode_segment<T, T_COMPLEX>(
                    &h[k], f0, f1, f2, f3, a0, a1, a2, a3, p0, p1, p2, p3,
                    delta_f, k, kmax, update_interval);
                break;
            
            case 4:
            default: // Catches degree >= 4
                f0 = (double)sample_frequencies[i-1];
                f3 = (double)sample_frequencies[i+2];
                f4 = (double)sample_frequencies[i+3];
                a0 = (double)amp[i-1]; a3 = (double)amp[i+2]; a4 = (double)amp[i+3];
                p0 = (double)phase[i-1]; p3 = (double)phase[i+2]; p4 = (double)phase[i+3];
                _decomp_Qcode_segment<T, T_COMPLEX>(
                    &h[k], f0, f1, f2, f3, f4, a0, a1, a2, a3, a4, p0, p1, p2, p3, p4,
                    delta_f, k, kmax, update_interval);
                break;
        }
        last_findex = kmax;
    }

    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(T_COMPLEX)*(hlen-last_findex));
    }
}


/* =====================================================================
 *
 * PUBLIC-FACING STUB FUNCTIONS (REDUCED BOILERPLATE)
 *
 * =====================================================================
 */

// --- LINEAR ---
void _decomp_ccode_double(std::complex<double> * h, double delta_f, const int64_t hlen, const int64_t start_index, double * sample_frequencies, double * amp, double * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<double, std::complex<double> >(1, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}
// Updated delta_f to double
void _decomp_ccode_float(std::complex<float> * h, double delta_f, const int64_t hlen, const int64_t start_index, float * sample_frequencies, float * amp, float * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<float, std::complex<float> >(1, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}

// --- QUADRATIC ---
void _decomp_qcode_double(std::complex<double> * h, double delta_f, const int64_t hlen, const int64_t start_index, double * sample_frequencies, double * amp, double * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<double, std::complex<double> >(2, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}
// Updated delta_f to double
void _decomp_qcode_float(std::complex<float> * h, double delta_f, const int64_t hlen, const int64_t start_index, float * sample_frequencies, float * amp, float * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<float, std::complex<float> >(2, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}

// --- CUBIC ---
void _decomp_tcode_double(std::complex<double> * h, double delta_f, const int64_t hlen, const int64_t start_index, double * sample_frequencies, double * amp, double * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<double, std::complex<double> >(3, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}
// Updated delta_f to double
void _decomp_tcode_float(std::complex<float> * h, double delta_f, const int64_t hlen, const int64_t start_index, float * sample_frequencies, float * amp, float * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<float, std::complex<float> >(3, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}

// --- QUARTIC ---
void _decomp_Qcode_double(std::complex<double> * h, double delta_f, const int64_t hlen, const int64_t start_index, double * sample_frequencies, double * amp, double * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<double, std::complex<double> >(4, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}
// Updated delta_f to double
void _decomp_Qcode_float(std::complex<float> * h, double delta_f, const int64_t hlen, const int64_t start_index, float * sample_frequencies, float * amp, float * phase, const int64_t sflen, const int64_t imin) {
    _decomp_main_loop<float, std::complex<float> >(4, h, delta_f, hlen, start_index, sample_frequencies, amp, phase, sflen, imin);
}
