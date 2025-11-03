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
// This code expects to be passed:
// h: array of complex doubles
//      the output array to write the results to.
// delta_f: double
//      the df of the output array
// hlen: int
//      the length of h
// start_index: int
//      the index to start the waveform in the output
//      frequency series; i.e., ceil(f_lower/df)
// sample_frequencies: array of real doubles
//      the frequencies at which the compressed waveform is sampled
// amp: array of real doubles
//      the amplitude of the waveform at the sample frequencies
// phase: array of real doubles
//      the phase of the waveform at the sample frequencies
// sflen: int
//      the length of the sample frequencies
// imin: int
//      the index to start at in the compressed series

// We will cast the output to a double array for faster processing.
// This takes advantage of the fact that complex arrays store
// their real and imaginary values next to each other in memory.

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
    double* outptr = (double*) h;

    // for keeping track of where in the output frequencies we are
    int findex, next_sfindex, kmax;

    // variables for computing the interpolation
    double df = (double) delta_f;
    double inv_df = 1./df;
    double f, inv_sdf;
    double sf, this_amp, this_phi;
    double next_sf = sample_frequencies[imin];
    double next_amp = amp[imin];
    double next_phi = phase[imin];
    double m_amp, b_amp;
    double m_phi, b_phi;
    double interp_amp, interp_phi;

    // variables for updating each interpolated frequency
    double h_re, h_im, incrh_re, incrh_im;
    double g_re, g_im, incrg_re, incrg_im;
    double dphi_re, dphi_im;

    // we will re-compute cos/sin of the phase at the following intervals:
    int update_interval = 128;

    // zero out the beginning
    memset(outptr, 0, sizeof(*outptr)*2*start_index);

    // move to the start position
    outptr += 2*start_index;
    findex = start_index;

    // cycle over the compressed samples
    for (int ii=imin; ii<(sflen-1); ii++){
        // update the linear interpolations
        sf = next_sf;
        next_sf = (double) sample_frequencies[ii+1];
        next_sfindex = (int) ceil(next_sf * inv_df);
        if (next_sfindex > hlen)
            next_sfindex = hlen;
        inv_sdf = 1./(next_sf - sf);
        this_amp = next_amp;
        next_amp = (double) amp[ii+1];
        this_phi = next_phi;
        next_phi = (double) phase[ii+1];
        m_amp = (next_amp - this_amp)*inv_sdf;
        b_amp = this_amp - m_amp*sf;
        m_phi = (next_phi - this_phi)*inv_sdf;
        b_phi = this_phi - m_phi*sf;

        // cycle over the interpolated points between this and the next
        // compressed sample
        while (findex < next_sfindex){
            // for the first step, compute the value of h from the interpolated
            // amplitude and phase
            f = findex*df;
            interp_amp = m_amp * f + b_amp;
            interp_phi = m_phi * f + b_phi;
            dphi_re = cos(m_phi * df);
            dphi_im = sin(m_phi * df);
            h_re = interp_amp * cos(interp_phi);
            h_im = interp_amp * sin(interp_phi);
            g_re = m_amp * df * cos(interp_phi);
            g_im = m_amp * df * sin(interp_phi);

            // save and update counters
            *outptr = h_re;
            *(outptr+1) = h_im;
            outptr += 2;
            findex++;

            // for the next update_interval steps, compute h by incrementing
            // the last h
            kmax = findex + update_interval;
            if (kmax > next_sfindex)
                kmax = next_sfindex;
            while (findex < kmax){
                incrh_re = h_re * dphi_re - h_im * dphi_im;
                incrh_im = h_re * dphi_im + h_im * dphi_re;
                incrg_re = g_re * dphi_re - g_im * dphi_im;
                incrg_im = g_re * dphi_im + g_im * dphi_re;
                h_re = incrh_re + incrg_re;
                h_im = incrh_im + incrg_im;
                g_re = incrg_re;
                g_im = incrg_im;

                // save and update counters
                *outptr = h_re;
                *(outptr+1) = h_im;
                outptr += 2;
                findex++;
            }
        }
        if (next_sfindex == hlen){
            break;
        }
    }

    // zero out the rest of the array
    memset(outptr, 0, sizeof(*outptr)*2*(hlen-findex));
}

/* I'm just copying this wholesale for float precision. If someone wants to
 * improve this function in the future, consider using C++'s dynamic types
 * where you can define the float and double functions together.
 */

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
    float* outptr = (float*) h;

    // for keeping track of where in the output frequencies we are
    int findex, next_sfindex, kmax;

    // variables for computing the interpolation
    float df = (float) delta_f;
    float inv_df = 1./df;
    float f, inv_sdf;
    float sf, this_amp, this_phi;
    float next_sf = sample_frequencies[imin];
    float next_amp = amp[imin];
    float next_phi = phase[imin];
    float m_amp, b_amp;
    float m_phi, b_phi;
    float interp_amp, interp_phi;

    // variables for updating each interpolated frequency
    float h_re, h_im, incrh_re, incrh_im;
    float g_re, g_im, incrg_re, incrg_im;
    float dphi_re, dphi_im;

    // we will re-compute cos/sin of the phase at the following intervals:
    int update_interval = 128;

    // zero out the beginning
    memset(outptr, 0, sizeof(*outptr)*2*start_index);

    // move to the start position
    outptr += 2*start_index;
    findex = start_index;

    // cycle over the compressed samples
    for (int ii=imin; ii<(sflen-1); ii++){
        // update the linear interpolations
        sf = next_sf;
        next_sf = (float) sample_frequencies[ii+1];
        next_sfindex = (int) ceil(next_sf * inv_df);
        if (next_sfindex > hlen)
            next_sfindex = hlen;
        inv_sdf = 1./(next_sf - sf);
        this_amp = next_amp;
        next_amp = (float) amp[ii+1];
        this_phi = next_phi;
        next_phi = (float) phase[ii+1];
        m_amp = (next_amp - this_amp)*inv_sdf;
        b_amp = this_amp - m_amp*sf;
        m_phi = (next_phi - this_phi)*inv_sdf;
        b_phi = this_phi - m_phi*sf;

        // cycle over the interpolated points between this and the next
        // compressed sample
        while (findex < next_sfindex){
            // for the first step, compute the value of h from the interpolated
            // amplitude and phase
            f = findex*df;
            interp_amp = m_amp * f + b_amp;
            interp_phi = m_phi * f + b_phi;
            dphi_re = cos(m_phi * df);
            dphi_im = sin(m_phi * df);
            h_re = interp_amp * cos(interp_phi);
            h_im = interp_amp * sin(interp_phi);
            g_re = m_amp * df * cos(interp_phi);
            g_im = m_amp * df * sin(interp_phi);

            // save and update counters
            *outptr = h_re;
            *(outptr+1) = h_im;
            outptr += 2;
            findex++;

            // for the next update_interval steps, compute h by incrementing
            // the last h
            kmax = findex + update_interval;
            if (kmax > next_sfindex)
                kmax = next_sfindex;
            while (findex < kmax){
                incrh_re = h_re * dphi_re - h_im * dphi_im;
                incrh_im = h_re * dphi_im + h_im * dphi_re;
                incrg_re = g_re * dphi_re - g_im * dphi_im;
                incrg_im = g_re * dphi_im + g_im * dphi_re;
                h_re = incrh_re + incrg_re;
                h_im = incrh_im + incrg_im;
                g_re = incrg_re;
                g_im = incrg_im;

                // save and update counters
                *outptr = h_re;
                *(outptr+1) = h_im;
                outptr += 2;
                findex++;
            }
        }
        if (next_sfindex == hlen){
            break;
        }
    }

    // zero out the rest of the array
    memset(outptr, 0, sizeof(*outptr)*2*(hlen-findex));
}


/* =====================================================================
 *
 * QUADRATIC INTERPOLATION FUNCTIONS ADDED BELOW
 *
 * =====================================================================
 */


/**
* Decompress a waveform using highly optimized quadratic interpolation
*
* We quadratically interpolate the amplitude and phase to fill the values
* between the sample frequencies.
*
* This uses a finite-difference stepper for both amplitude and phase,
* eliminating all transcendental function calls from the innermost loop.
*
* BUGFIX: The first segment (i == imin) uses the fast ITERATIVE LINEAR
* interpolation to avoid boundary conditions.
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
    int64_t k, kmax, k_sub_max;
    int64_t last_findex = start_index;
    int update_interval = 128; // For correcting FP error
    double f, f0, f1, f2, p0, p1, p2, a0, a1, a2;
    double df = delta_f;
    double df_inv = 1.0 / df;
    double h2 = df * df;

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

        if (kmax > hlen) {
            kmax = hlen;
        }

        // === BUGFIX: Use robust, FAST linear interpolation for the first segment ===
        if (i == imin) {
            double* outptr = (double*) &h[k];
            double inv_sdf = 1./(f2 - f1);
            double m_amp = (a2 - a1)*inv_sdf;
            double b_amp = a1 - m_amp*f1;
            double m_phi = (p2 - p1)*inv_sdf;
            double b_phi = p1 - m_phi*f1;

            double h_re, h_im, g_re, g_im, dphi_re, dphi_im;
            double incrh_re, incrh_im, incrg_re, incrg_im;

            int64_t findex = k;
            while (findex < kmax){
                f = findex*df;
                double interp_amp = m_amp * f + b_amp;
                double interp_phi = m_phi * f + b_phi;
                dphi_re = cos(m_phi * df);
                dphi_im = sin(m_phi * df);
                h_re = interp_amp * cos(interp_phi);
                h_im = interp_amp * sin(interp_phi);
                g_re = m_amp * df * cos(interp_phi);
                g_im = m_amp * df * sin(interp_phi);

                *outptr = h_re;
                *(outptr+1) = h_im;
                outptr += 2;
                findex++;

                k_sub_max = findex + update_interval;
                if (k_sub_max > kmax) k_sub_max = kmax;
                
                while (findex < k_sub_max){
                    incrh_re = h_re * dphi_re - h_im * dphi_im;
                    incrh_im = h_re * dphi_im + h_im * dphi_re;
                    incrg_re = g_re * dphi_re - g_im * dphi_im;
                    incrg_im = g_re * dphi_im + g_im * dphi_re;
                    h_re = incrh_re + incrg_re;
                    h_im = incrh_im + incrg_im;
                    g_re = incrg_re;
                    g_im = incrg_im;

                    *outptr = h_re;
                    *(outptr+1) = h_im;
                    outptr += 2;
                    findex++;
                }
            }
        }
        // === Use FAST quadratic interpolation for all subsequent segments ===
        else {
            // Get 3rd point for quadratic
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];

            // Denominators for Lagrange basis (constant for segment)
            double denom0 = (f0 - f1) * (f0 - f2);
            double denom1 = (f1 - f0) * (f1 - f2);
            double denom2 = (f2 - f0) * (f2 - f1);

            // --- Get power-basis coefficients: c2*f^2 + c1*f + c0 ---
            // c2 = v0/d0 + v1/d1 + v2/d2
            double c_a2 = a0/denom0 + a1/denom1 + a2/denom2;
            double c_p2 = p0/denom0 + p1/denom1 + p2/denom2;
            // c1 = -v0(f1+f2)/d0 - v1(f0+f2)/d1 - v2(f0+f1)/d2
            double c_a1 = -a0*(f1+f2)/denom0 - a1*(f0+f2)/denom1 - a2*(f0+f1)/denom2;
            double c_p1 = -p0*(f1+f2)/denom0 - p1*(f0+f2)/denom1 - p2*(f0+f1)/denom2;
            // c0 = v0*f1*f2/d0 + v1*f0*f2/d1 + v2*f0*f1/d2
            double c_a0 = a0*f1*f2/denom0 + a1*f0*f2/denom1 + a2*f0*f1/denom2;
            double c_p0 = p0*f1*f2/denom0 + p1*f0*f2/denom1 + p2*f0*f1/denom2;
            
            // --- Finite difference constants (constant for segment) ---
            double d2_a = 2 * c_a2 * h2;
            double d2_p = 2 * c_p2 * h2;
            std::complex<double> d2_phase = std::polar(1.0, d2_p);

            // Stepper variables
            double a, p, d1_a, d1_p;
            std::complex<double> phase, d_phase;

            for (; k < kmax; ) {
                // --- Re-calculate steppers to correct FP error ---
                f = k * df;
                k_sub_max = k + update_interval;
                if (k_sub_max > kmax) k_sub_max = kmax;

                // Initial values for this sub-block
                a = c_a2*f*f + c_a1*f + c_a0;
                p = c_p2*f*f + c_p1*f + c_p0;
                
                // First differences at f
                d1_a = c_a2*(2*f*df + h2) + c_a1*df;
                d1_p = c_p2*(2*f*df + h2) + c_p1*df;
                
                // Complex phase steppers
                phase = std::polar(1.0, p);
                d_phase = std::polar(1.0, d1_p);

                // --- Fast Inner Loop ---
                // (No sin/cos/polar calls)
                for (; k < k_sub_max; k++) {
                    h[k] = a * phase;
                    
                    // Step phase forward
                    phase = phase * d_phase;
                    d_phase = d_phase * d2_phase;
                    
                    // Step amplitude forward
                    a = a + d1_a;
                    d1_a = d1_a + d2_a;
                }
            }
        }
        last_findex = kmax;
    }

    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<double>)*(hlen-last_findex));
    }
}


/**
* Decompress a waveform using highly optimized quadratic interpolation (float version)
*/
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
    int64_t k, kmax, k_sub_max;
    int64_t last_findex = start_index;
    int update_interval = 128; // For correcting FP error
    float f, f0, f1, f2, p0, p1, p2, a0, a1, a2;
    float df = delta_f;
    float df_inv = 1.0f / df;
    float h2 = df * df;

    // zero out the beginning
    memset(h, 0, sizeof(std::complex<float>)*start_index);

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

        if (kmax > hlen) {
            kmax = hlen;
        }

        // === BUGFIX: Use robust, FAST linear interpolation for the first segment ===
        if (i == imin) {
            float* outptr = (float*) &h[k];
            float inv_sdf = 1.0f/(f2 - f1);
            float m_amp = (a2 - a1)*inv_sdf;
            float b_amp = a1 - m_amp*f1;
            float m_phi = (p2 - p1)*inv_sdf;
            float b_phi = p1 - m_phi*f1;

            float h_re, h_im, g_re, g_im, dphi_re, dphi_im;
            float incrh_re, incrh_im, incrg_re, incrg_im;

            int64_t findex = k;
            while (findex < kmax){
                f = findex*df;
                float interp_amp = m_amp * f + b_amp;
                float interp_phi = m_phi * f + b_phi;
                dphi_re = cosf(m_phi * df);
                dphi_im = sinf(m_phi * df);
                h_re = interp_amp * cosf(interp_phi);
                h_im = interp_amp * sinf(interp_phi);
                g_re = m_amp * df * cosf(interp_phi);
                g_im = m_amp * df * sinf(interp_phi);

                *outptr = h_re;
                *(outptr+1) = h_im;
                outptr += 2;
                findex++;

                k_sub_max = findex + update_interval;
                if (k_sub_max > kmax) k_sub_max = kmax;
                
                while (findex < k_sub_max){
                    incrh_re = h_re * dphi_re - h_im * dphi_im;
                    incrh_im = h_re * dphi_im + h_im * dphi_re;
                    incrg_re = g_re * dphi_re - g_im * dphi_im;
                    incrg_im = g_re * dphi_im + g_im * dphi_re;
                    h_re = incrh_re + incrg_re;
                    h_im = incrh_im + incrg_im;
                    g_re = incrg_re;
                    g_im = incrg_im;

                    *outptr = h_re;
                    *(outptr+1) = h_im;
                    outptr += 2;
                    findex++;
                }
            }
        }
        // === Use FAST quadratic interpolation for all subsequent segments ===
        else {
            // Get 3rd point for quadratic
            f0 = sample_frequencies[i-1];
            a0 = amp[i-1];
            p0 = phase[i-1];

            // Denominators for Lagrange basis (constant for segment)
            float denom0 = (f0 - f1) * (f0 - f2);
            float denom1 = (f1 - f0) * (f1 - f2);
            float denom2 = (f2 - f0) * (f2 - f1);

            // --- Get power-basis coefficients: c2*f^2 + c1*f + c0 ---
            // c2 = v0/d0 + v1/d1 + v2/d2
            float c_a2 = a0/denom0 + a1/denom1 + a2/denom2;
            float c_p2 = p0/denom0 + p1/denom1 + p2/denom2;
            // c1 = -v0(f1+f2)/d0 - v1(f0+f2)/d1 - v2(f0+f1)/d2
            float c_a1 = -a0*(f1+f2)/denom0 - a1*(f0+f2)/denom1 - a2*(f0+f1)/denom2;
            float c_p1 = -p0*(f1+f2)/denom0 - p1*(f0+f2)/denom1 - p2*(f0+f1)/denom2;
            // c0 = v0*f1*f2/d0 + v1*f0*f2/d1 + v2*f0*f1/d2
            float c_a0 = a0*f1*f2/denom0 + a1*f0*f2/denom1 + a2*f0*f1/denom2;
            float c_p0 = p0*f1*f2/denom0 + p1*f0*f2/denom1 + p2*f0*f1/denom2;
            
            // --- Finite difference constants (constant for segment) ---
            float d2_a = 2 * c_a2 * h2;
            float d2_p = 2 * c_p2 * h2;
            std::complex<float> d2_phase = std::polar(1.0f, d2_p);

            // Stepper variables
            float a, p, d1_a, d1_p;
            std::complex<float> phase, d_phase;

            for (; k < kmax; ) {
                // --- Re-calculate steppers to correct FP error ---
                f = k * df;
                k_sub_max = k + update_interval;
                if (k_sub_max > kmax) k_sub_max = kmax;

                // Initial values for this sub-block
                a = c_a2*f*f + c_a1*f + c_a0;
                p = c_p2*f*f + c_p1*f + c_p0;
                
                // First differences at f
                d1_a = c_a2*(2*f*df + h2) + c_a1*df;
                d1_p = c_p2*(2*f*df + h2) + c_p1*df;
                
                // Complex phase steppers
                phase = std::polar(1.0f, p);
                d_phase = std::polar(1.0f, d1_p);

                // --- Fast Inner Loop ---
                // (No sin/cos/polar calls)
                for (; k < k_sub_max; k++) {
                    h[k] = a * phase;
                    
                    // Step phase forward
                    phase = phase * d_phase;
                    d_phase = d_phase * d2_phase;
                    
                    // Step amplitude forward
                    a = a + d1_a;
                    d1_a = d1_a + d2_a;
                }
            }
        }
        last_findex = kmax;
    }

    // zero out the rest of the array
    if (last_findex < hlen) {
        memset(&h[last_findex], 0, sizeof(std::complex<float>)*(hlen-last_findex));
    }
}
