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

