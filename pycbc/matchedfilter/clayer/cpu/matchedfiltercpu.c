// Copyright (C) 2011 Badri Krishnan, Karsten Wiesner
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


//
// =============================================================================
//
//                                   Preamble
//
// =============================================================================
//
// matchedfiltercpu  

#include <stdio.h>
#include "datavectorcpu.h"
#include "pycbccpu.h"
#include "fftw.h" 
#include "matchedfiltercpu_private.h"
#include "matchedfiltercpu.h"

matched_filter_cpu_t* new_matched_filter_cpu_t(cpu_context_t* context)
{

    matched_filter_cpu_t* c;
    c = (matched_filter_cpu_t*) malloc( sizeof(matched_filter_cpu_t) );

    return c;
}


void delete_matched_filter_cpu_t( matched_filter_cpu_t* p )
{
    free( p );
}


void gen_snr_cpu( cpu_context_t* context,
		          real_vector_single_cpu_t* snr,
		          complex_vector_single_cpu_t* stilde, 
		          complex_vector_single_cpu_t* htilde,
		          complex_vector_single_cpu_t* q,
		          complex_vector_single_cpu_t* qtilde,
		          fftw_complex_single_plan* plan,
		          double f_min,
		          double sigma_sq)
{    

    static unsigned cnt=0;
    int j;
    const int N = 2 * (stilde->meta_data.vector_length - 1);
    double norm = 4.0 / ((double) N * (double) N * stilde->meta_data.delta_x); 
    norm *= norm;
    
    printf("%d: called gen_snr_cpu with context: %p, snr: %p s: %p h: %p q: %p qtilde %p fft_plan: %p, fmin: %f, sigma-sq: %f\n", 
           cnt++, context, snr, stilde, htilde, q, qtilde, plan, f_min, sigma_sq);
    
    /* perform the correlation */
    correlate_complex_freq_vectors( qtilde, stilde, htilde, f_min );

    /* execute complex single ifft (output, input, plan) */
    execute_complex_single_fft(q, qtilde, plan);

    /* normalize the snr */
    for ( j = 0; j < q->meta_data.vector_length; ++j )
      snr->data[j] = (norm/sigma_sq) * 
        (__real__ q->data[j] * __real__ q->data[j] + 
         __imag__ q->data[j] * __imag__ q->data[j]);
  
    snr->meta_data.delta_x = q->meta_data.delta_x;
  
    return;
}

void find_max_cpu( cpu_context_t* context,
                   double* max,
                   unsigned long* index,
                   real_vector_single_cpu_t* snr)
{    
    unsigned long i;

    *max = 0;
    for (i=0; i < snr->meta_data.vector_length; i++) {
        if (snr->data[i] > *max) {
            *max = snr->data[i];
            *index = i;
        }	 
    }
    
    return;
}


/* freq domain correlation of two vectors (suitably whitened) */
void correlate_complex_freq_vectors( complex_vector_single_cpu_t* out,
				     complex_vector_single_cpu_t* x, 
				     complex_vector_single_cpu_t* y, 
				     double f_min)
{
    int k;
    const double df = x->meta_data.delta_x;
    const int bin_min = f_min / df > 1 ? f_min / df : 1;
    const int len = x->meta_data.vector_length;

    for ( k = bin_min; k < len-1; k++ )
    {
        __real__ out->data[k] = __real__ x->data[k] * __real__ y->data[k] - 
                                __imag__ x->data[k] * 
                                ( 0.0 - __imag__ y->data[k] );
        __imag__ out->data[k] = __real__ x->data[k] * 
                                (0.0 - __imag__ y->data[k]) + 
                                __imag__ x->data[k] * __real__ y->data[k];
    }

    out->meta_data.delta_x = df;

    return;
}
