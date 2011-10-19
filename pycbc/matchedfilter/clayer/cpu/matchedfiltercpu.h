// Copyright (C) 2011 Karsten Wiesner
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
// matchedfiltercpu declarations 


#ifndef MATCHEDFILTERCPU_H
#define MATCHEDFILTERCPU_H

#include <stdlib.h>

typedef struct
{
    int nothing_to_define_yet;
}
matched_filter_cpu_t;

void gen_snr_cpu(cpu_context_t* context,
                 real_vector_single_cpu_t* snr,
                 complex_vector_single_cpu_t* stilde, 
                 complex_vector_single_cpu_t* htilde,
                 complex_vector_single_cpu_t* q,
                 complex_vector_single_cpu_t* qtilde,
                 /*complex_fft_plan_t* plan,*/
                 double f_min,
                 double sigma_sq);

#endif /* MATCHEDFILTERCPU_H */
