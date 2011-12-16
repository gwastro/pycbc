// Copyright (C) 2011 Karsten Wiesner, Drew Keppel
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
// waveformgeneratorcpu declarations

#ifndef WAVEFORMGENERATORCPU_H
#define WAVEFORMGENERATORCPU_H

#include <stdlib.h>

typedef struct
{
    int nothing_to_define_yet;
}
waveform_generator_cpu_t;

int gen_precon_vector_TaylorF2(
    real_vector_single_cpu_t*
    );

void gen_waveform_filter_TaylorF2(
    complex_vector_single_cpu_t* waveform_filter,
    double M,
    double eta,
    int order,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three
    );

void gen_waveform_strain_TaylorF2(
    complex_vector_single_cpu_t* waveform_strain,
    double M,
    double eta,
    int order,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three,
    real_vector_single_cpu_t* precon_vec
    );
#endif /* WAVEFORMGENERATORCPU_H */
