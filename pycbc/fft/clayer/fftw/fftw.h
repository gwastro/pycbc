// Copyright (C) 2011 Josh Willis
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
// fftw plan constructor and destructor prototypes for pycbc

#ifndef PYCBC_FFTW_H
#define PYCBC_FFTW_H

#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "datavectorcpu.h"

/* Struct common for single-precision plans */
#define FFTWF_PLAN_STRUCT_DATA \
  fftwf_plan    theplan; \
  unsigned long size; \
  int           fwdflag;

/* Struct common for double-precision plans */
#define FFTW_PLAN_STRUCT_DATA \
  fftw_plan     theplan; \
  unsigned long size; \
  int           fwdflag;

/* Now the actual structs our functions use */

typedef struct {
FFTWF_PLAN_STRUCT_DATA
} fftw_real_single_plan;

typedef struct {
FFTW_PLAN_STRUCT_DATA
} fftw_real_double_plan;

typedef struct {
FFTWF_PLAN_STRUCT_DATA
} fftw_complex_single_plan;

typedef struct {
FFTW_PLAN_STRUCT_DATA
} fftw_complex_double_plan;


/*
  Functions to execute FFTW plans
*/

void execute_complex_single_fft(complex_vector_single_cpu_t *output,
				complex_vector_single_cpu_t *input,
				fftw_complex_single_plan *plan);
void execute_real_single_forward_fft(complex_vector_single_cpu_t *output,
				     real_vector_single_cpu_t *input,
				     fftw_real_single_plan *plan);
void execute_real_single_reverse_fft(real_vector_single_cpu_t *output,
				     complex_vector_single_cpu_t *input,
				     fftw_real_single_plan *plan);

void execute_complex_double_fft(complex_vector_double_cpu_t *output,
				complex_vector_double_cpu_t *input,
				fftw_complex_double_plan *plan);
void execute_real_double_forward_fft(complex_vector_double_cpu_t *output,
				     real_vector_double_cpu_t *input,
				     fftw_real_double_plan *plan);
void execute_real_double_reverse_fft(real_vector_double_cpu_t *output,
				     complex_vector_double_cpu_t *input,
				     fftw_real_double_plan *plan);

#endif /* PYCBC_FFTW_H */
