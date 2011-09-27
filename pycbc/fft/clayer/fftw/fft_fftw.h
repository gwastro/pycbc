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

#ifndef FFT_FFTW_H
#define FFT_FFTW_H

#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "fft_fftw_private.h"
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"

/*
 Prototypes of constructors/destructors.  These are just SWIG wrapped into 
 functions; one layer above the interface file calls the appropriate function based 
 on where the memory will live.
*/

fft_real_single_plan_fftw *new_fft_real_single_plan_fftw(unsigned long size, 
							 int fwdflag, 
							 int measurelvl);
void delete_fft_real_single_plan_fftw(fft_real_single_plan_fftw *self);

fft_real_double_plan_fftw *new_fft_real_double_plan_fftw(unsigned long size, 
							 int fwdflag,
							 int measurelvl);
void delete_fft_real_double_plan_fftw(fft_real_double_plan_fftw *self);

fft_complex_single_plan_fftw *new_fft_complex_single_plan_fftw(unsigned long size,
							       int fwdflag, 
							       int measurelvl);
void delete_fft_complex_single_plan_fftw(fft_complex_single_plan_fftw *self);

fft_complex_double_plan_fftw *new_fft_complex_double_plan_fftw(unsigned long size, 
							       int fwdflag,
							       int measurelvl);
void delete_fft_complex_double_plan_fftw(fft_complex_double_plan_fftw *self);

/*
  Functions to execute FFTW plans
*/

void execute_complex_single_fft_fftw(complex_vector_single_cpu_t *output,
				     complex_vector_single_cpu_t *input,
				     fft_complex_single_plan_fftw *plan);
void execute_real_single_forward_fft_fftw(complex_vector_single_cpu_t *output,
					  real_vector_single_cpu_t *input,
					  fft_real_single_plan_fftw *plan);
void execute_real_single_reverse_fft_fftw(real_vector_single_cpu_t *output,
					  complex_vector_single_cpu_t *input,
					  fft_real_single_plan_fftw *plan);

void execute_complex_double_fft_fftw(complex_vector_double_cpu_t *output,
				     complex_vector_double_cpu_t *input,
				     fft_complex_double_plan_fftw *plan);
void execute_real_double_forward_fft_fftw(complex_vector_double_cpu_t *output,
					  real_vector_double_cpu_t *input,
					  fft_real_double_plan_fftw *plan);
void execute_real_double_reverse_fft_fftw(real_vector_double_cpu_t *output,
					  complex_vector_double_cpu_t *input,
					  fft_real_double_plan_fftw *plan);

#endif /* FFT_FFTW_H */
