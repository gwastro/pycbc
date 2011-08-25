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

#ifndef FFTCPU_H
#define FFTCPU_H

#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "../fft.h"
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"

/*
 Prototypes of constructors/destructors.  These are just SWIG wrapped into functions; one layer above
 the interface file calls the appropriate function based on where the memory will live.
*/

fft_real_single_plan *new_fft_real_single_plan_cpu(unsigned long size, int fwdflag, int measurelvl);
void delete_fft_real_single_plan_cpu(fft_real_single_plan *self);

fft_real_double_plan *new_fft_real_double_plan_cpu(unsigned long size, int fwdflag, int measurelvl);
void delete_fft_real_double_plan_cpu(fft_real_double_plan *self);

fft_complex_single_plan *new_fft_complex_single_plan_cpu(unsigned long size, int fwdflag, int measurelvl);
void delete_fft_complex_single_plan_cpu(fft_complex_single_plan *self);

fft_complex_double_plan *new_fft_complex_double_plan_cpu(unsigned long size, int fwdflag, int measurelvl);
void delete_fft_complex_double_plan_cpu(fft_complex_double_plan *self);



/*
  Functions to execute FFTW plans
*/
void execute_complex_single_fft_cpu(complex_vector_single_t *output, complex_vector_single_t *input,
				fft_complex_single_plan *plan);
void execute_real_single_forward_fft_cpu(complex_vector_single_t *output, real_vector_single_t *input,
				fft_real_single_plan *plan);
void execute_real_single_reverse_fft_cpu(real_vector_single_t *output, complex_vector_single_t *input,
				fft_real_single_plan *plan);

void execute_complex_double_fft_cpu(complex_vector_double_t *output, complex_vector_double_t *input,
				fft_complex_double_plan *plan);
void execute_real_double_forward_fft_cpu(complex_vector_double_t *output, real_vector_double_t *input,
				fft_real_double_plan *plan);
void execute_real_double_reverse_fft_cpu(real_vector_double_t *output, complex_vector_double_t *input,
				fft_real_double_plan *plan);

#endif /* FFTCPU_H */
