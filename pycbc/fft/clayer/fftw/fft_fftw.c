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


// FFTW plan constructor and destructor implementation for pycbc

#include <stdlib.h>
#include <stdio.h>

#include "fft_fftw_private.h"
#include "fft_fftw.h"
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"
/* Note that fft_fftw.h includes <complex.h> and <fftw3.h> (in that order, which is important!) */

fft_real_single_plan_fftw *new_fft_real_single_plan_fftw(unsigned long size,
							 int fwdflag,
							 int measurelvl)
{
  fft_real_single_plan_fftw *self;
  float *scratchin, *scratchout;
  int flags = 0;
  
  switch (measurelvl) {
  case 0:
    flags |= FFTW_ESTIMATE;
    break;
  default: /* Exhaustive */
    flags |= FFTW_EXHAUSTIVE;
    /* fall-through */
  case 2:
    flags |= FFTW_PATIENT;
    /* fall-through */
  case 1:
    flags |= FFTW_MEASURE;
    break;
  }
  
  self = (fft_real_single_plan_fftw *) malloc(sizeof(fft_real_single_plan_fftw));
  fftwf_import_system_wisdom();

  scratchin  = (float *) fftwf_malloc(size*sizeof(float));
  scratchout = (float *) fftwf_malloc(size*sizeof(float));
  self->theplan = fftwf_plan_r2r_1d(size,scratchin,scratchout,
				    fwdflag ? FFTW_R2HC : FFTW_HC2R,
				    flags);
  fftwf_free(scratchin);
  fftwf_free(scratchout);

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}


fft_real_double_plan_fftw *new_fft_real_double_plan_fftw(unsigned long size,
							 int fwdflag,
							 int measurelvl)
{
  fft_real_double_plan_fftw *self;
  double *scratchin, *scratchout;
  int flags = 0;
  
  switch (measurelvl) {
  case 0:
    flags |= FFTW_ESTIMATE;
    break;
  default: /* Exhaustive */
    flags |= FFTW_EXHAUSTIVE;
    /* fall-through */
  case 2:
    flags |= FFTW_PATIENT;
    /* fall-through */
  case 1:
    flags |= FFTW_MEASURE;
    break;
  }
  
  self = (fft_real_double_plan_fftw *) malloc(sizeof(fft_real_double_plan_fftw));
  fftw_import_system_wisdom();

  scratchin  = (double *) fftw_malloc(size*sizeof(double));
  scratchout = (double *) fftw_malloc(size*sizeof(double));
  self->theplan = fftw_plan_r2r_1d(size,scratchin,scratchout,
				   fwdflag ? FFTW_R2HC : FFTW_HC2R,
				   flags);
  fftw_free(scratchin);
  fftw_free(scratchout);

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}

fft_complex_single_plan_fftw *new_fft_complex_single_plan_fftw(unsigned long size,
							       int fwdflag,
							       int measurelvl)
{
  fft_complex_single_plan_fftw *self;
  fftwf_complex *scratchin, *scratchout;
  int flags = 0;
  
  switch (measurelvl) {
  case 0:
    flags |= FFTW_ESTIMATE;
    break;
  default: /* Exhaustive */
    flags |= FFTW_EXHAUSTIVE;
    /* fall-through */
  case 2:
    flags |= FFTW_PATIENT;
    /* fall-through */
  case 1:
    flags |= FFTW_MEASURE;
    break;
  }
  
  self = (fft_complex_single_plan_fftw *) malloc(sizeof(fft_complex_single_plan_fftw));
  fftwf_import_system_wisdom();

  scratchin  = (fftwf_complex *) fftwf_malloc(size*sizeof(fftwf_complex));
  scratchout = (fftwf_complex *) fftwf_malloc(size*sizeof(fftwf_complex));
  self->theplan = fftwf_plan_dft_1d(size,scratchin,scratchout,
				    fwdflag ? FFTW_FORWARD : FFTW_BACKWARD,
				    flags);

  fftwf_free(scratchin);
  fftwf_free(scratchout);

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}

fft_complex_double_plan_fftw *new_fft_complex_double_plan_fftw(unsigned long size,
							       int fwdflag,
							       int measurelvl){
  fft_complex_double_plan_fftw *self;
  fftw_complex *scratchin, *scratchout;
  int flags = 0;
  
  switch (measurelvl) {
  case 0:
    flags |= FFTW_ESTIMATE;
    break;
  default: /* Exhaustive */
    flags |= FFTW_EXHAUSTIVE;
    /* fall-through */
  case 2:
    flags |= FFTW_PATIENT;
    /* fall-through */
  case 1:
    flags |= FFTW_MEASURE;
    break;
  }
  
  self = (fft_complex_double_plan_fftw *) malloc(sizeof(fft_complex_double_plan_fftw));
  fftwf_import_system_wisdom();

  scratchin  = (fftw_complex *) fftw_malloc(size*sizeof(fftw_complex));
  scratchout = (fftw_complex *) fftw_malloc(size*sizeof(fftw_complex));
  self->theplan = fftw_plan_dft_1d(size,scratchin,scratchout,
				   fwdflag ? FFTW_FORWARD : FFTW_BACKWARD,
				   flags);

  fftw_free(scratchin);
  fftw_free(scratchout);

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}

void delete_real_single_plan_fftw(fft_real_single_plan_fftw *self){
  fftwf_destroy_plan(self->theplan);
  free(self);
  return;
}

void delete_complex_single_plan_fftw(fft_complex_single_plan_fftw *self){
  fftwf_destroy_plan(self->theplan);
  free(self);
  return;
}

void delete_real_double_plan_fftw(fft_real_double_plan_fftw *self){
  fftw_destroy_plan(self->theplan);
  free(self);
  return;
}

void delete_complex_double_plan_fftw(fft_complex_double_plan_fftw *self){
  fftw_destroy_plan(self->theplan);
  free(self);
  return;
}

/* Functions to execute single precision CPU FFTs with FFTW */
				
void execute_complex_single_fft_fftw(complex_vector_single_cpu_t *output,
				     complex_vector_single_cpu_t *input,
				     fft_complex_single_plan_fftw *plan)
{
  if ( !output || !input ){
    abort(); /* Need better error handling eventually */
  }

  if ( output->data == input->data ){ /* No in place transforms */
    abort(); /* Need better error handling eventually */
  }

  if ( !plan ){
    abort(); /* Need better error handling eventually */
  }

  if ( (output->meta_data.vector_length != plan->size) || 
       (input->meta_data.vector_length != plan->size) ){
    abort(); /* Need better error handling eventually */
  }

  fftwf_execute_dft(plan->theplan,(fftwf_complex *)input->data,(fftwf_complex *)output->data);

  output->meta_data.delta_x = 1.0/(input->meta_data.vector_length * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  return;
}

void execute_real_single_forward_fft_fftw(complex_vector_single_cpu_t *output,
					  real_vector_single_cpu_t *input,
					  fft_real_single_plan_fftw *plan)
{
  unsigned long k;
  float *tmp;

  if ( !output || !input ){
    abort(); /* Need better error handling eventually */
  }

  if ( !plan ){
    abort(); /* Need better error handling eventually */
  }

  if ( !(plan->fwdflag) ) {
    abort(); /* Need better error handling eventually */
  }

  if ( (input->meta_data.vector_length != plan->size) ||
       (output->meta_data.vector_length != (plan->size/2+1))){
    abort(); /* Need better error handling eventually */
  }

  tmp = (float *) fftwf_malloc(plan->size *sizeof(*tmp));

  fftwf_execute_r2r(plan->theplan,input->data,tmp);

  /* dc component */
  output->data[0] = tmp[0] + 0.0*I;

  /* everything else but Nyquist */

  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
    {
      output->data[k] = tmp[k] + tmp[plan->size - k]*I;
    }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
    {
      output->data[plan->size/2] = tmp[plan->size/2]+0.0*I;
    }

  /* vector step */
  output->meta_data.delta_x = 1.0 / (input->meta_data.vector_length * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  fftwf_free(tmp);
  return;
}

void execute_real_single_reverse_fft_fftw(real_vector_single_cpu_t *output,
					  complex_vector_single_cpu_t *input,
					  fft_real_single_plan_fftw *plan){
  unsigned long k;
  float *tmp;

  if ( !output || !input ){
    abort(); /* Need better error handling eventually */
  }

  if ( !plan ){
    abort(); /* Need better error handling eventually */
  }

  if ( plan->fwdflag ) {
    abort(); /* Need better error handling eventually */
  }

  if ( (output->meta_data.vector_length != plan->size) ||
       (input->meta_data.vector_length != (plan->size/2+1))){
    abort(); /* Need better error handling eventually */
  }

  if ( cimagf(input->data[0]) != 0.0 ) { // Imag part DC must be zero
    abort(); /* Need better error handling eventually */
  }

  if ( ((plan->size % 2) == 0) && cimagf(input->data[plan->size]) != 0.0 ) { // Imag part Nyquist must be zero
    abort(); /* Need better error handling eventually */
  }

  tmp = (float *) fftwf_malloc(plan->size *sizeof(*tmp));

  /* DC component */

  tmp[0] = crealf(input->data[0]);

  /* everything else but Nyquist */

  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
    {
      tmp[k]              = crealf(input->data[k]);
      tmp[plan->size - k] = cimagf(input->data[k]);
    }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
    {
      tmp[plan->size/2] = crealf(input->data[plan->size/2]);
    }

  fftwf_execute_r2r(plan->theplan,tmp,output->data);

  /* vector step */
  output->meta_data.delta_x = 1.0 / (2*(input->meta_data.vector_length-1) * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  fftwf_free(tmp);

  return;
}

/* Functions to execute double precision CPU FFTs with FFTW */
			
void execute_complex_double_fft_fftw(complex_vector_double_cpu_t *output,
				     complex_vector_double_cpu_t *input,
				     fft_complex_double_plan_fftw *plan){

  if ( !output || !input ){
    abort(); /* Need better error handling eventually */
  }

  if ( output->data == input->data ){ /* No in-place transforms */
    abort(); /* Need better error handling eventually */
  }

  if ( !plan ){
    abort(); /* Need better error handling eventually */
  }

  if ( (output->meta_data.vector_length != plan->size) || 
       (input->meta_data.vector_length != plan->size) ){
    abort(); /* Need better error handling eventually */
  }

  fftw_execute_dft(plan->theplan,(fftw_complex *)input->data,(fftw_complex *)output->data);

  output->meta_data.delta_x = 1.0/(input->meta_data.vector_length * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  return;
}

void execute_real_double_forward_fft_fftw(complex_vector_double_cpu_t *output,
					  real_vector_double_cpu_t *input,
					  fft_real_double_plan_fftw *plan){
  unsigned long k;
  double *tmp;

  if ( !output || !input ){
    abort(); /* Need better error handling eventually */
  }

  if ( !plan ){
    abort(); /* Need better error handling eventually */
  }

  if ( !(plan->fwdflag) ) {
    abort(); /* Need better error handling eventually */
  }

  if ( (input->meta_data.vector_length != plan->size) ||
       (output->meta_data.vector_length != (plan->size/2+1))){
    abort(); /* Need better error handling eventually */
  }

  tmp = (double *) fftw_malloc(plan->size *sizeof(*tmp));

  fftw_execute_r2r(plan->theplan,input->data,tmp);

  /* dc component */
  output->data[0] = tmp[0] + 0.0*I;

  /* everything else but Nyquist */

  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
    {
      output->data[k] = tmp[k] + tmp[plan->size - k]*I;
    }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
    {
      output->data[plan->size/2] = tmp[plan->size/2]+0.0*I;
    }

  /* vector step */
  output->meta_data.delta_x = 1.0 / (input->meta_data.vector_length * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  fftw_free(tmp);
  return;
}

void execute_real_double_reverse_fft_fftw(real_vector_double_cpu_t *output,
					  complex_vector_double_cpu_t *input,
					  fft_real_double_plan_fftw *plan){
  unsigned long k;
  double *tmp;

  if ( !output || !input ){
    abort(); /* Need better error handling eventually */
  }

  if ( !plan ){
    abort(); /* Need better error handling eventually */
  }

  if ( plan->fwdflag ) {
    abort(); /* Need better error handling eventually */
  }

  if ( (output->meta_data.vector_length != plan->size) ||
       (input->meta_data.vector_length != (plan->size/2+1))){
    abort(); /* Need better error handling eventually */
  }

  if ( cimag(input->data[0]) != 0.0 ) { // Imag part DC must be zero
    abort(); /* Need better error handling eventually */
  }

  if ( ((plan->size % 2) == 0) && cimag(input->data[plan->size]) != 0.0 ) { // Imag part Nyquist must be zero
    abort(); /* Need better error handling eventually */
  }

  tmp = (double *) fftw_malloc(plan->size *sizeof(*tmp));

  /* DC component */

  tmp[0] = creal(input->data[0]);

  /* everything else but Nyquist */

  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
    {
      tmp[k]              = creal(input->data[k]);
      tmp[plan->size - k] = cimag(input->data[k]);
    }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
    {
      tmp[plan->size/2] = creal(input->data[plan->size/2]);
    }

  fftw_execute_r2r(plan->theplan,tmp,output->data);

  /* vector step */
  output->meta_data.delta_x = 1.0 / (2*(input->meta_data.vector_length-1) * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  fftw_free(tmp);

  return;
}

