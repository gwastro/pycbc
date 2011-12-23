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

#include <except.h> /* For PyCBC exception handling */
#include "fftw_private.h"
#include "fftw.h"
#include "datavectorcpu.h"
/* Note that fftw.h includes <complex.h> and <fftw3.h> (in that order, which is important!) */

fftw_real_single_plan *new_fftw_real_single_plan(unsigned long size,
						 int fwdflag,
						 int measurelvl)
{
  fftw_real_single_plan *self;
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

  self = (fftw_real_single_plan *) malloc(sizeof(fftw_real_single_plan));
  if (!self) {
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate space for real single plan\n");
    return NULL;
  }

  fftwf_import_system_wisdom();

  scratchin  = (float *) fftwf_malloc(size*sizeof(float));
  scratchout = (float *) fftwf_malloc(size*sizeof(float));
  if ((!scratchin) || (!scratchout)){
    free(scratchin);
    free(scratchout);
    free(self);
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space needed to create real single plan\n");
    return NULL;
  }

  self->theplan = fftwf_plan_r2r_1d(size,scratchin,scratchout,
				    fwdflag ? FFTW_R2HC : FFTW_HC2R,
				    flags);
  fftwf_free(scratchin);
  fftwf_free(scratchout);

  if (!(self->theplan)){
    free(self);
    pycbc_throw_exception(PYCBC_RUNTIME_ERROR,"FFTW3 unable to create real single plan\n");
    return NULL;
  }

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}


fftw_real_double_plan *new_fftw_real_double_plan(unsigned long size,
						 int fwdflag,
						 int measurelvl)
{
  fftw_real_double_plan *self;
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

  self = (fftw_real_double_plan *) malloc(sizeof(fftw_real_double_plan));
  if (!self) {
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate space for real double plan\n");
    return NULL;
  }

  fftw_import_system_wisdom();

  scratchin  = (double *) fftw_malloc(size*sizeof(double));
  scratchout = (double *) fftw_malloc(size*sizeof(double));
  if ((!scratchin) || (!scratchout)){
    free(scratchin);
    free(scratchout);
    free(self);
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space needed to create real double plan\n");
    return NULL;
  }

  self->theplan = fftw_plan_r2r_1d(size,scratchin,scratchout,
				   fwdflag ? FFTW_R2HC : FFTW_HC2R,
				   flags);
  fftw_free(scratchin);
  fftw_free(scratchout);

  if (!(self->theplan)){
    free(self);
    pycbc_throw_exception(PYCBC_RUNTIME_ERROR,"FFTW3 unable to create real double plan\n");
    return NULL;
  }

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}

fftw_complex_single_plan *new_fftw_complex_single_plan(unsigned long size,
						       int fwdflag,
						       int measurelvl)
{
  fftw_complex_single_plan *self;
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

  self = (fftw_complex_single_plan *) malloc(sizeof(fftw_complex_single_plan));
  if (!self) {
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate space for complex single plan\n");
    return NULL;
  }

  fftwf_import_system_wisdom();

  scratchin  = (fftwf_complex *) fftwf_malloc(size*sizeof(fftwf_complex));
  scratchout = (fftwf_complex *) fftwf_malloc(size*sizeof(fftwf_complex));
  if ((!scratchin) || (!scratchout)){
    free(scratchin);
    free(scratchout);
    free(self);
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space needed to create complex single plan\n");
    return NULL;
  }

  self->theplan = fftwf_plan_dft_1d(size,scratchin,scratchout,
				    fwdflag ? FFTW_FORWARD : FFTW_BACKWARD,
				    flags);

  fftwf_free(scratchin);
  fftwf_free(scratchout);

  if (!(self->theplan)){
    free(self);
    pycbc_throw_exception(PYCBC_RUNTIME_ERROR,"FFTW3 unable to create complex single plan\n");
    return NULL;
  }

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}

fftw_complex_double_plan *new_fftw_complex_double_plan(unsigned long size,
						       int fwdflag,
						       int measurelvl){
  fftw_complex_double_plan *self;
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

  self = (fftw_complex_double_plan *) malloc(sizeof(fftw_complex_double_plan));
  if (!self) {
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate space for complex double plan\n");
    return NULL;
  }

  fftw_import_system_wisdom();

  scratchin  = (fftw_complex *) fftw_malloc(size*sizeof(fftw_complex));
  scratchout = (fftw_complex *) fftw_malloc(size*sizeof(fftw_complex));
  if ((!scratchin) || (!scratchout)){
    free(scratchin);
    free(scratchout);
    free(self);
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space needed to create complex double plan\n");
    return NULL;
  }

  self->theplan = fftw_plan_dft_1d(size,scratchin,scratchout,
				   fwdflag ? FFTW_FORWARD : FFTW_BACKWARD,
				   flags);

  fftw_free(scratchin);
  fftw_free(scratchout);

  if (!(self->theplan)){
    free(self);
    pycbc_throw_exception(PYCBC_RUNTIME_ERROR,"FFTW3 unable to create complex double plan\n");
    return NULL;
  }

  self->size = size;
  self->fwdflag = fwdflag;
  return self;
}

void delete_fftw_real_single_plan(fftw_real_single_plan *self){
  fftwf_destroy_plan(self->theplan);
  free(self);
  return;
}

void delete_fftw_complex_single_plan(fftw_complex_single_plan *self){
  fftwf_destroy_plan(self->theplan);
  free(self);
  return;
}

void delete_fftw_real_double_plan(fftw_real_double_plan *self){
  fftw_destroy_plan(self->theplan);
  free(self);
  return;
}

void delete_fftw_complex_double_plan(fftw_complex_double_plan *self){
  fftw_destroy_plan(self->theplan);
  free(self);
  return;
}

/* Functions to execute single precision CPU FFTs with FFTW */

void execute_complex_single_fft(complex_vector_single_cpu_t *output,
				complex_vector_single_cpu_t *input,
				fftw_complex_single_plan *plan)
{
  if ( !output || !input ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector to execute_complex_single_fft().\n");
    return;
  }

  if ( !(output->data) || !(input->data) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector->data to execute_complex_single_fft().\n");
    return;
  }

  if ( output->data == input->data ){ /* No in place transforms */
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"No in-place transforms for execute_complex_single_fft().\n");
    return;
  }

  if ( !plan ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty plan given to execute_complex_single_fft().\n");
    return;
  }

  if ( (output->meta_data.vector_length != plan->size) ||
       (input->meta_data.vector_length != plan->size) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Input or output detavector length does not match plan length in execute_complex_single_fft().\n");
    return;
  }

  fftwf_execute_dft(plan->theplan,(fftwf_complex *)input->data,(fftwf_complex *)output->data);

  output->meta_data.delta_x = 1.0/(input->meta_data.vector_length * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  return;
}

void execute_real_single_forward_fft(complex_vector_single_cpu_t *output,
				     real_vector_single_cpu_t *input,
				     fftw_real_single_plan *plan)
{
  unsigned long k;
  float *tmp;

  if ( !output || !input ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector to execute_real_single_forward_fft().\n");
    return;
  }

  if ( !(output->data) || !(input->data) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector->data to execute_real_single_forward_fft().\n");
    return;
  }

  if ( output->data == input->data ){ /* No in place transforms */
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"No in-place transforms for execute_real_single_forward_fft().\n");
    return;
  }

  if ( !plan ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty plan given to execute_real_single_forward_fft().\n");
    return;
  }

  if ( !(plan->fwdflag) ) {
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Called execute_real_single_forward_fft() withOUT fwdflag set.\n");
    return;
  }

  if ( (input->meta_data.vector_length != plan->size) ||
       (output->meta_data.vector_length != (plan->size/2+1))){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Input or output detavector length does not match plan length in execute_real_single_forward_fft().\n");
    return;
  }

  tmp = (float *) fftwf_malloc(plan->size *sizeof(*tmp));
  if (!tmp){
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space in execute_real_single_forward_fft().\n");
    return;
  }


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

void execute_real_single_reverse_fft(real_vector_single_cpu_t *output,
				     complex_vector_single_cpu_t *input,
				     fftw_real_single_plan *plan){
  unsigned long k;
  float *tmp;

  if ( !output || !input ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector to execute_real_single_reverse_fft().\n");
    return;
  }

  if ( !(output->data) || !(input->data) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector->data to execute_real_single_reverse_fft().\n");
    return;
  }

  if ( output->data == input->data ){ /* No in place transforms */
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"No in-place transforms for execute_real_single_reverse_fft().\n");
    return;
  }

  if ( !plan ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty plan given to execute_real_single_reverse_fft().\n");
    return;
  }

  if (plan->fwdflag) {
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Called execute_real_single_reverse_fft() WITH fwdflag set.\n");
    return;
  }

  if ( (output->meta_data.vector_length != plan->size) ||
       (input->meta_data.vector_length != (plan->size/2+1))){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Input or output detavector length does not match plan length in execute_real_single_reverse_fft().\n");
    return;
  }

  if ( cimagf(input->data[0]) != 0.0 ) { // Imag part DC must be zero
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Imaginary part of input DC nonzero in execute_real_single_reverse_fft().\n");
    return;
  }

  if ( ((plan->size % 2) == 0) && cimagf(input->data[plan->size]) != 0.0 ) { // Imag part Nyquist must be zero
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Imaginary part of input Nyquist nonzero for even length in execute_real_single_reverse_fft().\n");
    return;
  }

  tmp = (float *) fftwf_malloc(plan->size *sizeof(*tmp));
  if (!tmp){
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space in execute_real_single_reverse_fft().\n");
    return;
  }

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

void execute_complex_double_fft(complex_vector_double_cpu_t *output,
				complex_vector_double_cpu_t *input,
				fftw_complex_double_plan *plan){

  if ( !output || !input ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector to execute_complex_double_fft().\n");
    return;
  }

  if ( !(output->data) || !(input->data) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector->data to execute_complex_double_fft().\n");
    return;
  }

  if ( output->data == input->data ){ /* No in-place transforms */
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"No in-place transforms for execute_complex_double_fft().\n");
    return;
  }

  if ( !plan ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty plan given to execute_complex_double_fft().\n");
    return;
  }

  if ( (output->meta_data.vector_length != plan->size) ||
       (input->meta_data.vector_length != plan->size) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Input or output detavector length does not match plan length in execute_complex_double_fft().\n");
    return;
  }

  fftw_execute_dft(plan->theplan,(fftw_complex *)input->data,(fftw_complex *)output->data);

  output->meta_data.delta_x = 1.0/(input->meta_data.vector_length * input->meta_data.delta_x);
  output->meta_data.start = input->meta_data.start;

  return;
}

void execute_real_double_forward_fft(complex_vector_double_cpu_t *output,
				     real_vector_double_cpu_t *input,
				     fftw_real_double_plan *plan){
  unsigned long k;
  double *tmp;

  if ( !output || !input ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector to execute_real_double_forward_fft().\n");
    return;
  }

  if ( !(output->data) || !(input->data) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector->data to execute_real_double_forward_fft().\n");
    return;
  }

  if ( output->data == input->data ){ /* No in place transforms */
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"No in-place transforms for execute_real_double_forward_fft().\n");
    return;
  }

  if ( !plan ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty plan given to execute_real_double_forward_fft().\n");
    return;
  }

  if ( !(plan->fwdflag) ) {
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Called execute_real_double_forward_fft() withOUT fwdflag set.\n");
    return;
  }

  if ( (input->meta_data.vector_length != plan->size) ||
       (output->meta_data.vector_length != (plan->size/2+1))){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Input or output detavector length does not match plan length in execute_real_double_forward_fft().\n");
    return;
  }

  tmp = (double *) fftw_malloc(plan->size *sizeof(*tmp));
  if (!tmp){
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space in execute_real_single_forward_fft().\n");
    return;
  }

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

void execute_real_double_reverse_fft(real_vector_double_cpu_t *output,
				     complex_vector_double_cpu_t *input,
				     fftw_real_double_plan *plan){
  unsigned long k;
  double *tmp;

  if ( !output || !input ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector to execute_real_double_reverse_fft().\n");
    return;
  }

  if ( !(output->data) || !(input->data) ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty input or output datavector->data to execute_real_double_reverse_fft().\n");
    return;
  }

  if ( output->data == input->data ){ /* No in place transforms */
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"No in-place transforms for execute_real_double_reverse_fft().\n");
    return;
  }

  if ( !plan ){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Empty plan given to execute_real_double_reverse_fft().\n");
    return;
  }

  if ( plan->fwdflag ) {
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Called execute_real_double_reverse_fft() WITH fwdflag set.\n");
    return;
  }

  if ( (output->meta_data.vector_length != plan->size) ||
       (input->meta_data.vector_length != (plan->size/2+1))){
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Input or output detavector length does not match plan length in execute_real_double_reverse_fft().\n");
    return;
  }

  if ( cimag(input->data[0]) != 0.0 ) { // Imag part DC must be zero
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Imaginary part of input DC nonzero in execute_real_double_reverse_fft().\n");
    return;
  }

  if ( ((plan->size % 2) == 0) && cimag(input->data[plan->size]) != 0.0 ) { // Imag part Nyquist must be zero
    pycbc_throw_exception(PYCBC_VALUE_ERROR,"Imaginary part of input Nyquist nonzero for even length in execute_real_double_reverse_fft().\n");
    return;
  }

  tmp = (double *) fftw_malloc(plan->size *sizeof(*tmp));
  if (!tmp){
    pycbc_throw_exception(PYCBC_MEMORY_ERROR,"Could not allocate scratch space in execute_real_double_reverse_fft().\n");
    return;
  }

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

