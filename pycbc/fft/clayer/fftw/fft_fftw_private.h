#ifndef FFT_FFTW_PRIVATE_H
#define FFT_FFTW_PRIVATE_H

#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "fft_fftw.h"

/*
 Prototypes of constructors/destructors.  
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

#endif  /* #ifndef FFT_FFTW_PRIVATE_H */
