#ifndef FFT_FFTW_PRIVATE_H
#define FFT_FFTW_PRIVATE_H

#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>

/* Struct common for single-precision plans */
#define FFT_PLAN_STRUCT_DATA_FFTWF \
  fftwf_plan    *theplan; \
  unsigned long size; \
  int           fwdflag; \

/* Struct common for double-precision plans */
#define FFT_PLAN_STRUCT_DATA_FFTW \
  fftw_plan    *theplan; \
  unsigned long size; \
  int           fwdflag; \

/* Now the actual structs our functions use */

typedef struct {
  //FFT_PLAN_STRUCT_DATA_FFTWF
  fftwf_plan    *theplan; 
  unsigned long size; 
  int           fwdflag; 
} fft_real_single_plan_fftw;

typedef struct {
FFT_PLAN_STRUCT_DATA_FFTW
} fft_real_double_plan_fftw;

typedef struct {
FFT_PLAN_STRUCT_DATA_FFTWF
} fft_complex_single_plan_fftw;

typedef struct {
FFT_PLAN_STRUCT_DATA_FFTW
} fft_complex_double_plan_fftw;


#endif  /* #ifndef FFT_FFTW_PRIVATE_H */
