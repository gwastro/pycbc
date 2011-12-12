#ifndef PYCBC_FFTW_PRIVATE_H
#define PYCBC_FFTW_PRIVATE_H

#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "fftw.h"

/*
 Prototypes of constructors/destructors.
*/

fftw_real_single_plan *new_fftw_real_single_plan(unsigned long size,
						 int fwdflag,
						 int measurelvl);
void delete_fftw_real_single_plan(fftw_real_single_plan *self);

fftw_real_double_plan *new_fftw_real_double_plan(unsigned long size,
						 int fwdflag,
						 int measurelvl);
void delete_fftw_real_double_plan(fftw_real_double_plan *self);

fftw_complex_single_plan *new_fftw_complex_single_plan(unsigned long size,
						       int fwdflag,
						       int measurelvl);
void delete_fftw_complex_single_plan(fftw_complex_single_plan *self);

fftw_complex_double_plan *new_fftw_complex_double_plan(unsigned long size,
						       int fwdflag,
						       int measurelvl);
void delete_fftw_complex_double_plan(fftw_complex_double_plan *self);

#endif  /* #ifndef PYCBC_FFTW_PRIVATE_H */
