#ifndef MATCHEDFILTERCPU_PROTOTYPES_H
#define MATCHEDFILTERCPU_PROTOTYPES_H

#include <stdlib.h>
#include "../../../clayer/cpu/pycbccpu_types.h"


//prototypes of all functions that swig wraps to methods

void gen_snr_cpu(cpu_context_t* context,
		 real_vector_single_t* snr,
		 complex_vector_single_t* stilde, 
		 complex_vector_single_t* htilde,
		 complex_vector_single_t* q,
		 complex_vector_single_t* qtilde,
		 /*complex_fft_plan_t* plan,*/
		 double f_min,
		 double sigma_sq);

void correlate_complex_freq_vectors( complex_vector_single_t* out,
				     complex_vector_single_t* x, 
				     complex_vector_single_t* y, 
				     double f_min);

#endif /* MATCHEDFILTERCPU_PROTOTYPES_H */
