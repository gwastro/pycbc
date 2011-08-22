#ifndef MATCHEDFILTERCPU_PROTOTYPES_H
#define MATCHEDFILTERCPU_PROTOTYPES_H

#include <stdlib.h>

//prototypes of all functions that swig wraps to methods
int gen_snr_cpu(complex_vector_single_t* stilde, 
                complex_vector_single_t* htilde,
                real_vector_single_t* snr);

#endif /* MATCHEDFILTERCPU_PROTOTYPES_H */
