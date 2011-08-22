#ifndef MATCHEDFILTERCPU_PROTOTYPES_H
#define MATCHEDFILTERCPU_PROTOTYPES_H

#include <stdlib.h>

//prototypes of all functions that swig wraps to methods
void gen_snr_cpu(cpu_context_t* context,
                real_vector_single_t* snr,
                complex_vector_single_t* stilde, 
                complex_vector_single_t* htilde);

#endif /* MATCHEDFILTERCPU_PROTOTYPES_H */
