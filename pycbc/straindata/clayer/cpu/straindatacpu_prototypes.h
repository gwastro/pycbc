#ifndef STRAINDATACPU_PROTOTYPES_H
#define STRAINDATACPU_PROTOTYPES_H

#include <stdlib.h>

//prototypes of all functions that swig wraps to methods

void* fftw_generate_plan(unsigned long length, real_vector_single_t* in_tmp,
                       complex_vector_single_t* out_tmp, char* sign, char* style);

int fftw_transform_segments(void* plan,
                            real_vector_single_t* in_buf, 
                            complex_vector_single_t* out_buf);


#endif /* STRAINDATACPU_PROTOTYPES_H */
