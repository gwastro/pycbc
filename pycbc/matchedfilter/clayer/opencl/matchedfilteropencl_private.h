#ifndef MATCHEDFILTEROPENCL_PROTOTYPES_H
#define MATCHEDFILTEROPENCL_PROTOTYPES_H

#include <stdlib.h>


// prototypes of all methodes that will extend pure c typedefs
extern "C" matched_filter_opencl_t* new_matched_filter_opencl_t(cl_context_t*);
void delete_matched_filter_opencl_t( matched_filter_opencl_t* );

//prototypes of all functions that swig wraps to methods
extern "C" void gen_snr_opencl(cl_context_t* context,
                    matched_filter_opencl_t * matchedfilter,
                    real_vector_single_opencl_t* snr,
                    complex_vector_single_opencl_t* stilde, 
                    complex_vector_single_opencl_t* htilde);

#endif /* MATCHEDFILTEROPENCL_PROTOTYPES_H */
