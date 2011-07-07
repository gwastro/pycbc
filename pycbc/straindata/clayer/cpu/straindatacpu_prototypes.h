#ifndef STRAINDATACPU_PROTOTYPES_H
#define STRAINDATACPU_PROTOTYPES_H

#include <stdlib.h>

//prototypes of all functions that swig wraps to methods

void* fftw_generate_plan(unsigned long length, real_vector_single_t* in_tmp,
                    complex_vector_single_t* out_tmp, char* sign, char* style);

int fftw_transform_segments(void* plan, real_vector_single_t* in_buf, 
                            unsigned long input_buf_offset, 
                            complex_vector_single_t* out_buf);

int frame_cpp_read_frames(real_vector_double_t* out_buf, char* channel_name, 
                          unsigned long gps_start_time, unsigned long gps_end_time, 
                          char* cache_url);

#endif /* STRAINDATACPU_PROTOTYPES_H */
