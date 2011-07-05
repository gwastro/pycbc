#include <stdio.h>
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"

// 
// Prototype implementation of the ffts for segmenting straindata 
//

void* fftw_generate_plan(unsigned long length, real_vector_single_t* in_tmp,
                    complex_vector_single_t* out_tmp, char* sign, char* style)
{
    void* plan;
    
    // for testing 3 * 4 byte buffer as "plan" prototype object
    plan = calloc( 3 , 4 );
    
    printf("fftw_generate_plan: length= %ld, in_tmp = %p, out_tmp = %p, sign= %s, style= %s ==> plan= %p\n",
           length, in_tmp, out_tmp, sign, style, plan);
    
    return plan;
    
}

int fftw_transform_segments(void* plan, real_vector_single_t* in_buf, 
                            unsigned long input_buf_offset,
                            complex_vector_single_t* out_buf)
{

    printf("fftw_transform_segments: plan= %p, in_buf + offset = %p, out_buf = %p\n",
           plan, in_buf->data + input_buf_offset, out_buf->data);
    
    // free(plan); // temporarily for testing
    // The plan lives in python as:
    // <Swig Object of type 'void *' at 0x1004ebfc0>
    // When will it be destroyed ??? 
    
    return 0;
}



