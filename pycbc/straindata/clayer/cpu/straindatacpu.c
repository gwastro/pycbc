#include <stdio.h>
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"


// 
// Prototyping implementation of the ffts for segmenting straindata for 
// Josh (aka Mr. FFTW) 
//


void* fftw_generate_plan(unsigned long length, real_vector_single_t* in_tmp,
                       complex_vector_single_t* out_tmp, char* sign, char* style)
{
    void* plan;
    
    plan = calloc( 3 , 4 );  // testing 3 * 4 byte buffer as "plan" object
    
    printf("fftw_generate_plan: length= %d, in_tmp = %p, out_tmp = %p, sign= %s, style= %s ==> plan= %p\n",
           length, in_tmp, out_tmp, sign, style, plan);
    
    return plan;
    
}

int fftw_transform_segments(void* plan,     
                            real_vector_single_t* in_buf, 
                            complex_vector_single_t* out_buf)
{

    printf("fftw_transform_segments called w/: plan= %p, in_buf = %p, out_buf = %p\n",
           plan, in_buf, out_buf);
    
    free(plan);
    
    return 0;
}



