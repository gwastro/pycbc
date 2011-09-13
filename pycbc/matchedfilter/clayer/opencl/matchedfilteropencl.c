#include <stdio.h>
#include "../../../datavector/clayer/opencl/datavectoropencl_types.h"
#include "../../../clayer/opencl/pycbcopencl_types.h"
#include "matchedfilteropencl_types.h"

// constructor/destructor of matched filter ------------------------------------

matched_filter_opencl_t* new_matched_filter_opencl_t(void)
{
    matched_filter_opencl_t* c;
    
    c = (matched_filter_opencl_t*) malloc( sizeof(matched_filter_opencl_t) );
    
    c->nothing_to_define_yet = 1;
    
    return c;
}

void delete_matched_filter_opencl_t( matched_filter_opencl_t* p )
{
    //free( p-> ... whatever was allocated );
    free( p );
}


// processing functions of matched filter --------------------------------------

void gen_snr_opencl(cl_context_t* context,
                   real_vector_single_opencl_t* snr,
                   complex_vector_single_opencl_t* stilde, 
                   complex_vector_single_opencl_t* htilde)
{    

    unsigned dev;

    printf("gen_snr_opencl in C layer called with context->device_id: %d\n", 
            context->device_id);

     dev= context->device_id;

    //context->set_error(1);
    
    return;

}
