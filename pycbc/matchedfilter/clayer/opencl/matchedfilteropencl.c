#include <stdio.h>
#include "../../../datavector/clayer/opencl/datavectoropencl_types.h"
#include "../../../clayer/opencl/pycbcopencl_types.h"


void gen_snr_opencl(cl_context_t* context,
                   real_vector_single_opencl_t* snr,
                   complex_vector_single_opencl_t* stilde, 
                   complex_vector_single_opencl_t* htilde)
{    

    printf("gen_snr_opencl in C layer called with context_ptr: %p, device_id: %d\n", 
           context, context->device_id);
    
    //context->set_error(4);
    
    return;

}
