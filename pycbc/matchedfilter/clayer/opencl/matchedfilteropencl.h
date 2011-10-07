#ifndef MATCHEDFILTEROPENCL_TYPES_H
#define MATCHEDFILTEROPENCL_TYPES_H

#include <stdlib.h>
#include <CL/opencl.h>

typedef struct
{
    cl_program program;

    cl_mem cl_stilde;
    cl_mem cl_htilde;

    cl_kernel gpu_snr_product;
}
matched_filter_opencl_t;

#endif /* MATCHEDFILTEROPENCL_TYPES_H */
