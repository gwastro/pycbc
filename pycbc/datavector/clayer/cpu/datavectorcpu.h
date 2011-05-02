#ifndef DATAVECTORCPU_H
#define DATAVECTORCPU_H

#include <stdlib.h>


enum cbc_memory_meta_types_t {
    gpu_cuda_global_memory,
    gpu_cuda_constant_memory,
    gpu_cuda_texture_memory,
    gpu_opencl_global_memory,
    gpu_opencl_constant_memory,
    gpu_opencl_local_memory,
    gpu_opencl_private_memory,
    cpugpu_cuda_zero_latency_memory,
    cpu_generic_memory,
    cpu_pinned_memory,
    cpu_fftw_aligned_memory,
    non_specified_memory
};

typedef struct
{
    unsigned long int t_start;
    double dx;
    unsigned int vector_length;
    size_t element_size_bytes;
    enum cbc_memory_meta_types_t memory_type;
    void *data;
}
real_vector_t;


#endif /* DATAVECTORCPU_H */
