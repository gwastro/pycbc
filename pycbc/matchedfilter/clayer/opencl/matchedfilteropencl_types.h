#ifndef MATCHEDFILTEROPENCL_TYPES_H
#define MATCHEDFILTEROPENCL_TYPES_H

#include <stdlib.h>
#include <CL/opencl.h>

typedef struct
{
    cl_program program
    char * kernel_source;
    size_t kernel_source_size;
    cl_kernel gpu_snr_product;
    cl_kernel gpu_snr_normalize;
    cl_kernel gpu_template_product;
    cl_kernel gpu_template_product_chi2;
    cl_kernel gpu_init_buffer_real;
    cl_kernel gpu_init_peak_num_buffers;
    cl_kernel gpu_sumreduce_variance;
    cl_kernel gpu_sumreduce_variance_final;
    cl_kernel gpu_sumreduce_chi2;
    cl_kernel gpu_sumreduce_chi2_final;
    cl_kernel gpu_peak_search;
    cl_kernel gpu_gen_templates;
    cl_kernel gpu_hamm_weight;
    cl_kernel gpu_avrg_psd;
}
matched_filter_opencl_t;

#endif /* MATCHEDFILTEROPENCL_TYPES_H */
