// Copyright (C) 2011 Gergely Deberczeni (Gergely.Debreczeni@rmki.kfki.hu),
//                    Karsten Wiesner  (karsten.wiesner@aei.mpg.de)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


//
// =============================================================================
//
//                                   Preamble
//
// =============================================================================
//
// The matched_filter_opencl_t initialization, adn definitions


#include <stdio.h>
#include "../../../datavector/clayer/opencl/datavectoropencl_types.h"
#include "../../../clayer/opencl/pycbcopencl_types.h"
#include "matchedfilteropencl_types.h"

// constructor/destructor of matched filter ------------------------------------

matched_filter_opencl_t* new_matched_filter_opencl_t(void)
{
    matched_filter_opencl_t* c;
    int err;

    c = (matched_filter_opencl_t*) malloc( sizeof(matched_filter_opencl_t) );

//Loading the kernel source
    std::ifstream kernel_source_file("./gpu_inspiral_matchedfilter.cl");
    std::string prog( std::istreambuf_iterator<char>(kernel_source_file), (std::istreambuf_iterator<char>()));
    const char* kernel_source = prog.c_str();
    size_t kernel_source_size = prog.size();

//Creating the program from the kernel source
    c->program = clCreateProgramWithSource(context, 1, &kernel_source, &kernel_source_size, &err);
    if (gpuinsp_checkError(err,"Creating program from source file") !=0)  goto cleanup;

//Building the kernel
    err = clBuildProgram(c->program, 0 , NULL, NULL, NULL, NULL);

//Trying to fetch the build log if build failed
    if (err != CL_BUILD_SUCCESS) {
          size_t deviceNum;
          cl_device_id *cdDevices;

          err = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceNum);
          cdDevices = (cl_device_id *)malloc(deviceNum * sizeof(cl_device_id));
          err  = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceNum * sizeof(cl_device_id), cdDevices, NULL);

          size_t LogSize;

          err = clGetProgramBuildInfo(c -> program, cdDevices[0], CL_PROGRAM_BUILD_LOG , 0 , NULL, &LogSize);
          char *  LogMessage =  (char *) malloc(LogSize);
          err = clGetProgramBuildInfo(c -> program, cdDevices[0], CL_PROGRAM_BUILD_LOG , LogSize, LogMessage, NULL);

          printf("The build log message is: %s \n",LogMessage);
          goto cleanup;
        }
//Ok, now we can create the kernels
    c->gpu_snr_product = clCreateKernel(c->program, "gpuSnrProduct", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuSnrProduct)") != CL_SUCCESS) goto cleanup;
    c->gpu_snr_normalize = clCreateKernel(c->program, "gpuSnrNormalize", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuSnrNormalize)") != CL_SUCCESS) goto cleanup;
    c->gpu_template_product = clCreateKernel(c->program, "gpuTemplateProduct", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuTemplateProduct)") != CL_SUCCESS) goto cleanup;
    c->gpu_template_product_chi2 = clCreateKernel(c->program, "gpuTemplateProductChi2", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuTemplateProductChi2)") != CL_SUCCESS) goto cleanup;
    c->gpu_init_buffer_real = clCreateKernel(c->program, "gpuInitBufferReal", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuInitBufferReal)") != CL_SUCCESS) goto cleanup;
    c->gpu_init_peak_num_buffers  = clCreateKernel(c->program, "gpuInitPeakNumBuffers", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuInitPeakNumBuffers)") != CL_SUCCESS) goto cleanup;
    c->gpu_sumreduce_variance = clCreateKernel(c->program, "gpuSumReduceVariance", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuSumReduceVariance)") != CL_SUCCESS) goto cleanup;
    c->gpu_sumreduce_variance_final = clCreateKernel(c->program, "gpuSumReduceVarianceFinal", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuSumReduceVarianceFinal)") != CL_SUCCESS) goto cleanup;
    c->gpu_sumreduce_chi2 = clCreateKernel(c->program, "gpuSumReduceChi2", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuSumReduceChi2)") != CL_SUCCESS) goto cleanup;
    c->gpu_sumreduce_chi2_final = clCreateKernel(c->program, "gpuSumReduceChi2Final", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuSumReduceChi2Final)") != CL_SUCCESS) goto cleanup;
    c->gpu_peak_search = clCreateKernel(c->program, "gpuPeakSearch", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuPeakSearch)") != CL_SUCCESS) goto cleanup;
    c->gpu_gen_templates = clCreateKernel(c->program, "gpuGenerateTemplates", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuGenerateTemplates)") != CL_SUCCESS) goto cleanup;
    c->gpu_hamm_weight = clCreateKernel(c->program, "gpuHammingWeight", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuHammingWeight)") != CL_SUCCESS) goto cleanup;
    c->gpu_avrg_psd = clCreateKernel(c->program, "gpuAvrgPSD", &err);
    if(gpuinsp_checkError(err, "clCreateKernel(gpuAvrgPSD)") != CL_SUCCESS) goto cleanup;

//Normal termination
    return c;

//Freeing up the temporaly pointers and quiting
cleanup:
    if (c) free(c);
    if (cdDevices) free(cdDevices);
    if (LogMessage) free(LogMessage);
    return NULL;
}

void delete_matched_filter_opencl_t( matched_filter_opencl_t* p )
{
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
