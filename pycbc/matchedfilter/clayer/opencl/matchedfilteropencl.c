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
#include <string.h>
#include <iostream>
#include <fstream>
#include "../../../datavector/clayer/opencl/datavectoropencl_types.h"
#include "../../../clayer/opencl/pycbcopencl_types.h"
#include "../../../clayer/opencl/gpu_inspiral_gpuutils.h"
#include "matchedfilteropencl_types.h"

using namespace std;
// constructor/destructor of matched filter ------------------------------------

extern "C" matched_filter_opencl_t* new_matched_filter_opencl_t(cl_context_t * context)
{
    matched_filter_opencl_t* c;

//Auxiliary variables for getting the build log
    size_t deviceNum;
    size_t LogSize;
    char * BuildLog = NULL;
    cl_device_id *cdDevices;

    int err;

    c = (matched_filter_opencl_t*) malloc( sizeof(matched_filter_opencl_t) );

//Loading the kernel source
    std::ifstream kernel_source_file("./gpu_inspiral_matchedfilter.cl");
    std::string prog( std::istreambuf_iterator<char>(kernel_source_file), (std::istreambuf_iterator<char>()));
    const char * kernel_source = prog.c_str();
    size_t kernel_source_size = prog.size();


//Creating the program from the kernel source
    c->program = clCreateProgramWithSource(context->context, 1, &kernel_source, &kernel_source_size, &err);
    if (gpuinsp_checkError(err,"Creating program from source file") !=0)  
      goto cleanup;

//Building the kernel
    err = clBuildProgram(c->program, 0 , NULL, NULL, NULL, NULL);

//Trying to fetch the build log if build failed
    if (err != CL_BUILD_SUCCESS) {

          err = clGetContextInfo(context->context, CL_CONTEXT_DEVICES, 0, NULL, &deviceNum);
          cdDevices = (cl_device_id *)malloc(deviceNum * sizeof(cl_device_id));
          err  = clGetContextInfo(context->context, CL_CONTEXT_DEVICES, deviceNum * sizeof(cl_device_id), cdDevices, NULL);


          err = clGetProgramBuildInfo(c -> program, cdDevices[0], CL_PROGRAM_BUILD_LOG , 0 , NULL, &LogSize);
          BuildLog =  (char *) malloc(LogSize);
          err = clGetProgramBuildInfo(c -> program, cdDevices[0], CL_PROGRAM_BUILD_LOG , LogSize, BuildLog, NULL);

          printf("The build log message is: %s \n",BuildLog);
          goto cleanup;
        }

//Ok, now we can create the kernels

    printf("in mf constructor c->program %p  ", c->program);


     c->gpu_snr_product = clCreateKernel(c->program, "gpuSnrProduct", &err);

    if(gpuinsp_checkError(err, "clCreateKernel(gpuSnrProduct)") != CL_SUCCESS) goto cleanup;

//Normal termination
    return c;

//Freeing up the temporaly pointers and quiting
cleanup:
    if (c) free(c);
    if (cdDevices) free(cdDevices);
    if (BuildLog) free(BuildLog);
    return NULL;

}

void delete_matched_filter_opencl_t( matched_filter_opencl_t* p )
{
    free( p );
}


// processing functions of matched filter --------------------------------------

extern "C" void gen_snr_opencl(cl_context_t* context,
                   matched_filter_opencl_t * matchedfilter,
                   real_vector_single_opencl_t* snr,
                   complex_vector_single_opencl_t* stilde,
                   complex_vector_single_opencl_t* htilde)
{

    size_t local_size   = 256;
    size_t global_size  = stilde->meta_data.vector_length/2;
    int err          = 0;

    err = clSetKernelArg(matchedfilter->gpu_snr_product, 0, sizeof(size_t), &stilde->meta_data.vector_length);
    err |= clSetKernelArg(matchedfilter->gpu_snr_product, 1, sizeof(cl_mem), &stilde->real_data);
    err |= clSetKernelArg(matchedfilter->gpu_snr_product, 2, sizeof(cl_mem), &stilde->imag_data);
    err |= clSetKernelArg(matchedfilter->gpu_snr_product, 3, sizeof(cl_mem), &htilde->real_data);
    err |= clSetKernelArg(matchedfilter->gpu_snr_product, 4, sizeof(cl_mem), &htilde->imag_data);
    err |= clSetKernelArg(matchedfilter->gpu_snr_product, 5, sizeof(cl_mem), &snr->data);
    if(gpuinsp_checkError(err, "clSetKernelArg(gpu_snr_product)") != CL_SUCCESS) return;


    err = clEnqueueNDRangeKernel(context->kernel_queue, matchedfilter->gpu_snr_product, 1, 0, &global_size, &local_size, 0, NULL, NULL);
    if (gpuinsp_checkError(err, "clEnqueueNDRangeKernel(gpu_snr_product)") != CL_SUCCESS) return;
    clFinish(context->kernel_queue);

}
