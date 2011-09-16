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
    if (gpuinsp_checkError(err,"Creating program from source file") !=0)  goto cleanup;

//Building the kernel
    err = clBuildProgram(c->program, 0 , NULL, NULL, NULL, NULL);
    if (gpuinsp_checkError(err,"Building program source") !=0)  goto cleanup;

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

    printf("gen_snr_opencl in C layer called with context->device_id: %d\n", context->device_id);

/*
   err = clSetKernelArg(_plan->gpu_gen_templates, 0, sizeof(int), &_plan->signal_batch_size);
        err |= clSetKernelArg(_plan->gpu_gen_templates, 1, sizeof(int), &_plan->signal_batch_number);
        err |= clSetKernelArg(_plan->gpu_gen_templates, 2, sizeof(int), &_plan->flow);
        err |= clSetKernelArg(_plan->gpu_gen_templates, 3, sizeof(cl_mem), &_plan->cl_mass_array);
        err |= clSetKernelArg(_plan->gpu_gen_templates, 4, sizeof(int), &_plan->template_length);
       	err |= clSetKernelArg(_plan->gpu_gen_templates, 5, sizeof(cl_mem), &_plan->cl_signals_batch_real[i]);
        err |= clSetKernelArg(_plan->gpu_gen_templates, 6, sizeof(cl_mem), &_plan->cl_signals_batch_imag[i]);
        if(checkError(err, "clSetKernelArg(gpu_gen_templates)") != CL_SUCCESS) return err;

        local_size = 256;
        global_size = _plan->template_length*_plan->signal_batch_size/2;
        err = clEnqueueNDRangeKernel(_plan->kernel_queue, _plan->gpu_gen_templates, 1, 0, &global_size, &local_size, 0, NULL, NULL);
        if (checkError(err, "clEnqueueNDRangeKernel(gpu_gen_templates)") != CL_SUCCESS) return err;
    }
    clFinish(_plan->kernel_queue);
*/

    return;

}
