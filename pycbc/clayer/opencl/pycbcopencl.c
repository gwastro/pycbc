// Copyright (C) 2011 Karsten Wiesner
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
// pycbc's constructors and destructors implementation for pycbc

#include <stdio.h>
#include <string.h>
#include "pycbcopencl_types.h"
#include "pycbcopencl_prototypes.h"
#include "gpu_inspiral_gpuutils.h"

unsigned pycbcopencl_error_stash = 0;
char pycbcopencl_error_message[ERROR_STRING_LEN];

int pycbc_opencl_check_error()
{
    //printf("debug: pycbc_opencl_check_error\n");
    return pycbcopencl_error_stash;
}

char* pycbc_opencl_get_error_message()
{
    //printf("debug: pycbc_opencl_get_error_message\n");

    return pycbcopencl_error_message;
    //return pycbcopencl_err_map[pycbcopencl_err_stash];
}

void pycbc_opencl_set_error(int opencl_error_id, char* generic_err_message)
{    
    char opencl_err_message[256];    
    
    switch (opencl_error_id) {
        case CL_SUCCESS:
            strcpy(opencl_err_message,"CL_SUCCESS");
            pycbcopencl_error_stash = PYCBC_NO_ERROR;
            break;
        case CL_DEVICE_NOT_FOUND:
            strcpy(opencl_err_message,"CL_DEVICE_NOT_FOUND");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_DEVICE_NOT_AVAILABLE:
            strcpy(opencl_err_message,"CL_DEVICE_NOT_AVAILABLE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_COMPILER_NOT_AVAILABLE:
            strcpy(opencl_err_message,"CL_COMPILER_NOT_AVAILABLE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            strcpy(opencl_err_message,"CL_MEM_OBJECT_ALLOCATION_FAILURE");
            pycbcopencl_error_stash = PYCBC_MEMORY_ERROR;
            break;
        case CL_OUT_OF_RESOURCES:
            strcpy(opencl_err_message,"CL_OUT_OF_RESOURCES");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_OUT_OF_HOST_MEMORY:
            strcpy(opencl_err_message,"CL_OUT_OF_HOST_MEMORY");
            pycbcopencl_error_stash = PYCBC_MEMORY_ERROR;
            break;
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            strcpy(opencl_err_message,"CL_PROFILING_INFO_NOT_AVAILABLE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_MEM_COPY_OVERLAP:
            strcpy(opencl_err_message,"CL_MEM_COPY_OVERLAP");
            break;
        case CL_IMAGE_FORMAT_MISMATCH:
            strcpy(opencl_err_message,"CL_IMAGE_FORMAT_MISMATCH");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            strcpy(opencl_err_message,"CL_IMAGE_FORMAT_NOT_SUPPORTED");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_BUILD_PROGRAM_FAILURE:
            strcpy(opencl_err_message,"CL_BUILD_PROGRAM_FAILURE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_MAP_FAILURE:
            strcpy(opencl_err_message,"CL_MAP_FAILURE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
            /* These are opencl 1.1 features not yet defined
             case CL_MISALIGNED_SUB_BUFFER_OFFSET:
             strcpy(opencl_err_message,"CL_MISALIGNED_SUB_BUFFER_OFFSET");
             pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
             break;
             case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
             strcpy(opencl_err_message,"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST");
             pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
             break;
             */
        case CL_INVALID_VALUE:
            strcpy(opencl_err_message,"CL_INVALID_VALUE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_DEVICE_TYPE:
            strcpy(opencl_err_message,"CL_INVALID_DEVICE_TYPE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_PLATFORM:
            strcpy(opencl_err_message,"CL_INVALID_PLATFORM");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_DEVICE:
            strcpy(opencl_err_message,"CL_INVALID_DEVICE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_CONTEXT:
            strcpy(opencl_err_message,"CL_INVALID_CONTEXT");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_QUEUE_PROPERTIES:
            strcpy(opencl_err_message,"CL_INVALID_QUEUE_PROPERTIES");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_COMMAND_QUEUE:
            strcpy(opencl_err_message,"CL_INVALID_COMMAND_QUEUE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_HOST_PTR:
            strcpy(opencl_err_message,"CL_INVALID_HOST_PTR");
            pycbcopencl_error_stash = PYCBC_MEMORY_ERROR;
            break;
        case CL_INVALID_MEM_OBJECT:
            strcpy(opencl_err_message,"CL_INVALID_MEM_OBJECT");
            pycbcopencl_error_stash = PYCBC_MEMORY_ERROR;
            break;
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            strcpy(opencl_err_message,"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_IMAGE_SIZE:
            strcpy(opencl_err_message,"CL_INVALID_IMAGE_SIZE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_SAMPLER:
            strcpy(opencl_err_message,"CL_INVALID_SAMPLER");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_BINARY:
            strcpy(opencl_err_message,"CL_INVALID_BINARY");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_BUILD_OPTIONS:
            strcpy(opencl_err_message,"CL_INVALID_BUILD_OPTIONS");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_PROGRAM:
            strcpy(opencl_err_message,"CL_INVALID_PROGRAM");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_PROGRAM_EXECUTABLE:
            strcpy(opencl_err_message,"CL_INVALID_PROGRAM_EXECUTABLE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL_NAME:
            strcpy(opencl_err_message,"CL_INVALID_KERNEL_NAME");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL_DEFINITION:
            strcpy(opencl_err_message,"CL_INVALID_KERNEL_DEFINITION");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL:
            strcpy(opencl_err_message,"CL_INVALID_KERNEL");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_ARG_INDEX:
            strcpy(opencl_err_message,"CL_INVALID_ARG_INDEX");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_ARG_VALUE:
            strcpy(opencl_err_message,"CL_INVALID_ARG_VALUE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_ARG_SIZE:
            strcpy(opencl_err_message,"CL_INVALID_ARG_SIZE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL_ARGS:
            strcpy(opencl_err_message,"CL_INVALID_KERNEL_ARGS");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_WORK_DIMENSION:
            strcpy(opencl_err_message,"CL_INVALID_WORK_DIMENSION");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_WORK_GROUP_SIZE:
            strcpy(opencl_err_message,"CL_INVALID_WORK_GROUP_SIZE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_WORK_ITEM_SIZE:
            strcpy(opencl_err_message,"CL_INVALID_WORK_ITEM_SIZE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_GLOBAL_OFFSET:
            strcpy(opencl_err_message,"CL_INVALID_GLOBAL_OFFSET");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_EVENT_WAIT_LIST:
            strcpy(opencl_err_message,"CL_INVALID_EVENT_WAIT_LIST");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_EVENT:
            strcpy(opencl_err_message,"CL_INVALID_EVENT");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_OPERATION:
            strcpy(opencl_err_message,"CL_INVALID_OPERATION");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_GL_OBJECT:
            strcpy(opencl_err_message,"CL_INVALID_GL_OBJECT");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_BUFFER_SIZE:
            strcpy(opencl_err_message,"CL_INVALID_BUFFER_SIZE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_MIP_LEVEL:
            strcpy(opencl_err_message,"CL_INVALID_MIP_LEVEL");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_GLOBAL_WORK_SIZE:
            strcpy(opencl_err_message,"CL_INVALID_GLOBAL_WORK_SIZE");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
        default :
            strcpy(opencl_err_message,"Unknown error");
            pycbcopencl_error_stash = PYCBC_RUNTIME_ERROR;
            break;
    }
    
    sprintf(pycbcopencl_error_message, " %s %s", opencl_err_message, 
	    generic_err_message);
    return;
}

void pycbc_opencl_clear_error()
{
    //printf("debug: pycbc_opencl_clear_error\n"); 
    pycbcopencl_error_stash = 0;
}


cl_context_t* new_cl_context_t(unsigned device_id)
{
  cl_context_t* c  = NULL;;
    int opencl_err = 0;

    c = (cl_context_t*) malloc(sizeof(cl_context_t));

    c->device_id = device_id;
    c->set_error = pycbc_opencl_set_error;
    
    // this will update c with all OpenCl context elements
    opencl_err = gpuinsp_InitGPU(c, device_id);
    
    if(opencl_err != CL_SUCCESS)
      c->set_error( opencl_err, "gpuinsp_InitGPU failed");

    return c;
}

void delete_cl_context_t( cl_context_t* p )
{
    free( p );
}

