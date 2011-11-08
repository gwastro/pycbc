#include <stdio.h>
#include <string.h>
#include <pycbc/clayer/except.h>
#include <gpu_inspiral_gpuutils.h>
#include <stdarg.h>

void pycbc_throw_opencl_exception(int opencl_error_id, ...)
{    
    char error_message[EXCEPTION_MESSAGE_SIZE+1];
    char generic_error_message[EXCEPTION_MESSAGE_SIZE+1]; 
    unsigned int error_status;
    
    va_list argp;
    va_start(argp, opencl_error_id);
    char* fmt=va_arg(argp,char*);
    if (fmt != (char*)0)
      vsprintf(generic_error_message, fmt, argp);
    va_end(argp);   
    
    //Map the opencl error codes to pycbc error codes
    switch (opencl_error_id) {
        case CL_SUCCESS:
            strncpy(error_message,"CL_SUCCESS",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_NO_ERROR;
            break;
        case CL_DEVICE_NOT_FOUND:
            strncpy(error_message,"CL_DEVICE_NOT_FOUND",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_DEVICE_NOT_AVAILABLE:
            strncpy(error_message,"CL_DEVICE_NOT_AVAILABLE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_COMPILER_NOT_AVAILABLE:
            strncpy(error_message,"CL_COMPILER_NOT_AVAILABLE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            strncpy(error_message,"CL_MEM_OBJECT_ALLOCATION_FAILURE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_MEMORY_ERROR;
            break;
        case CL_OUT_OF_RESOURCES:
            strncpy(error_message,"CL_OUT_OF_RESOURCES",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_OUT_OF_HOST_MEMORY:
            strncpy(error_message,"CL_OUT_OF_HOST_MEMORY",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_MEMORY_ERROR;
            break;
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            strncpy(error_message,"CL_PROFILING_INFO_NOT_AVAILABLE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_MEM_COPY_OVERLAP:
            strncpy(error_message,"CL_MEM_COPY_OVERLAP",EXCEPTION_MESSAGE_SIZE);
            break;
        case CL_IMAGE_FORMAT_MISMATCH:
            strncpy(error_message,"CL_IMAGE_FORMAT_MISMATCH",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            strncpy(error_message,"CL_IMAGE_FORMAT_NOT_SUPPORTED",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_BUILD_PROGRAM_FAILURE:
            strncpy(error_message,"CL_BUILD_PROGRAM_FAILURE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_MAP_FAILURE:
            strncpy(error_message,"CL_MAP_FAILURE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
            /* These are opencl 1.1 features not yet defined
             case CL_MISALIGNED_SUB_BUFFER_OFFSET:
             strncpy(error_message,"CL_MISALIGNED_SUB_BUFFER_OFFSET");
             error_status = PYCBC_RUNTIME_ERROR;
             break;
             case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
             strncpy(error_message,"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST");
             error_status = PYCBC_RUNTIME_ERROR;
             break;
             */
        case CL_INVALID_VALUE:
            strncpy(error_message,"CL_INVALID_VALUE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_DEVICE_TYPE:
            strncpy(error_message,"CL_INVALID_DEVICE_TYPE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_PLATFORM:
            strncpy(error_message,"CL_INVALID_PLATFORM",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_DEVICE:
            strncpy(error_message,"CL_INVALID_DEVICE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_CONTEXT:
            strncpy(error_message,"CL_INVALID_CONTEXT",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_QUEUE_PROPERTIES:
            strncpy(error_message,"CL_INVALID_QUEUE_PROPERTIES",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_COMMAND_QUEUE:
            strncpy(error_message,"CL_INVALID_COMMAND_QUEUE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_HOST_PTR:
            strncpy(error_message,"CL_INVALID_HOST_PTR",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_MEMORY_ERROR;
            break;
        case CL_INVALID_MEM_OBJECT:
            strncpy(error_message,"CL_INVALID_MEM_OBJECT",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_MEMORY_ERROR;
            break;
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            strncpy(error_message,"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_IMAGE_SIZE:
            strncpy(error_message,"CL_INVALID_IMAGE_SIZE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_SAMPLER:
            strncpy(error_message,"CL_INVALID_SAMPLER",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_BINARY:
            strncpy(error_message,"CL_INVALID_BINARY",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_BUILD_OPTIONS:
            strncpy(error_message,"CL_INVALID_BUILD_OPTIONS",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_PROGRAM:
            strncpy(error_message,"CL_INVALID_PROGRAM",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_PROGRAM_EXECUTABLE:
            strncpy(error_message,"CL_INVALID_PROGRAM_EXECUTABLE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL_NAME:
            strncpy(error_message,"CL_INVALID_KERNEL_NAME",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL_DEFINITION:
            strncpy(error_message,"CL_INVALID_KERNEL_DEFINITION",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL:
            strncpy(error_message,"CL_INVALID_KERNEL",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_ARG_INDEX:
            strncpy(error_message,"CL_INVALID_ARG_INDEX",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_ARG_VALUE:
            strncpy(error_message,"CL_INVALID_ARG_VALUE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_ARG_SIZE:
            strncpy(error_message,"CL_INVALID_ARG_SIZE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_KERNEL_ARGS:
            strncpy(error_message,"CL_INVALID_KERNEL_ARGS",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_WORK_DIMENSION:
            strncpy(error_message,"CL_INVALID_WORK_DIMENSION",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_WORK_GROUP_SIZE:
            strncpy(error_message,"CL_INVALID_WORK_GROUP_SIZE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_WORK_ITEM_SIZE:
            strncpy(error_message,"CL_INVALID_WORK_ITEM_SIZE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_GLOBAL_OFFSET:
            strncpy(error_message,"CL_INVALID_GLOBAL_OFFSET",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_EVENT_WAIT_LIST:
            strncpy(error_message,"CL_INVALID_EVENT_WAIT_LIST",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_EVENT:
            strncpy(error_message,"CL_INVALID_EVENT",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_OPERATION:
            strncpy(error_message,"CL_INVALID_OPERATION",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_GL_OBJECT:
            strncpy(error_message,"CL_INVALID_GL_OBJECT",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_BUFFER_SIZE:
            strncpy(error_message,"CL_INVALID_BUFFER_SIZE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_MIP_LEVEL:
            strncpy(error_message,"CL_INVALID_MIP_LEVEL",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        case CL_INVALID_GLOBAL_WORK_SIZE:
            strncpy(error_message,"CL_INVALID_GLOBAL_WORK_SIZE",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
        default :
            strncpy(error_message,"Unknown error",EXCEPTION_MESSAGE_SIZE);
            error_status = PYCBC_RUNTIME_ERROR;
            break;
    }

    //append the generic error message to the opencl specific one
    strncat(error_message,":",EXCEPTION_MESSAGE_SIZE-strlen(error_message));
    strncat(error_message,generic_error_message,EXCEPTION_MESSAGE_SIZE-strlen(error_message));
    
    throw_exception_bare(error_status,error_message);
    return;
}
