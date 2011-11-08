#include <stdio.h>
#include <string.h>
#include <pycbc/clayer/except.h>
#include <cuda_runtime.h>
#include <stdarg.h>

void pycbc_throw_cuda_exception(cudaError_t cuda_error_t, ...)
{    
    char error_message[EXCEPTION_MESSAGE_SIZE+1];
    char generic_error_message[EXCEPTION_MESSAGE_SIZE+1]; 
    unsigned int error_status;
    
    va_list argp;
    va_start(argp, cuda_error_t);
    char* fmt=va_arg(argp,char*);
    if (fmt != (char*)0)
      vsprintf(generic_error_message, fmt, argp);
    va_end(argp); 
    
    //Map cuda error codes to pycbc error codes
    switch (cuda_error_t) {
      case cudaSuccess:                     error_status = PYCBC_NO_ERROR; break;
      case cudaErrorMissingConfiguration:   error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorMemoryAllocation:       error_status = PYCBC_MEMORY_ERROR; break;
      case cudaErrorInitializationError:    error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorLaunchFailure:          error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorPriorLaunchFailure:     error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorLaunchTimeout:          error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorLaunchOutOfResources:   error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidDeviceFunction:  error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidConfiguration:   error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidDevice :         error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidValue :          error_status = PYCBC_VALUE_ERROR; break;
      case cudaErrorInvalidPitchValue:      error_status = PYCBC_VALUE_ERROR; break;
      case cudaErrorInvalidSymbol :         error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorMapBufferObjectFailed : error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorUnmapBufferObjectFailed : error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidHostPointer:     error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidDevicePointer:   error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidTexture 	:       error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidTextureBinding : error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidChannelDescriptor: error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidMemcpyDirection: error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorAddressOfConstant :     error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorTextureFetchFailed :    error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorTextureNotBound 	:     error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorSynchronizationError:   error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInvalidFilterSetting :  error_status = PYCBC_VALUE_ERROR; break;
      case cudaErrorInvalidNormSetting :    error_status = PYCBC_VALUE_ERROR; break;
      case cudaErrorMixedDeviceExecution 	: error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorCudartUnloading 	:     error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorUnknown 	 :            error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorNotYetImplemented :     error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorMemoryValueTooLarge :   error_status = PYCBC_VALUE_ERROR; break;
      case cudaErrorInvalidResourceHandle:  error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorNotReady 	:             error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorInsufficientDriver 	 :  error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorSetOnActiveProcess 	:   error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorNoDevice 	 :            error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorStartupFailure 	:       error_status = PYCBC_RUNTIME_ERROR; break;
      case cudaErrorApiFailureBase 	:       error_status = PYCBC_RUNTIME_ERROR; break;
    }
    
    char* cuda_message=cudaGetErrorString(cuda_error_t);

    //append the generic error message to the cuda specific one
    strncat(error_message,cuda_message,EXCEPTION_MESSAGE_SIZE);
    strncat(error_message,":",EXCEPTION_MESSAGE_SIZE-strlen(error_message));
    strncat(error_message,generic_error_message,EXCEPTION_MESSAGE_SIZE-strlen(error_message));
    
    pycbc_throw_exception_bare(error_status,error_message);
    return;
}
