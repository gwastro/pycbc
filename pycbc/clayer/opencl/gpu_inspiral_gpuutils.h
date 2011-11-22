#include <CL/opencl.h>

gpuinsp_checkError(cl_int err, const char* errstr);
gpuinsp_getErrMessage(char * errMessage, const cl_int err);
gpuinsp_DestroyGPU(cl_context_t * c);
gpuinsp_InitGPU(cl_context_t* c, unsigned device_id);
