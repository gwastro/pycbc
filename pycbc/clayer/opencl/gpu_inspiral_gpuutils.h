extern "C" cl_int        gpuinsp_checkError(cl_int err, const char* errstr);
extern "C" void          gpuinsp_getErrMessage(char * errMessage, const cl_int err);
extern "C" cl_int        gpuinsp_DestroyGPU(cl_context_t * c);
extern "C" cl_int        gpuinsp_InitGPU(cl_context_t* c, unsigned device_id);

