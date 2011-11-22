cl_int        gpuinsp_checkError(cl_int err, const char* errstr);
void          gpuinsp_getErrMessage(char * errMessage, const cl_int err);
cl_int        gpuinsp_DestroyGPU(cl_context_t * c);
cl_int        gpuinsp_InitGPU(cl_context_t* c, unsigned device_id);
