#include <cuda_runtime.h>

#ifndef PYCBC_CUDA_EXCEPT_H
#define PYCBC_CUDA_EXCEPT_H

void pycbc_throw_cuda_exception(cudaError_t cuda_error_t, ...);

#endif
