#!/usr/bin/python
import pycbc, os, subprocess, ctypes

preamble = """
    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <math.h>
    #include <cuda_runtime.h>
    #include <cufft.h>
    #include <cufftXt.h>

    #define checkCudaErrors(val)  __checkCudaErrors__ ( (val), #val, __FILE__, __LINE__ )

    template <typename T>
    inline void __checkCudaErrors__(T code, const char *func, const char *file, int line) 
    {
        if (code) {
            fprintf(stderr, "CUDA error at %s:%d code=%d \\"%s\\" \\n",
                    file, line, (unsigned int)code, func);
            cudaDeviceReset();
            exit(EXIT_FAILURE);
        }
    }
"""

zero_callback = """
    __device__ cufftComplex correlate(void* input, size_t offset, 
                                void* caller_info, void* shared) {
        return (cufftComplex){0, 0};
    }
"""

noop_callback = """
    __device__ cufftComplex correlate(void* input, size_t offset, 
                                void* caller_info, void* shared) {
        return ((cufftComplex*)input)[offset];
    }
"""

fftsrc = """
    __device__ cufftCallbackLoadC input_callback = correlate; 
    
    extern "C"  cufftHandle* create_plan(unsigned int size){
        cufftHandle* plan = new cufftHandle;
        size_t work_size;
        cufftCreate(plan);
        checkCudaErrors(cufftMakePlan1d(*plan, size, CUFFT_C2C, 1, &work_size));
        
        cufftCallbackLoadC h_input_callback;
        checkCudaErrors(cudaMemcpyFromSymbol(&h_input_callback, input_callback, 
                                             sizeof(h_input_callback)));       
        checkCudaErrors(cufftXtSetCallback(*plan, (void **)&h_input_callback,
                                          CUFFT_CB_LD_COMPLEX, 0));
        return plan;
    }
       
    extern "C" void execute(cufftHandle* plan, cufftComplex* in, cufftComplex* out){   
         checkCudaErrors(cufftExecC2C(*plan, in, out, CUFFT_INVERSE));
    }    
"""

def compile(source, name):
    """ Compile the string source code into a shared object linked against
    the static version of cufft for callback support.
    """
    cache = os.path.join(pycbc._cache_dir_path, name)
    hash_file = cache + ".hash"
    lib_file = cache + ".so"
    obj_file = cache + ".o"
    
    try:
        if int(open(hash_file, "r").read()) == hash(source):
            return lib_file
        raise ValueError
    except:
        pass

    src_file = cache + ".cu"
    fsrc = open(src_file, "w")
    fsrc.write(source)
    fsrc.close()
    
    cmd = ["nvcc", "-ccbin", "g++", "-dc", "-m64", 
           "--compiler-options", "'-fPIC'", 
           "-o", obj_file, 
           "-c", src_file] 
    print " ".join(cmd)
    subprocess.check_call(cmd)
    
    cmd = ["nvcc", "-shared", "-ccbin", "g++", "-m64", 
       "-o", lib_file, obj_file, "-lcufft_static", "-lculibos"]  
    print " ".join(cmd)
    subprocess.check_call(cmd)
    
    hash_file = cache + ".hash"
    fhash = open(hash_file, "w")
    fhash.write(str(hash(source)))
    return lib_file
 
def get_fn_plan():
    """ Get the IFFT execute and plan functions
    """
    source = preamble + noop_callback + fftsrc
    path = compile(source, "pycbc_cufft")
    lib = ctypes.cdll.LoadLibrary(path)
    fn = lib.execute
    fn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    plan = lib.create_plan
    plan.restype = ctypes.c_void_p
    plan.argyptes = [ctypes.c_uint]
    return fn, plan
 
_plans = {}
def c2c_ifft(invec, outvec):
    if len(invec) not in _plans:
        fn, pfn = get_fn_plan()
        plan = pfn(len(invec))
        _plans[len(invec)] = (fn, plan)
    fn, plan = _plans[len(invec)]
    fn(plan, int(invec.data.gpudata), int(outvec.data.gpudata))
