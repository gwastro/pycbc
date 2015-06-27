#!/usr/bin/python
import pycbc, os, subprocess, ctypes
from mako.template import Template

full_corr = """
    __device__ cufftComplex in_call(void* input, size_t offset, 
                                void* caller_info, void* shared) {                                
        cufftComplex r;
        
        cufftComplex s = ((cufftComplex*) input)[offset];
        cufftComplex h = ((cufftComplex*) %s)[offset];
        
        r.x = h.x * s.x + h.y * s.y;
        r.y = h.x * s.y - h.y * s.x;
        
        return r;
    }
"""

half_zero_corr = """
    __device__ cufftComplex in_call(void* input, size_t offset, 
                                void* caller_info, void* shared) {
        if (offset > %s)                        
            return (cufftComplex){0, 0};
        else{                                   
            cufftComplex r;
            
            cufftComplex s = ((cufftComplex*) input)[offset];
            cufftComplex h = ((cufftComplex*) %s)[offset];
            
            r.x = h.x * s.x + h.y * s.y;
            r.y = h.x * s.y - h.y * s.x;
            
            return r;
        }
        
    }
"""

copy_callback = """
    __device__ cufftComplex correlate(void* input, size_t offset, 
                                void* caller_info, void* shared) {
        return ((cufftComplex*)input)[offset];
    }
"""

copy_out = """
    __device__ void out_call(void *out, size_t offset, cufftComplex element, 
                             void *caller_info, void *shared){
           ((cufftComplex*) out)[offset] = element;
    }

"""


no_out = """
    __device__ void out_call(void *out, size_t offset, cufftComplex element, 
                             void *caller_info, void *shared){
    }

"""

fftsrc = Template("""
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
    
    % if input_callback:
        ${input_callback}
        __device__ cufftCallbackLoadC input_callback = in_call; 
    % endif
    
    % if output_callback:
        ${output_callback}
        __device__ cufftCallbackStoreC output_callback = out_call; 
    % endif
    
    extern "C"  cufftHandle* create_plan(unsigned int size, void* other){
        cufftHandle* plan = new cufftHandle;
        size_t work_size;
        cufftCreate(plan);
        checkCudaErrors(cufftMakePlan1d(*plan, size, CUFFT_C2C, 1, &work_size));
        
        
        % if input_callback:
            cufftCallbackLoadC h_input_callback;
            checkCudaErrors(cudaMemcpyFromSymbol(&h_input_callback, input_callback, 
                                             sizeof(h_input_callback)));         
            checkCudaErrors(cufftXtSetCallback(*plan, (void **) &h_input_callback,
                                          CUFFT_CB_LD_COMPLEX, 
                                          &other));
        % endif
        
        % if output_callback:
            cufftCallbackStoreC h_output_callback;                                         
            checkCudaErrors(cudaMemcpyFromSymbol(&h_output_callback, output_callback, 
                                                 sizeof(h_output_callback)));                                            
            checkCudaErrors(cufftXtSetCallback(*plan, (void **) &h_output_callback,
                                              CUFFT_CB_ST_COMPLEX, 
                                              &other));
        % endif

        return plan;
    }
       
    extern "C" void execute(cufftHandle* plan, cufftComplex* in, cufftComplex* out){   
         checkCudaErrors(cufftExecC2C(*plan, in, out, CUFFT_INVERSE));
    }    
""")

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
 
def get_fn_plan(callback=None, out_callback=None, name='pycbc_cufft'):
    """ Get the IFFT execute and plan functions
    """
    source = fftsrc.render(input_callback=callback, output_callback=out_callback)
    path = compile(source, name)
    lib = ctypes.cdll.LoadLibrary(path)
    fn = lib.execute
    fn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    plan = lib.create_plan
    plan.restype = ctypes.c_void_p
    plan.argyptes = [ctypes.c_uint, ctypes.c_void_p]
    return fn, plan

_plans = {}
# This relies on the fact that we don't change the memory locations
def c2c_correlate_ifft(htilde, stilde, outvec):
    key = 'cnf'
    if key not in _plans:
        fn, pfn = get_fn_plan(callback=full_corr % (int(htilde.data.gpudata)))
        plan = pfn(len(outvec), int(htilde.data.gpudata))
        _plans[key] = (fn, plan, int(htilde.data.gpudata))
    fn, plan, h = _plans[key]
    fn(plan, int(stilde.data.gpudata), int(outvec.data.gpudata))

# This relies on the fact that we don't change the memory locations
def c2c_half_correlate_ifft(htilde, stilde, outvec):
    key = 'cn'
    if key not in _plans:
        fn, pfn = get_fn_plan(callback=half_zero_corr % (len(outvec) / 2, int(htilde.data.gpudata)),
                              out_callback=None)
        plan = pfn(len(outvec), int(htilde.data.gpudata))
        _plans[key] = (fn, plan, int(htilde.data.gpudata))
    fn, plan, h = _plans[key]
    fn(plan, int(stilde.data.gpudata), int(outvec.data.gpudata))

# This relies on the fact that we don't change the memory locations
def c2c_half_correlate_ifft2(htilde, stilde, outvec):
    key = 'cn'
    if key not in _plans:
        fn, pfn = get_fn_plan(callback=half_zero_corr % (len(outvec) / 2, int(htilde.data.gpudata)),
                              out_callback=copy_out)
        plan = pfn(len(outvec), int(htilde.data.gpudata))
        _plans[key] = (fn, plan, int(htilde.data.gpudata))
    fn, plan, h = _plans[key]
    fn(plan, int(stilde.data.gpudata), int(outvec.data.gpudata))
    
    # This relies on the fact that we don't change the memory locations
def c2c_half_correlate_ifft3(htilde, stilde, outvec):
    key = 'cn'
    if key not in _plans:
        fn, pfn = get_fn_plan(callback=half_zero_corr % (len(outvec) / 2, int(htilde.data.gpudata)),
                              out_callback=no_out)
        plan = pfn(len(outvec), int(htilde.data.gpudata))
        _plans[key] = (fn, plan, int(htilde.data.gpudata))
    fn, plan, h = _plans[key]
    fn(plan, int(stilde.data.gpudata), int(outvec.data.gpudata))
    
# This function will only work with a single vector size in the entire program
def c2c_zeros_ifft(corr, outvec):
    key = len(corr)
    if key not in _plans:
        fn, pfn = get_fn_plan(callback=half_zero % (len(corr) / 2))
        plan = pfn(len(corr), 0)
        _plans[key] = (fn, plan)
    fn, plan = _plans[key]
    fn(plan, int(corr.data.gpudata), int(outvec.data.gpudata))
