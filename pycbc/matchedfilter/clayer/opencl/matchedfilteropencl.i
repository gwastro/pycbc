/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) matchedfilteropencl

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "../../../datavector/clayer/opencl/datavectoropencl_types.h"
#include "../../../clayer/opencl/pycbcopencl_types.h"
#include "matchedfilteropencl_types.h"    
#include "matchedfilteropencl_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

// To extend the c-types by methodes they have to be defined here
// but to declare function prototypes as well would raise a 
// "is multiply defined error". That is the reason for splitting 
// the headerfiles

%include "matchedfilteropencl_types.h"
%include "exception.i"

%exception matched_filter_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "matched_filter_opencl_t allocation fails");
        return NULL;
    }
}
%extend matched_filter_opencl_t {
    
    matched_filter_opencl_t(cl_context_t*);
    ~matched_filter_opencl_t();
    
    //%typemap(check) unsigned some_MF_constructor_parameter {
    //    if ($1 < 0) {
    //        SWIG_exception(SWIG_ValueError, "Invalid some_MF_constructor_parameter");
    //    }   
    //}
    
    char* __str__() {
        static char a[512];
        snprintf( a, sizeof(a)/sizeof(*a), 
                 "<matched_filter_opencl_t at %p, some_MF_parameter %d>", 
                 self, 5); // self->device_id);
        return a;
    }
}

// prototype declaration of a function to wrap (has to be impl. in a c-file)
void gen_snr_opencl(cl_context_t* context,
                    matched_filter_opencl_t * matchedfilter,
                    real_vector_single_opencl_t* snr,
                    complex_vector_single_opencl_t* stilde, 
                    complex_vector_single_opencl_t* htilde);
