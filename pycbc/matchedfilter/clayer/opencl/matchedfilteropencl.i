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

// inline definition of a function to wrap:
//%include "../../../datavector/clayer/opencl/datavectoropencl_types.h"
//%include "matchedfilteropencl_types.h"
//%inline %{
//real_vector_t* gen_snr_opencl(real_vector_t* stilde, real_vector_t* htilde)
//{
//    return stilde;
//
//}
//%}

// prototype declaration of a function to wrap (has to be impl. in a c-file)
void gen_snr_opencl(cl_context_t* context,
                    real_vector_single_opencl_t* snr,
                    complex_vector_single_opencl_t* stilde, 
                    complex_vector_single_opencl_t* htilde);
