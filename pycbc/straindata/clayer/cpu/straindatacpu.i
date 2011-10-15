// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:syntax=c

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING,module="pycbc.straindata.clayer") cpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include <datavectorcpu.h>
#include "straindatacpu_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

// To extend the c-types by methodes they have to be defined here
// but to declare function prototypes as well would raise a 
// "is multiply defined error". That is the reason for splitting 
// the headerfiles

// prototype declaration of a function to wrap (has to be impl. in a c-file)
void* fftw_generate_plan(unsigned long length, real_vector_single_cpu_t* in_tmp,
                     complex_vector_single_cpu_t* out_tmp, char* sign, char* style);

int fftw_transform_segments(void* plan, real_vector_single_cpu_t* in_buf, 
                            unsigned long input_buf_offset,
                            complex_vector_single_cpu_t* out_buf);

int frame_cpp_read_frames(real_vector_double_cpu_t* out_buf, char* channel_name, 
                          unsigned long gps_start_time, unsigned long gps_end_time, 
                          char* cache_url);

