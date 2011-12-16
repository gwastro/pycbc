// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:syntax=c

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING,module="pycbc.waveformgenerator.clayer") waveformgeneratorcpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "datavectorcpu.h"
#include "pycbccpu.h"    
#include "waveformgeneratorcpu.h"
#include "waveformgeneratorcpu_private.h"
    %}

%pythoncode %{
    import types
    import warnings
    %}

%include "waveformgeneratorcpu.h"

%extend waveform_generator_cpu_t {
    
    waveform_generator_cpu_t(cpu_context_t*);
    ~waveform_generator_cpu_t();
    
    //%typemap(check) unsigned some_XYZ_constructor_parameter {
    //    if ($1 < 0) {
    //        SWIG_exception(SWIG_ValueError, "Invalid some_XYZ_constructor_parameter");
    //    }   
    //}
    
    char* __str__() {
        static char a[512];
        snprintf( a, sizeof(a)/sizeof(*a), 
                 "<waveform_generator_cpu_t nothing_to_define_yet %d>", 
                 self->nothing_to_define_yet);
        return a;
    }
}
