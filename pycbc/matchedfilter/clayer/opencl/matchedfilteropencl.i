/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) matchedfilteropencl

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "datavectoropencl.h"
#include "pycbcopencl.h"
#include "pycbcopencl_private.h"
#include "matchedfilteropencl.h"    
#include "matchedfilteropencl_private.h"
%}

%pythoncode %{
    import types
    import warnings
    %}


%include "matchedfilteropencl.h"
%include "exception.i"


%exception {
    char* error_message;
    int   error_id;
    $action
    if ( (error_id = pycbc_opencl_check_error()) )
    {
	    error_message = pycbc_opencl_get_error_message();
	    pycbc_opencl_clear_error();
	    switch (error_id)
        {
            case PYCBC_ATTRIBUTE_ERROR:
                PyErr_SetString(PyExc_AttributeError, error_message);
                break;
            case PYCBC_EOF_ERROR:
                PyErr_SetString(PyExc_EOFError, error_message);
                break;
            case PYCBC_IO_ERROR:
                PyErr_SetString(PyExc_IOError, error_message);
                break;
            case PYCBC_INDEX_ERROR:
                PyErr_SetString(PyExc_IndexError, error_message);
                break;
            case PYCBC_TYPE_ERROR:
                PyErr_SetString(PyExc_TypeError, error_message);
                break;
            case PYCBC_VALUE_ERROR:
                PyErr_SetString(PyExc_ValueError, error_message);
                break;
            case PYCBC_MEMORY_ERROR:
                PyErr_SetString(PyExc_MemoryError, error_message);
                break;
            case PYCBC_NAME_ERROR:
                PyErr_SetString(PyExc_NameError, error_message);
                break;
            case PYCBC_OVERFLOW_ERROR:
                PyErr_SetString(PyExc_OverflowError, error_message);
                break;
            case PYCBC_ZERO_DIVISION_ERROR:
                PyErr_SetString(PyExc_ZeroDivisionError, error_message);
                break;
            case PYCBC_RUNTIME_ERROR:
            default:
                PyErr_SetString(PyExc_RuntimeError, error_message);
        }
        return NULL;
    }
}

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
                 "<matched_filter_opencl_t at %p, prog %p, cl_stilde %p, cl_htilde %p, gpu_snr_prod %p>", 
		  self, self->program, self->cl_stilde, self->cl_htilde, self->gpu_snr_product);
        return a;
    }
}

// prototype declaration of a function to wrap (has to be impl. in a c-file)
void gen_snr_opencl(cl_context_t* context,
                    matched_filter_opencl_t * matchedfilter,
                    real_vector_single_opencl_t* snr,
                    complex_vector_single_opencl_t* stilde, 
                    complex_vector_single_opencl_t* htilde);
