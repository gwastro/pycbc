// Copyright (C) 2011 Karsten Wiesner
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


//
// =============================================================================
//
//                                   Preamble
//
// =============================================================================
//
// pycbc swig file for pycbc

/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) pycbcopencl

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
%{
#include "pycbcopencl.h"
#include "pycbcopencl_private.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

%include "pycbcopencl.h"
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

// Function specific errorhandling
%exception cl_context_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "cl_context_t allocation fails");
        return NULL;
    }
}
%extend cl_context_t {

    cl_context_t(unsigned device_id);
    ~cl_context_t();

    %typemap(check) unsigned device_id {
        if ($1 < 0) {
            SWIG_exception(SWIG_ValueError, "Invalid device id");
        }   
    }

    char* __str__() {
        static char a[512];
        snprintf( a, sizeof(a)/sizeof(*a), 
             "<cl_context_t at %p, device_id %d>", 
             self, self->device_id);
        return a;
    }
}
