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

%module(docstring=DOCSTRING) pycbccpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
%{
#include "pycbccpu_types.h"
#include "pycbccpu_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

// To extend the c-types by methodes they have to be defined here
// but to declare function prototypes as well would raise a 
// "is multiply defined error". That is the reason for splitting 
// the headerfiles into _types and _prototypes
// prototype declaration of a function to wrap (has to be impl. in a c-file)

%include "pycbccpu_types.h"
%include "exception.i"

// Generic errorhandling 
%exception {
    char* err_message;
    $action
    if (pycbc_cpu_check_err_occurred()) {
        err_message= pycbc_cpu_get_err_message();
        pycbc_cpu_clear_error();
        SWIG_exception(SWIG_RuntimeError, err_message);
        return NULL;
    }
}

// Function specific errorhandling
%exception cpu_context_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "cpu_context_t allocation fails");
        return NULL;
    }
}
%extend cpu_context_t {

    cpu_context_t(unsigned device_id);
    ~cpu_context_t();

    %typemap(check) unsigned device_id {
        if ($1 < 0) {
            SWIG_exception(SWIG_ValueError, "Invalid device id");
        }   
    }

    char* __str__() {
        static char a[512];
        snprintf( a, sizeof(a)/sizeof(*a), 
             "<cpu_context_t at %p, device_id %d>", 
             self, self->device_id);
        return a;
    }
}