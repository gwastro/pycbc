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
// datavectorcpu swig file for pycbc

/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING,module="pycbc.datavector.clayer") cpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
%{
#include "datavectorcpu.h"
#include "datavectorcpu_private.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

%include "../datavector.h"
%include "../datavector.i"
%include "datavectorcpu.h"
%include "exception.i"

%exception real_vector_single_cpu_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "real_vector_single_cpu_t allocation fails");
        return NULL;
    }
}
%extend real_vector_single_cpu_t {
    TYPE_INTERFACE_TEMPLATE(real_vector_single_cpu_t,float)
}

%exception real_vector_double_cpu_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "real_vector_double_cpu_t allocation fails");
        return NULL;
    }
}
%extend real_vector_double_cpu_t {
    TYPE_INTERFACE_TEMPLATE(real_vector_double_cpu_t,double)
}

%exception complex_vector_single_cpu_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "complex_vector_single_cpu_t allocation fails");
        return NULL;
    }
}
%extend complex_vector_single_cpu_t {
    TYPE_INTERFACE_TEMPLATE(complex_vector_single_cpu_t,complex float)
}

%exception complex_vector_double_cpu_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "complex_vector_double_cpu_t allocation fails");
        return NULL;
    }
}
%extend complex_vector_double_cpu_t {
    TYPE_INTERFACE_TEMPLATE(complex_vector_double_cpu_t,complex double )
}
