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
#include "pycbccpu.h"
#include "datavectorcpu.h"
#include "datavectorcpu_private.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

%include "datavector.h"
%include "datavectorcpu.h"
%include "exception.i"

// add new element properties for CPU vectors here:
%define TYPE_INTERFACE_TEMPLATE(name,type)
name(cpu_context_t *, unsigned long vector_length, double delta_x);
~name();

%typemap(in) complex float {
    $1 = (float)PyComplex_RealAsDouble($input) + (float)PyComplex_ImagAsDouble($input)*I;
}

%typemap(in) complex double {
    $1 = PyComplex_RealAsDouble($input) + PyComplex_ImagAsDouble($input)*I;
}

%typemap(out) complex float {
    $result = PyComplex_FromDoubles((double) crealf($1), (double)cimagf($1));    
}

%typemap(out) complex double {
    $result = PyComplex_FromDoubles(creal($1), cimag($1)); 
}

%typemap(check) unsigned long vector_index {
    if ($1 >= arg1->meta_data.vector_length) {
        SWIG_exception(SWIG_ValueError, "Index for datavector access out of range");
    }
}


char* __str__() {
    static char a[512];
    snprintf( a, sizeof(a)/sizeof(*a), 
             "<name, length %ld, data ptr %p>", 
             self->meta_data.vector_length, self->data );
    return a;
}

unsigned __len__() {
    return self->meta_data.vector_length;
}

type __getitem__(unsigned long vector_index) {
    type* data = (type*) self->data; 
    return (type) data[vector_index];
}

void __setitem__(unsigned long vector_index, type value) {
    type* data = (type*) self->data; 
    data[vector_index] = value;
}

void set_start( unsigned long vector_index ) {
    self->meta_data.start = vector_index;
}

unsigned get_start( void ) {
    return self->meta_data.start;
}

void set_delta_x( double delta_x ) {
    self->meta_data.delta_x = delta_x;
}

double get_delta_x( void ) {
    return self->meta_data.delta_x;
}

%enddef

// exceptions section ----------------------------------------------------------

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
