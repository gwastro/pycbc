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
// datavectoropencl swig file for pycbc

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) datavectoropencl

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
%{
#include "pycbcopencl.h"
#include "datavectoropencl.h"
#include "datavectoropencl_private.h"
%}

%pythoncode %{
    import types
    import warnings
    %}


%include "datavector.h"
%include "datavectoropencl.h"
%include "exception.i"

// swig properties for datavector elements

%{
#include <complex.h>
%}

// add new element properties for OpenCl real vectors here:
%define TYPE_INTERFACE_TEMPLATE_REAL(name)
name(cl_context_t *, unsigned long vector_length, double delta_x);
~name();

char* __str__() {
    static char a[512];
    snprintf( a, sizeof(a)/sizeof(*a), 
             "<name, length %ld, data ptr %p", 
             self->meta_data.vector_length, self->data);

    return a;
}

unsigned __len__() {
    return self->meta_data.vector_length;
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

// add new element properties for OpenCl complex vectors here:
%define TYPE_INTERFACE_TEMPLATE_CMPLX(name)
name(cl_context_t*, unsigned long vector_length, double delta_x);
~name();

char* __str__() {
    static char a[512];
    snprintf( a, sizeof(a)/sizeof(*a), 
             "<name, length %ld, real_data ptr %p, imag_data ptr %p>", 
             self->meta_data.vector_length, self->real_data, self->imag_data);
    
    return a;
}

unsigned __len__() {
    return self->meta_data.vector_length;
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

%exception real_vector_single_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "real_vector_single_opencl_t allocation fails");
        return NULL;
    }
}
%extend real_vector_single_opencl_t {
    TYPE_INTERFACE_TEMPLATE_REAL(real_vector_single_opencl_t)
}

%exception real_vector_double_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "real_vector_double_opencl_t allocation fails");
        return NULL;
    }
}
%extend real_vector_double_opencl_t {
    TYPE_INTERFACE_TEMPLATE_REAL(real_vector_double_opencl_t)
}

%exception complex_vector_single_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "complex_vector_single_opencl_t allocation fails");
        return NULL;
    }
}
%extend complex_vector_single_opencl_t {
    TYPE_INTERFACE_TEMPLATE_CMPLX(complex_vector_single_opencl_t)

}

%exception complex_vector_double_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "complex_vector_double_opencl_t allocation fails");
        return NULL;
    }
}
%extend complex_vector_double_opencl_t {
    TYPE_INTERFACE_TEMPLATE_CMPLX(complex_vector_double_opencl_t)
}

// transfer functions ----------------------------------------------------------
void transfer_real_vector_single_from_cpu(cl_context_t *,
                                          real_vector_single_opencl_t,
                                          real_vector_single_cpu_t);

void transfer_real_vector_single_to_cpu( cl_context_t *,
                                         real_vector_single_cpu_t,
                                         real_vector_single_opencl_t);

void transfer_real_vector_double_from_cpu(cl_context_t *,
                                          real_vector_double_opencl_t,
                                          real_vector_double_cpu_t);

void transfer_real_vector_double_to_cpu( cl_context_t *,
                                         real_vector_double_cpu_t,
                                         real_vector_double_opencl_t);

void transfer_complex_vector_single_from_cpu(cl_context_t*,
                                             complex_vector_single_opencl_t,
                                             complex_vector_single_cpu_t);

void transfer_complex_vector_single_to_cpu( cl_context_t*,
                                            complex_vector_single_cpu_t,
                                            complex_vector_single_opencl_t);

void transfer_complex_vector_double_from_cpu(cl_context_t*,
                                             complex_vector_double_opencl_t,
                                             complex_vector_double_cpu_t);

void transfer_complex_vector_double_to_cpu( cl_context_t*,
                                            complex_vector_double_cpu_t,
                                            complex_vector_double_opencl_t);
