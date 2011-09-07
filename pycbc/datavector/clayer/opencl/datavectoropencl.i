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
// datavector swig file for pycbc

/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) datavectoropencl

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
%{
#include "datavectoropencl_types.h"
#include "datavectoropencl_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

// To extend the c-types by methodes they have to be defined here
// but to declare function prototypes as well would raise a 
// "is multiply defined error". That is the reason for splitting 
// the headerfiles into _types and _prototypes
%include "../datavector_types.h"
// opencl datavectors have two pointers so 
// we have to implement the swigging different than -> %include "../datavector_types.i"
%include "datavectoropencl_types.h"
%include "exception.i"

// swig properties for datavector elements

%{
#include <complex.h>
%}

// add new element properties here:
%define TYPE_INTERFACE_TEMPLATE(name,type)
name(unsigned long vector_length, double delta_x);
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
//    snprintf( a, sizeof(a)/sizeof(*a), 
//             "<name, length %ld, real_data ptr %p, imag_data ptr %p>", 
//             self->meta_data.vector_length, self->real_data, self->imag_data);
    
// obvoiusly we have the ptoblem here of real_vecs have only a *data ptr
// and complex ones have the two real and imag data ptrs    

    return a;
}

unsigned __len__() {
    return self->meta_data.vector_length;
}


// ToDo
// As well it might make trouble to NOT use the C99 complex type here! 
// BTW.: Also check all the complex to python complex stuff at the other
// <arch> datavectors !!!!!!   

///
/// ToDo find out to get/set complex data ////////////////////////////////
///

// same problem applies here:
// obvoiusly we have the ptoblem here of real_vecs have only a *data ptr
// and complex ones have the two real and imag data ptrs    


/*type __getitem__(unsigned long vector_index) {
    type* data = (type*) self->real_data; 
    return (type) real_data[vector_index];
}

void __setitem__(unsigned long vector_index, type value) {
    type* data = (type*) self->real_data; 
    real_data[vector_index] = value;
}
 
*/

//////////////////////////////////////////////////////////////////////////



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
    TYPE_INTERFACE_TEMPLATE(real_vector_single_opencl_t,float)
}

%exception real_vector_double_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "real_vector_double_opencl_t allocation fails");
        return NULL;
    }
}
%extend real_vector_double_opencl_t {
    TYPE_INTERFACE_TEMPLATE(real_vector_double_opencl_t,double)
}

%exception complex_vector_single_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "complex_vector_single_opencl_t allocation fails");
        return NULL;
    }
}
%extend complex_vector_single_opencl_t {
    TYPE_INTERFACE_TEMPLATE(complex_vector_single_opencl_t,complex_float_t)
}

%exception complex_vector_double_opencl_t {
    $action
    if (!result) {
        SWIG_exception(SWIG_MemoryError, "complex_vector_double_opencl_t allocation fails");
        return NULL;
    }
}
%extend complex_vector_double_opencl_t {
    TYPE_INTERFACE_TEMPLATE(complex_vector_double_opencl_t,complex_double_opencl_t)
}

