// Copyright (C) 2011 Karsten Wiesner, Josh Willis
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
