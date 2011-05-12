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

%module(docstring=DOCSTRING) datavectorcpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "datavectorcpu_types.h"
#include "datavectorcpu_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

// To extend the c-types by methodes they have to be defined here
// but to declare function prototypes as well would raise a 
// "is multiply defined error". That is the reason for splitting 
// the headerfiles
%include "datavectorcpu_types.h"

%extend real_vector_single_t {
    real_vector_single_t(int vector_length);
    ~real_vector_single_t();
    
    char *__str__() {
        static char a[1024];
        snprintf( a, sizeof(a)/sizeof(*a), 
                     "<real_vector_single_t, in cpu memory, length %d, data ptr %p>", 
                     self->meta_data.vector_length, self->data );
        return a;
    }
    
    int __len__() {
        return self->meta_data.vector_length;
    }
    
    double __getitem__(int i) {
        float* data = (float*) self->data; 
        return (float) data[i];
    }

    void __setitem__(int i, double value) {
        float* data = (float*) self->data; 
        data[i] = value;
    }
    
    void set_start( unsigned long int start ) {
        self->meta_data.start = start;
    }
    
    double get_start( void ) {
        return self->meta_data.start;
    }
    
    void set_dx( double dx ) {
        self->meta_data.dx = dx;
    }
    
    double get_dx( void ) {
        return self->meta_data.dx;
    }
    
    void reinitialize (float init_value) {
        int i;
        float* data = (float*) self->data;
        int length = self->meta_data.vector_length;
        for (i=0; i <= length; i++) {
            *data++ = init_value;
        }
    }    
}

%extend real_vector_double_t {
    real_vector_double_t(int vector_length);
    ~real_vector_double_t();
    
    char *__str__() {
        static char a[1024];
        snprintf( a, sizeof(a)/sizeof(*a), 
                 "<real_vector_double_t, in cpu memory, length %d, data ptr %p>", 
                 self->meta_data.vector_length, self->data );
        return a;
    }
    
    int __len__() {
        return self->meta_data.vector_length;
    }
    
    double __getitem__(int i) {
        double* data = (double*) self->data; 
        return (double) data[i];
    }
    
    void __setitem__(int i, double value) {
        double* data = (double*) self->data; 
        data[i] = value;
    }
    
    void set_start( unsigned long int start ) {
        self->meta_data.start = start;
    }
    
    double get_start( void ) {
        return self->meta_data.start;
    }
    
    void set_dx( double dx ) {
        self->meta_data.dx = dx;
    }
    
    double get_dx( void ) {
        return self->meta_data.dx;
    }
    
    void reinitialize (double init_value) {
        int i;
        double* data = (double*) self->data;
        int length = self->meta_data.vector_length;
        for (i=0; i <= length; i++) {
            *data++ = init_value;
        }
    }    
}

%extend complex_vector_single_t {
    complex_vector_single_t(int vector_length);
    ~complex_vector_single_t();
    
    char *__str__() {
        static char a[1024];
        snprintf( a, sizeof(a)/sizeof(*a), 
                 "<complex_vector_single_t, in cpu memory, length %d, data ptr %p>", 
                 self->meta_data.vector_length, self->data );
        return a;
    }
    
    int __len__() {
        return self->meta_data.vector_length;
    }
    
    double __getitem__(int i) {
        float* data = (float*) self->data; 
        return (float) data[i];
    }
    
    void __setitem__(int i, double value) {
        float* data = (float*) self->data; 
        data[i] = value;
    }
    
    void set_start( unsigned long int start ) {
        self->meta_data.start = start;
    }
    
    double get_start( void ) {
        return self->meta_data.start;
    }
    
    void set_dx( double dx ) {
        self->meta_data.dx = dx;
    }
    
    double get_dx( void ) {
        return self->meta_data.dx;
    }
    
    void reinitialize (float init_value) {
        int i;
        float* data = (float*) self->data;
        int length = self->meta_data.vector_length;
        for (i=0; i <= length; i++) {
            *data++ = init_value;
        }
    }    
}

%extend complex_vector_double_t {
    complex_vector_double_t(int vector_length);
    ~complex_vector_double_t();
    
    char *__str__() {
        static char a[1024];
        snprintf( a, sizeof(a)/sizeof(*a), 
                 "<complex_vector_double_t, in cpu memory, length %d, data ptr %p>", 
                 self->meta_data.vector_length, self->data );
        return a;
    }
    
    int __len__() {
        return self->meta_data.vector_length;
    }
    
    double __getitem__(int i) {
        double* data = (double*) self->data; 
        return (double) data[i];
    }
    
    void __setitem__(int i, double value) {
        double* data = (double*) self->data; 
        data[i] = value;
    }
    
    void set_start( unsigned long int start ) {
        self->meta_data.start = start;
    }
    
    double get_start( void ) {
        return self->meta_data.start;
    }
    
    void set_dx( double dx ) {
        self->meta_data.dx = dx;
    }
    
    double get_dx( void ) {
        return self->meta_data.dx;
    }
    
    void reinitialize (double init_value) {
        int i;
        double* data = (double*) self->data;
        int length = self->meta_data.vector_length;
        for (i=0; i <= length; i++) {
            *data++ = init_value;
        }
    }    
}



