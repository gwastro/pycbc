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
// swig properties for datavector elements


// add new element properties here:
%define TYPE_INTERFACE_TEMPLATE(name,type)
name(unsigned long vector_length);
~name();

%typemap(in) complex_float_t {
    $1.re = (float)PyComplex_RealAsDouble($input);
    $1.im = (float)PyComplex_ImagAsDouble($input);
}

%typemap(in) complex_double_t {
    $1.re = PyComplex_RealAsDouble($input);
    $1.im = PyComplex_ImagAsDouble($input);
}

%typemap(out) complex_float_t {
    $result = PyComplex_FromDoubles((double)$1.re, (double)$1.im);    
}

%typemap(out) complex_double_t {
    $result = PyComplex_FromDoubles($1.re, $1.im); 
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

    /* Even though the exception appeares in python it 
       can't be catched there by "except". Because the python interpreter
       get's confused if the funcion returns a none-NULL 
    try {
        if (i >= self->meta_data.vector_length)
            throw(RangeError);
        type* data = (type*) self->data; 
        data[i] = value;
    } catch (RangeError) {
        PyErr_SetString(PyExc_IndexError,"Index for datavector access out of range");
    } finally 
        PyErr_SetString(PyExc_MemoryError,"Unknown exception while datavector access");
    */
}

void set_start( unsigned long vector_index ) {
    self->meta_data.start = vector_index;
}

unsigned get_start( void ) {
    return self->meta_data.start;
}

void set_dx( double dx ) {            // ToDo better rename dx to delta_t
    self->meta_data.dx = dx;
}

double get_dx( void ) {               // "  "
    return self->meta_data.dx;
}

void set_generic_new_element_for_testing( int value ) {
    self->meta_data.generic_new_element_for_testing = value;
}

int get_generic_new_element_for_testing( void ) {
    return self->meta_data.generic_new_element_for_testing;
}

%enddef

