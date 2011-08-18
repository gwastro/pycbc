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
// datavector meta data type for pycbc

#ifndef DATAVECTOR_TYPES_H
#define DATAVECTOR_TYPES_H

#include <stdlib.h>

// To add new elements add them in the typedef meta_data_t  
// and in the CONSTRUCTOR_TEMPLATE macro
// and in the TYPE_INTERFACE_TEMPLATE macro in the datavector_types.i file

typedef struct
{
    
    // add "epoch" -> segment start time accociated w/ initial time series 
    unsigned long start;    // better rename to x0;
                            // offset of the data from origin series applies 
                            // only for segmenting ????
    double        delta_x;  // Depending on the data either sample intervall
                            // in time domain or frequency domain
    
    unsigned long vector_length;
    size_t        element_size_bytes;
}
meta_data_t;

typedef struct
{
    float re;
    float im;
}
complex_float_t;

typedef struct
{
    double re;
    double im;
}
complex_double_t;


// TODO start needs to be added to constructor

#define CONSTRUCTOR_TEMPLATE(name,type)\
name* c;\
c = (name*) malloc( sizeof(name) );\
c->meta_data.start = 0;\
c->meta_data.delta_x = delta_x;\
c->meta_data.vector_length = length;\
c->meta_data.element_size_bytes = sizeof(type);\

#endif /* DATAVECTOR_TYPES_H */
