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

// To add new elements add them here in the typedef meta_data_t, 
// in the CONSTRUCTOR_TEMPLATE 
// and in the TYPE_INTERFACE_TEMPLATE macro in the datavector_types.i file

typedef struct
{
    unsigned long int start;
    double            dx;
    unsigned int      vector_length;
    size_t            element_size_bytes;
    int               generic_new_element_for_testing;
}
meta_data_t;

#define CONSTRUCTOR_TEMPLATE(name,type)\
name* c;\
c = (name*) malloc( sizeof(name) );\
c->meta_data.start = 0;\
c->meta_data.dx = 1;\
c->meta_data.vector_length = length;\
c->meta_data.element_size_bytes = sizeof(type);\
c->meta_data.generic_new_element_for_testing =10;\

#endif /* DATAVECTOR_TYPES_H */
