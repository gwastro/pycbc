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
// datavector constructors and destructors implementation for pycbc

#include <stdio.h>
#include "datavectorcpu_types.h"
#include "datavectorcpu_prototypes.h"

real_vector_single_t* new_real_vector_single_t(int length)
{

    real_vector_single_t* c;
    c = (real_vector_single_t*) malloc( sizeof(real_vector_single_t) );

//    snprintf( c->meta_data.module_name, sizeof(c->meta_data.module_name)/sizeof(*c->meta_data.module_name), 
//             "datavectorcpu");
    
    c->meta_data.start = 0;
    c->meta_data.dx = 1;
    c->meta_data.vector_length = length;
    
    c->meta_data.element_size_bytes = sizeof(float); // single precision specified here
    c->data = calloc( c->meta_data.vector_length , c->meta_data.element_size_bytes );

    printf("created real_vector_single_t at %p\n", c );
    return c;
}

void delete_real_vector_single_t( real_vector_single_t* p )
{
    printf("deleting real_vector_single_t at %p\n", p );
    free( p->data );
    free( p );
}

real_vector_double_t* new_real_vector_double_t(int length)
{
    
    real_vector_double_t* c;
    c = (real_vector_double_t*) malloc( sizeof(real_vector_double_t) );
    
    c->meta_data.start = 0;
    c->meta_data.dx = 1;
    c->meta_data.vector_length = length;
    
    c->meta_data.element_size_bytes = sizeof(double); // double precision specified here
    c->data = calloc( c->meta_data.vector_length , c->meta_data.element_size_bytes );
    
    printf("created real_vector_double_t at %p\n", c );
    return c;
}

void delete_real_vector_double_t( real_vector_double_t* p )
{
    printf("deleting real_vector_double_t at %p\n", p );
    free( p->data );
    free( p );
}

complex_vector_single_t* new_complex_vector_single_t(int length)
{
    
    complex_vector_single_t* c;
    c = (complex_vector_single_t*) malloc( sizeof(complex_vector_single_t) );
    
    c->meta_data.start = 0;
    c->meta_data.dx = 1;
    c->meta_data.vector_length = length;
    
    c->meta_data.element_size_bytes = sizeof(float);   // single precision specified here
    c->data = calloc( 2 * c->meta_data.vector_length , // complex => 2 * length
                      c->meta_data.element_size_bytes );
    
    printf("created complex_vector_single_t at %p\n", c );
    return c;
}

void delete_complex_vector_single_t( complex_vector_single_t* p )
{
    printf("deleting complex_vector_single_t at %p\n", p );
    free( p->data );
    free( p );
}

complex_vector_double_t* new_complex_vector_double_t(int length)
{
    
    complex_vector_double_t* c;
    c = (complex_vector_double_t*) malloc( sizeof(complex_vector_double_t) );
    
    c->meta_data.start = 0;
    c->meta_data.dx = 1;
    c->meta_data.vector_length = length;
    
    c->meta_data.element_size_bytes = sizeof(double);  // double precision specified here
    c->data = calloc( 2 * c->meta_data.vector_length , // complex => 2 * length
                      c->meta_data.element_size_bytes );
    
    printf("created complex_vector_double_t at %p\n", c );
    return c;
}

void delete_complex_vector_double_t( complex_vector_double_t* p )
{
    printf("deleting complex_vector_double_t at %p\n", p );
    free( p->data );
    free( p );
}

