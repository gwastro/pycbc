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
//opencl
// =============================================================================
//
// datavector constructors and destructors implementation for pycbc

#include <stdio.h>
#include "datavectoropencl_types.h"
#include "datavectoropencl_prototypes.h"

real_vector_single_t* new_real_vector_single_t(unsigned length)
{
    CONSTRUCTOR_TEMPLATE(real_vector_single_t, float)
    c->data = calloc( c->meta_data.vector_length , 
                      c->meta_data.element_size_bytes );
    printf("created real_vector_single_t at %p\n", c );
    return c;
}

void delete_real_vector_single_t( real_vector_single_t* p )
{
    printf("deleting real_vector_single_t at %p\n", p );
    free( p->data );
    free( p );
}

real_vector_double_t* new_real_vector_double_t(unsigned length)
{
    
    CONSTRUCTOR_TEMPLATE(real_vector_double_t, double)
    c->data = calloc( c->meta_data.vector_length , 
                      c->meta_data.element_size_bytes );
    printf("created real_vector_double_t at %p\n", c );
    return c;
}

void delete_real_vector_double_t( real_vector_double_t* p )
{
    printf("deleting real_vector_double_t at %p\n", p );
    free( p->data );
    free( p );
}

complex_vector_single_t* new_complex_vector_single_t(unsigned length)
{
    
    CONSTRUCTOR_TEMPLATE(complex_vector_single_t, float)
    c->data = calloc( 2 * c->meta_data.vector_length ,
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

complex_vector_double_t* new_complex_vector_double_t(unsigned length)
{
    
    CONSTRUCTOR_TEMPLATE(complex_vector_double_t, double)    
    c->data = calloc( 2 * c->meta_data.vector_length,
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

