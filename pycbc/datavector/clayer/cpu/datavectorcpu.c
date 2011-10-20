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

#include "pycbccpu.h"
#include "datavectorcpu.h"
#include "datavectorcpu_private.h"

real_vector_single_cpu_t* new_real_vector_single_cpu_t(cpu_context_t* context,
                                                       unsigned long length, 
                                                       double delta_x)
{
    CONSTRUCTOR_TEMPLATE(real_vector_single_cpu_t, float)

    c->data = (float*)calloc( c->meta_data.vector_length , 
                      c->meta_data.element_size_bytes );
    return c;
}

void delete_real_vector_single_cpu_t( real_vector_single_cpu_t* p )
{
    free( p->data );
    free( p );
}

real_vector_double_cpu_t* new_real_vector_double_cpu_t(cpu_context_t* context,
                                                       unsigned long length, 
                                                       double delta_x)
{
    CONSTRUCTOR_TEMPLATE(real_vector_double_cpu_t, double)
    c->data = (double*)calloc( c->meta_data.vector_length , 
                      c->meta_data.element_size_bytes );
    return c;
}

void delete_real_vector_double_cpu_t( real_vector_double_cpu_t* p )
{
    free( p->data );
    free( p );
}

complex_vector_single_cpu_t* new_complex_vector_single_cpu_t(cpu_context_t* context,
                                                             unsigned long length, 
                                                             double delta_x)
{
    CONSTRUCTOR_TEMPLATE(complex_vector_single_cpu_t, complex float)
    c->data = (complex float *)calloc(c->meta_data.vector_length ,
                          c->meta_data.element_size_bytes );
    return c;
}

void delete_complex_vector_single_cpu_t( complex_vector_single_cpu_t* p )
{

    free( p->data );
    free( p );
}

complex_vector_double_cpu_t* new_complex_vector_double_cpu_t(cpu_context_t* context,
                                                             unsigned long length, 
                                                             double delta_x)
{
    CONSTRUCTOR_TEMPLATE(complex_vector_double_cpu_t, complex double)    
    c->data = (complex double *)calloc(c->meta_data.vector_length,
                          c->meta_data.element_size_bytes );
    return c;
}

void delete_complex_vector_double_cpu_t( complex_vector_double_cpu_t* p )
{
    free( p->data );
    free( p );
}

