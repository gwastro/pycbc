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
// datavectorcpu declarations that are not going to be swig wrapped 
// thus they are private property of the clayer


#ifndef DATAVECTORCPU_PRIVATE_H
#define DATAVECTORCPU_PRIVATE_H

#include <stdlib.h>


// prototypes of constructors/destructors.  

real_vector_single_cpu_t* new_real_vector_single_cpu_t(unsigned long, double);
void delete_real_vector_single_cpu_t( real_vector_single_cpu_t* );

real_vector_double_cpu_t* new_real_vector_double_cpu_t(unsigned long, double);
void delete_real_vector_double_cpu_t( real_vector_double_cpu_t* );

complex_vector_single_cpu_t* new_complex_vector_single_cpu_t(unsigned long, double);
void delete_complex_vector_single_cpu_t( complex_vector_single_cpu_t* );

complex_vector_double_cpu_t* new_complex_vector_double_cpu_t(unsigned long, double);
void delete_complex_vector_double_cpu_t( complex_vector_double_cpu_t* );

#endif /* DATAVECTORCPU_PRIVATE_H */
