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
// datavector constructor destructor prototypes for pycbc

#ifndef DATAVECTOROPENCL_PROTOTYPES_H
#define DATAVECTOROPENCL_PROTOTYPES_H

#include <stdlib.h>

// prototypes of all methodes that will extend pure c typedefs
real_vector_single_opencl_t* new_real_vector_single_opencl_t(unsigned long, double);
void delete_real_vector_single_opencl_t( real_vector_single_opencl_t* );

real_vector_double_opencl_t* new_real_vector_double_opencl_t(unsigned long, double);
void delete_real_vector_double_opencl_t( real_vector_double_opencl_t* );

complex_vector_single_opencl_t* new_complex_vector_single_opencl_t(unsigned long, double);
void delete_complex_vector_single_opencl_t( complex_vector_single_opencl_t* );

complex_vector_double_opencl_t* new_complex_vector_double_opencl_t(unsigned long, double);
void delete_complex_vector_double_opencl_t( complex_vector_double_opencl_t* );

#endif /* DATAVECTOROPENCL_PROTOTYPES_H */
