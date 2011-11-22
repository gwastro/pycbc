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
// datavector type definitions for pycbc

#ifndef DATAVECTOROPENCL_H
#define DATAVECTOROPENCL_H

#include <stdlib.h>
#include "../datavector.h"
#include <datavectorcpu.h>
#include "datavectoropencl.h"
#include "pycbcopencl.h"
#include <CL/opencl.h>

typedef struct
{
    meta_data_t meta_data;
    cl_mem data;
}
real_vector_single_opencl_t;

typedef struct
{
    meta_data_t meta_data;
    cl_mem data;
}
real_vector_double_opencl_t;

typedef struct
{
    meta_data_t meta_data;
    cl_mem real_data;
    cl_mem imag_data;
}
complex_vector_single_opencl_t;

typedef struct
{
    meta_data_t meta_data;
    cl_mem real_data;
    cl_mem imag_data;
}
complex_vector_double_opencl_t;


// transfer functions
void transfer_real_vector_single_from_cpu(cl_context_t *,
                                          real_vector_single_opencl_t,
                                          real_vector_single_cpu_t);

void transfer_real_vector_single_to_cpu(cl_context_t *,
                                         real_vector_single_cpu_t,
                                         real_vector_single_opencl_t);

void transfer_real_vector_double_from_cpu(cl_context_t *,
                                          real_vector_double_opencl_t,
                                          real_vector_double_cpu_t);

void transfer_real_vector_double_to_cpu( cl_context_t *,
                                         real_vector_double_cpu_t,
                                         real_vector_double_opencl_t);



void transfer_complex_vector_single_from_cpu(cl_context_t *,
                                             complex_vector_single_opencl_t,
                                             complex_vector_single_cpu_t);
void transfer_complex_vector_single_to_cpu( cl_context_t*,
                                            complex_vector_single_cpu_t,
                                            complex_vector_single_opencl_t);

void transfer_complex_vector_double_from_cpu(cl_context_t*,
                                             complex_vector_double_opencl_t,
                                             complex_vector_double_cpu_t);
void transfer_complex_vector_double_to_cpu( cl_context_t*,
                                           complex_vector_double_cpu_t,
                                           complex_vector_double_opencl_t);


#endif /* DATAVECTOROPENCL_H */
