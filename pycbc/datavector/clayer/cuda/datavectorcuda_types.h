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

#ifndef DATAVECTORCUDA_TYPES_H
#define DATAVECTORCUDA_TYPES_H

#include <stdlib.h>
#include "../datavector_types.h"

typedef struct
{
    meta_data_t meta_data;
    float *data;
}
real_vector_single_cuda_t;

typedef struct
{
    meta_data_t meta_data;
    double *data;
}
real_vector_double_cuda_t;

typedef struct
{
    meta_data_t meta_data;
    float *data;
}
complex_vector_single_cuda_t;

typedef struct
{
    meta_data_t meta_data;
    double *data;
}
complex_vector_double_cuda_t;

#endif /* DATAVECTORCUDA_TYPES_H */
