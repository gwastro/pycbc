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
// matchedfilteropencl declarations 

#ifndef MATCHEDFILTEROPENCL_H
#define MATCHEDFILTEROPENCL_H

#include <stdlib.h>
#include <CL/opencl.h>

typedef struct
{
    cl_program program;

    cl_mem cl_stilde;
    cl_mem cl_htilde;

    cl_kernel gpu_snr_product;
}
matched_filter_opencl_t;


extern "C" void gen_snr_opencl(cl_context_t* context,
                    matched_filter_opencl_t * matchedfilter,
                    real_vector_single_opencl_t* snr,
                    complex_vector_single_opencl_t* stilde, 
                    complex_vector_single_opencl_t* htilde);


#endif /* MATCHEDFILTEROPENCL_H */
