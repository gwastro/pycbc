// Copyright (C) 2011 Karsten Wiesner, Josh Willis
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
// straindatacpu declarations 

#ifndef STRAINDATACPU_H
#define STRAINDATACPU_H

#include <stdlib.h>

typedef struct
{
    int nothing_to_define_yet;
}
strain_data_cpu_t;

//prototypes of all functions that swig wraps to methods

int copy_subvector_complex_double_cpu(complex_vector_double_cpu_t *src, complex_vector_double_cpu_t *dst,
				      unsigned long offset, unsigned long length);
int copy_subvector_complex_single_cpu(complex_vector_double_cpu_t *src, complex_vector_double_cpu_t *dst,
				      unsigned long offset, unsigned long length);
int copy_subvector_real_double_cpu(complex_vector_double_cpu_t *src, complex_vector_double_cpu_t *dst,
				   unsigned long offset, unsigned long length);
int copy_subvector_real_single_cpu(complex_vector_double_cpu_t *src, complex_vector_double_cpu_t *dst,
				   unsigned long offset, unsigned long length);

void* fftw_generate_plan(unsigned long length, real_vector_single_cpu_t* in_tmp,
                         complex_vector_single_cpu_t* out_tmp, char* sign, char* style);

int fftw_transform_segments(void* plan, real_vector_single_cpu_t* in_buf, 
                            unsigned long input_buf_offset, 
                            complex_vector_single_cpu_t* out_buf);

int frame_cpp_read_frames(real_vector_double_cpu_t* out_buf, char* channel_name, 
                          unsigned long gps_start_time, unsigned long gps_end_time, 
                          char* cache_url);


#endif /* STRAINDATACPU_H */
