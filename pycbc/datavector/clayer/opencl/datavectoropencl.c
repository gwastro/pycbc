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
#include <pycbcopencl_types.h>
#include "gpu_inspiral_gpuutils.h"
#include "datavectoropencl_types.h"
#include "datavectoropencl_prototypes.h"

// real_vector_single methods --------------------------------------------------

real_vector_single_opencl_t* new_real_vector_single_opencl_t(cl_context_t * context,
                                                        unsigned long length,
                                                        double delta_x)
{
    int err;

    CONSTRUCTOR_TEMPLATE(real_vector_single_opencl_t, float)
    c->data = clCreateBuffer(context->context, CL_MEM_READ_WRITE, sizeof(float)*length, NULL, &err);
    if (gpuinsp_checkError(err,"Allocating real float opencl vector on GPU") !=0 )
           return NULL;
    return c;
}

void delete_real_vector_single_opencl_t( real_vector_single_opencl_t* p )
{

    if (p->data) clReleaseMemObject(p->data);
    free( p );
}

void transfer_real_vector_single_from_cpu( cl_context_t * context,
                                           real_vector_single_opencl_t dest,
                                           real_vector_single_cpu_t src )
{
    int err;
    err = clEnqueueWriteBuffer(context->io_queue, dest.data, CL_TRUE, 0, src.meta_data.element_size_bytes*src.meta_data.vector_length, src.data, 0, NULL, NULL);
    clFinish(context->io_queue);
}

void transfer_real_vector_single_to_cpu( cl_context_t * context,
                                         real_vector_single_cpu_t dest,
                                         real_vector_single_opencl_t src )
{
    int err;
    err = clEnqueueReadBuffer(context->io_queue, src.data, CL_TRUE, 0, src.meta_data.element_size_bytes*src.meta_data.vector_length, dest.data, 0, NULL, NULL);
    clFinish(context->io_queue);
}


// real_vector_double methods --------------------------------------------------

real_vector_double_opencl_t* new_real_vector_double_opencl_t(cl_context_t * context,
                                                        unsigned long length, 
                                                        double delta_x)
{
    int err;

    CONSTRUCTOR_TEMPLATE(real_vector_double_opencl_t, double)
    c->data = clCreateBuffer(context->context, CL_MEM_READ_WRITE, c->meta_data.element_size_bytes*c->meta_data.vector_length, NULL, &err);
    if (gpuinsp_checkError(err,"Allocating real float opencl vector on GPU") !=0 )
           return NULL;
    return c;
}

void delete_double_vector_double_opencl_t( real_vector_single_opencl_t* p )
{

    if (p->data) clReleaseMemObject(p->data);
    free( p );
}

void transfer_double_vector_double_from_cpu( cl_context_t * context,
                                           real_vector_double_opencl_t dest,
                                           real_vector_double_cpu_t src )
{
    int err;
    err = clEnqueueWriteBuffer(context->io_queue, dest.data, CL_TRUE, 0, src.meta_data.element_size_bytes*src.meta_data.vector_length, src.data, 0, NULL, NULL);
    clFinish(context->io_queue);
}

void transfer_double_vector_double_to_cpu( cl_context_t * context,
                                         real_vector_double_cpu_t dest,
                                         real_vector_double_opencl_t src )
{
    int err;
    err = clEnqueueReadBuffer(context->io_queue, src.data, CL_TRUE, 0, src.meta_data.element_size_bytes*src.meta_data.vector_length, dest.data, 0, NULL, NULL);
    clFinish(context->io_queue);
}


// complex_vector_single methods -----------------------------------------------

complex_vector_single_opencl_t* new_complex_vector_single_opencl_t(
                                                        unsigned long length, 
                                                        double delta_x)
{
    CONSTRUCTOR_TEMPLATE(complex_vector_single_opencl_t, float)

    c->real_data = (float*)calloc(c->meta_data.vector_length ,
                                  c->meta_data.element_size_bytes );
    c->imag_data = (float*)calloc(c->meta_data.vector_length ,
                                  c->meta_data.element_size_bytes );
    return c;
}

void delete_complex_vector_single_opencl_t( complex_vector_single_opencl_t* p )
{
    
    free( p->real_data );
    free( p->imag_data );
    free( p );
}

void transfer_complex_vector_single_from_cpu( 
                                            complex_vector_single_opencl_t dest, 
                                            complex_vector_single_cpu_t src )
{
    unsigned long i=0;
    
    if (dest.meta_data.vector_length != src.meta_data.vector_length)
        return; // ERROR!
    
    // prototyping the transfer (currently opencldata is cpudata)
    for (i=0; i < dest.meta_data.vector_length; i++) {
        *dest.real_data++ = __real__ *src.data;
        *dest.imag_data++ = __imag__ *src.data++;
    }
}

void transfer_complex_vector_single_to_cpu( complex_vector_single_cpu_t dest, 
                                            complex_vector_single_opencl_t src )
{
    unsigned long i=0;
    
    if (dest.meta_data.vector_length != src.meta_data.vector_length)
        return; // ERROR!
    
    // prototyping the transfer (currently opencldata is cpudata)
    for (i=0; i < dest.meta_data.vector_length; i++) {
        __real__ *dest.data   = *src.real_data++;
        __imag__ *dest.data++ = *src.imag_data++;
    }
}


// complex_vector_double methods -----------------------------------------------

complex_vector_double_opencl_t* new_complex_vector_double_opencl_t(
                                                        unsigned long length, 
                                                        double delta_x)
{
    CONSTRUCTOR_TEMPLATE(complex_vector_double_opencl_t, double)    

    c->real_data = (double*)calloc(c->meta_data.vector_length,
                                   c->meta_data.element_size_bytes );
    c->imag_data = (double*)calloc(c->meta_data.vector_length,
                                   c->meta_data.element_size_bytes );
    return c;
}

void delete_complex_vector_double_opencl_t( complex_vector_double_opencl_t* p )
{
    free( p->real_data );
    free( p->imag_data );
    free( p );
}

void transfer_complex_vector_double_from_cpu(
                                            complex_vector_double_opencl_t dest, 
                                            complex_vector_double_cpu_t src )
{
    unsigned long i=0;
    
    if (dest.meta_data.vector_length != src.meta_data.vector_length)
        return; // ERROR!
    
    // prototyping the transfer (currently opencldata is cpudata)
    for (i=0; i < dest.meta_data.vector_length; i++) {
        *dest.real_data++ = __real__ *src.data;
        *dest.imag_data++ = __imag__ *src.data++;
    }
}

void transfer_complex_vector_double_to_cpu( complex_vector_double_cpu_t dest, 
                                           complex_vector_double_opencl_t src )
{
    unsigned long i=0;
    
    if (dest.meta_data.vector_length != src.meta_data.vector_length)
        return; // ERROR!
    
    // prototyping the transfer (currently opencldata is cpudata)
    for (i=0; i < dest.meta_data.vector_length; i++) {
        __real__ *dest.data   = *src.real_data++;
        __imag__ *dest.data++ = *src.imag_data++;
    }
}

