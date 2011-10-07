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
#include <complex.h>
#include "pycbcopencl_types.h"
#include "gpu_inspiral_gpuutils.h"
#include "datavectoropencl.h"
#include "datavectoropencl_private.h"

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

void delete_real_vector_double_opencl_t( real_vector_double_opencl_t* p )
{

    if (p->data) clReleaseMemObject(p->data);
    free( p );
}

void transfer_real_vector_double_from_cpu( cl_context_t * context,
                                           real_vector_double_opencl_t dest,
                                           real_vector_double_cpu_t src )
{
    int err;
    err = clEnqueueWriteBuffer(context->io_queue, dest.data, CL_TRUE, 0, src.meta_data.element_size_bytes*src.meta_data.vector_length, src.data, 0, NULL, NULL);
    clFinish(context->io_queue);
}

void transfer_real_vector_double_to_cpu( cl_context_t * context,
                                         real_vector_double_cpu_t dest,
                                         real_vector_double_opencl_t src )
{
    int err;
    err = clEnqueueReadBuffer(context->io_queue, src.data, CL_TRUE, 0, src.meta_data.element_size_bytes*src.meta_data.vector_length, dest.data, 0, NULL, NULL);
    clFinish(context->io_queue);
}


// complex_vector_single methods -----------------------------------------------

complex_vector_single_opencl_t* new_complex_vector_single_opencl_t(cl_context_t * context,
                                                        unsigned long length,
                                                        double delta_x)
{
    int err;
    CONSTRUCTOR_TEMPLATE(complex_vector_single_opencl_t, float)

    c->real_data = clCreateBuffer(context->context, CL_MEM_READ_WRITE, sizeof(float)*length, NULL, &err);
    if (gpuinsp_checkError(err,"Allocating real part of complex float opencl vector on GPU") !=0 )
           return NULL;

    c->imag_data = clCreateBuffer(context->context, CL_MEM_READ_WRITE, sizeof(float)*length, NULL, &err);
    if (gpuinsp_checkError(err,"Allocating imag part of complex float opencl vector on GPU") !=0 )
           return NULL;
    return c;
}

void delete_complex_vector_single_opencl_t( complex_vector_single_opencl_t* p )
{
    if (p->real_data) clReleaseMemObject(p->real_data);
    if (p->imag_data) clReleaseMemObject(p->imag_data);
    free( p );
}

void transfer_complex_vector_single_from_cpu( cl_context_t* context,
                                            complex_vector_single_opencl_t dest,
                                            complex_vector_single_cpu_t src )
{

   int err;
   if (dest.meta_data.vector_length != src.meta_data.vector_length)
        return; // ERROR!

   size_t buffersize = dest.meta_data.element_size_bytes*dest.meta_data.vector_length;
   float*   cpureal = (float*) malloc(buffersize);
   float*   cpuimag = (float*) malloc(buffersize);

   if (cpureal == NULL || cpuimag == NULL ) {
       printf("Cannot allocate temporary buffer on CPU.\n");
       goto cleanup;
   }
   for (unsigned long i = 0; i < src.meta_data.vector_length; i++) {
     cpureal[i] = __real__ src.data[i];
     cpuimag[i] = __imag__ src.data[i];
   }

   err = clEnqueueWriteBuffer(context->io_queue, dest.real_data, CL_TRUE, 0, buffersize, cpureal, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex single vector real part from CPU") !=0 ) goto cleanup;

   err = clEnqueueWriteBuffer(context->io_queue, dest.imag_data, CL_TRUE, 0, buffersize, cpuimag, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex single vector imag part from CPU") !=0 ) goto cleanup;

   clFinish(context->io_queue);

cleanup:
   free(cpureal);
   free(cpuimag);
}


void transfer_complex_vector_single_to_cpu( cl_context_t* context,
                                            complex_vector_single_cpu_t dest,
                                            complex_vector_single_opencl_t src )
{

   int err;
   if (dest.meta_data.vector_length != src.meta_data.vector_length)
       return; // ERROR!

   long int buffersize = src.meta_data.element_size_bytes*src.meta_data.vector_length;
   float*   cpureal = (float*) malloc(buffersize);
   float*   cpuimag = (float*) malloc(buffersize);

   if (cpureal == NULL || cpuimag == NULL ) {
       printf("Cannot allocate temporary buffer on CPU.\n");
       goto cleanup;
   }

   err = clEnqueueReadBuffer(context->io_queue, src.real_data, CL_TRUE, 0, buffersize, cpureal, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex single vector real part from GPU") !=0 ) goto cleanup;

   err |= clEnqueueReadBuffer(context->io_queue, src.imag_data, CL_TRUE, 0, buffersize, cpuimag, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex single vector imag part from GPU") !=0 ) goto cleanup;
   clFinish(context->io_queue);

   for (unsigned long i = 0; i < src.meta_data.vector_length; i++) {
     __real__ dest.data[i] = cpureal[i];
     __imag__ dest.data[i] = cpuimag[i];
   }

cleanup:
   free(cpureal);
   free(cpuimag);
}


// complex_vector_double methods -----------------------------------------------

complex_vector_double_opencl_t* new_complex_vector_double_opencl_t(cl_context_t * context,
                                                        unsigned long length,
                                                        double delta_x)
{

    int err;
    CONSTRUCTOR_TEMPLATE(complex_vector_double_opencl_t, double)

    c->real_data = clCreateBuffer(context->context, CL_MEM_READ_WRITE, sizeof(double)*length, NULL, &err);
    if (gpuinsp_checkError(err,"Allocating real part of complex double opencl vector on GPU") !=0 )
           return NULL;

    c->imag_data = clCreateBuffer(context->context, CL_MEM_READ_WRITE, sizeof(double)*length, NULL, &err);
    if (gpuinsp_checkError(err,"Allocating imag part of complex double opencl vector on GPU") !=0 )
           return NULL;
    return c;
}

void delete_complex_vector_double_opencl_t( complex_vector_double_opencl_t* p )
{
    if (p->real_data) clReleaseMemObject(p->real_data);
    if (p->imag_data) clReleaseMemObject(p->imag_data);
    free( p );
}

void transfer_complex_vector_double_from_cpu( cl_context_t* context,
                                            complex_vector_double_opencl_t dest,
                                            complex_vector_double_cpu_t src )
{
   int err;
   if (dest.meta_data.vector_length != src.meta_data.vector_length)
        return; // ERROR!

   size_t buffersize = dest.meta_data.element_size_bytes*dest.meta_data.vector_length;
   double*   cpureal = (double*) malloc(buffersize);
   double*   cpuimag = (double*) malloc(buffersize);

   if (cpureal == NULL || cpuimag == NULL ) {
       printf("Cannot allocate temporary buffer on CPU.\n");
       goto cleanup;
   }
   for (unsigned long i = 0; i < src.meta_data.vector_length; i++) {
     cpureal[i] = __real__ src.data[i];
     cpuimag[i] = __imag__ src.data[i];
   }

   err = clEnqueueWriteBuffer(context->io_queue, dest.real_data, CL_TRUE, 0, buffersize, cpureal, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex double vector real part from CPU") !=0 ) goto cleanup;

   err = clEnqueueWriteBuffer(context->io_queue, dest.imag_data, CL_TRUE, 0, buffersize, cpuimag, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex double vector imag part from CPU") !=0 ) goto cleanup;

   clFinish(context->io_queue);

cleanup:
   free(cpureal);
   free(cpuimag);
}

void transfer_complex_vector_double_to_cpu(cl_context_t* context,
                                           complex_vector_double_cpu_t dest,
                                           complex_vector_double_opencl_t src )
{

   int err;
   if (dest.meta_data.vector_length != src.meta_data.vector_length)
       return; // ERROR!

   long int buffersize = src.meta_data.element_size_bytes*src.meta_data.vector_length;
   double*   cpureal = (double*) malloc(buffersize);
   double*   cpuimag = (double*) malloc(buffersize);

   if (cpureal == NULL || cpuimag == NULL ) {
       printf("Cannot allocate temporary buffer on CPU.\n");
       goto cleanup;
   }

   err = clEnqueueReadBuffer(context->io_queue, src.real_data, CL_TRUE, 0, buffersize, cpureal, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex double vector real part from GPU") !=0 ) goto cleanup;

   err |= clEnqueueReadBuffer(context->io_queue, src.imag_data, CL_TRUE, 0, buffersize, cpuimag, 0, NULL, NULL);
   if (gpuinsp_checkError(err,"Transfering complex double vector imag part from GPU") !=0 ) goto cleanup;
   clFinish(context->io_queue);

   for (unsigned long i = 0; i < src.meta_data.vector_length; i++) {
     __real__ dest.data[i] = cpureal[i];
     __imag__ dest.data[i] = cpuimag[i];
   }

cleanup:
   free(cpureal);
   free(cpuimag);
}

