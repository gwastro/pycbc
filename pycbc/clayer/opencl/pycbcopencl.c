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
// pycbc's constructors and destructors implementation for pycbc

#include <stdio.h>
#include <string.h>
#include "pycbcopencl_types.h"
#include "pycbcopencl_prototypes.h"
#include "gpu_inspiral_gpuutils.h"

//#include "gclUtils.h"

unsigned pycbcopencl_error_stash = 0;
char pycbcopencl_error_message[ERROR_STRING_LEN];

int pycbc_opencl_check_error()
{
    //printf("debug: pycbc_opencl_check_error\n");
    return pycbcopencl_error_stash;
}

char* pycbc_opencl_get_error_message()
{
    //printf("debug: pycbc_opencl_get_error_message\n");

    return pycbcopencl_error_message;
    //return pycbcopencl_err_map[pycbcopencl_err_stash];
}

void pycbc_opencl_set_error(int error_id, char* error_message)
{
    //printf("debug: pycbc_opencl_set_error\n"); 
    pycbcopencl_error_stash = error_id;
    strncpy(pycbcopencl_error_message, error_message, ERROR_STRING_LEN);
}

void pycbc_opencl_clear_error()
{
    //printf("debug: pycbc_opencl_clear_error\n"); 
    pycbcopencl_error_stash = 0;
}


cl_context_t* new_cl_context_t(unsigned device_id)
{
    cl_context_t* c;

    c = (cl_context_t*) malloc(sizeof(cl_context_t));

    c->device_id = device_id;
    c->set_error = pycbc_opencl_set_error;
    
    // this will update c with all OpenCl context elements
    int err = gpuinsp_InitGPU(c, device_id);
    
    printf("Error%d",err);
    
    if(err)
      c->set_error( PYCBC_RUNTIME_ERROR, "gpuinsp_InitGPU failed");

    // just testing the error handling
    // c->set_error( PYCBC_RUNTIME_ERROR, "gpuinsp_InitGPU failed");
    
    return c;
}

void delete_cl_context_t( cl_context_t* p )
{
    free( p );
}

