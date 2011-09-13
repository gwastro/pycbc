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
#include "pycbcopencl_types.h"
#include "pycbcopencl_prototypes.h"
#include "gpu_inspiral_gpuutils.h"

//#include "gclUtils.h"

unsigned pycbcopencl_err_stash = 0;
char pycbcopencl_err_map[][ERR_STRING_LEN] = {
    "No Error", 
    "Memory Error", 
    "Undefined Error"};
// ... extend error map on demand

int pycbc_opencl_check_err_occurred()
{
    //printf("debug: pycbc_opencl_check_err_occurred\n");
    if (!pycbcopencl_err_stash) {
        return 0;
    }
    else {
        return 1;
    }
}

char* pycbc_opencl_get_err_message()
{
    //printf("debug: pycbc_opencl_get_err_message\n");  
    return pycbcopencl_err_map[pycbcopencl_err_stash];
}

void pycbc_opencl_set_error(unsigned err)
{
    //printf("debug: pycbc_opencl_set_error\n"); 
    pycbcopencl_err_stash = err;
}

void pycbc_opencl_clear_error()
{
    //printf("debug: pycbc_opencl_clear_error\n"); 
    pycbcopencl_err_stash = 0;
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
      c->set_error(err);
    
    return c;
}

void delete_cl_context_t( cl_context_t* p )
{
    free( p );
}

