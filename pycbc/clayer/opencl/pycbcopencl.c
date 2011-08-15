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

#include "../except.h"

#include "pycbcopencl_types.h"
#include "pycbcopencl_prototypes.h"

//#include "gclUtils.h"


cl_context_t* new_cl_context_t(unsigned device_id)
{
    int err;
    cl_context_t* c;
    
    c = (cl_context_t*) malloc(sizeof(cl_context_t));
    
  //  gclInitErrorMessages();
  //  gclInitLogMessages();
  //  gclLogLevel = gclINFO;
  //  gclLog(gclINFO, "Initializing GPU environment...");
  //  err = gclInitGPU(gclFirstAvailable);
    
    
    // all these globals should become an element of cl_context_t
    /*
    cl_int                 gcl_err             = CL_SUCCESS;
    cl_uint                gcl_numPlatforms;
    cl_platform_id         gcl_platform        = NULL;
    cl_platform_id*        gcl_platforms       = new cl_platform_id[gclMaxPlatforms];
    size_t                 gcl_devicenumber;
    cl_device_id*          gcl_devices;
    cl_device_id           gcl_available_device= NULL;
    cl_context_properties  gcl_cps[3]          = {CL_CONTEXT_PLATFORM, (cl_context_properties) gcl_platform, 0};
    cl_context_properties* gcl_cprops;
    cl_context             gcl_context;
    cl_program             gcl_programs[gclMaxPrograms]; 
    */
    
   // if(err != CL_SUCCESS)
   //     printf("OpenCl init error: %d", err);
    
    
    // *** currently buggy - I'm working on it. K. Wiesner *** 
    //if(c == NULL)
    //    throw(MemoryError, "Could not allocate cl_context_t")
    
    return c;
}

void delete_cl_context_t( cl_context_t* p )
{
    free( p );
}

