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
// pycbc type definitions for pycbc

#ifndef PYCBCOPENCL_TYPES_H
#define PYCBCOPENCL_TYPES_H

#include <stdlib.h>
#include <CL/opencl.h>

typedef struct
{
    unsigned       device_id;
    //cl_platform_id platform;
    //cl_device_id   device;
    //cl_context     context;
 
    // ... opencl context elements 

    int (*err_occurred)(void);
    char* (*err_message)(void);
    void (*set_error)(unsigned);

}
cl_context_t;

#endif /* PYCBCOPENCL_TYPES_H */
