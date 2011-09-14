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

/* pycbc error ids - python exceptions will be raised accordingly */
#define PYCBC_NO_ERROR            0        /* no error */
#define PYCBC_ATTRIBUTE_ERROR     1        /* invalid attribute accessed */
#define PYCBC_EOF_ERROR           2        /* end of file on io */
#define PYCBC_IO_ERROR            3        /* io error */
#define PYCBC_INDEX_ERROR         4        /* subscript out of range */
#define PYCBC_TYPE_ERROR          5        /* invalid type */
#define PYCBC_VALUE_ERROR         6        /* right type, wrong value */
#define PYCBC_MEMORY_ERROR        7        /* out of memory */
#define PYCBC_NAME_ERROR          8        /* name not found */
#define PYCBC_OVERFLOW_ERROR      9        /* artithmetic operation too large */
#define PYCBC_ZERO_DIVISION_ERROR 10       /* dive by zero */
#define PYCBC_RUNTIME_ERROR       11       /* everything else */

typedef struct
{
    unsigned         device_id;
    cl_platform_id   platform;
    cl_device_id     device;
    cl_context       context;
    cl_command_queue kernel_queue;
    cl_command_queue io_queue;
    cl_program       program;

  void (*set_error)(int, char*);
}
cl_context_t;

#endif /* PYCBCOPENCL_TYPES_H */
