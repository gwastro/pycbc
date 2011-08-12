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

//#include "../except.h"

#include "pycbcopencl_types.h"
#include "pycbcopencl_prototypes.h"

cl_context_t* new_cl_context_t(unsigned device_id)
{
    cl_context_t* c;
    
    c       = (cl_context_t*) malloc(sizeof(cl_context_t));
    
    printf("<#message#>");
    
    //if(c == NULL)
    //    throw(MemoryError);
    
    return c;
}

void delete_cl_context_t( cl_context_t* p )
{
    printf("destructed");
    free( p );
}

