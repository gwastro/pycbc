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

#include "pycbccpu.h"
#include "pycbccpu_private.h"

unsigned pycbccpu_err_stash = 0;
char pycbccpu_err_map[][ERR_STRING_LEN] = {
    "No Error", 
    "Memory Error", 
    "Undefined Error"};
// ... extend error map on demand

int pycbc_cpu_check_err_occurred()
{
    //printf("debug: pycbc_cpu_check_err_occurred\n");
    if (!pycbccpu_err_stash) {
        return 0;
    }
    else {
        return 1;
    }
}

char* pycbc_cpu_get_err_message() 
{
    //printf("debug: pycbc_cpu_get_err_message\n");    
    return pycbccpu_err_map[pycbccpu_err_stash];
}

void pycbc_cpu_set_error(unsigned err)
{
    //printf("debug: pycbc_cpu_set_error\n"); 
    pycbccpu_err_stash = err;
}

void pycbc_cpu_clear_error()
{
    //printf("debug: pycbc_cpu_clear_error\n"); 
    pycbccpu_err_stash = 0;
}




cpu_context_t* new_cpu_context_t(unsigned device_id)
{
    cpu_context_t* c;
    
    c = (cpu_context_t*) malloc(sizeof(cpu_context_t));

    c->device_id = device_id;
    c->set_error = pycbc_cpu_set_error;

    // testing error handler
    // c->set_error(1);
    
    return c;
}

void delete_cpu_context_t( cpu_context_t* p )
{
    free( p );
}

