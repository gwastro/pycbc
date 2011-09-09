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

#include "pycbccpu_types.h"
#include "pycbccpu_prototypes.h"

//#include "gclUtils.h"

unsigned pycbccpu_err_stash = 0;
char pycbccpu_err_map[][ERR_STRING_LEN] = {
    "No Error", 
    "Memory Error", 
    "berr", 
    "cerr", 
    "Undefined Error"};


int pycbc_err_occurred()
{
    printf("\ndebug: pycbc_err_occurred\n");
    
    if (!pycbccpu_err_stash) {
        return 0;
    }
    else {
        return 1;
    }
}

char* pycbc_err_message() 
{
    unsigned tmp;
    printf("debug: pycbc_err_message");    

    tmp = pycbccpu_err_stash;
    pycbccpu_err_stash = 0;
    
    return pycbccpu_err_map[tmp];
}

void pycbc_set_error(unsigned err)
{
    printf("debug: pycbc_set_error"); 
    pycbccpu_err_stash = err;
}

cpu_context_t* new_cpu_context_t(unsigned device_id)
{
    cpu_context_t* c;
    
    c = (cpu_context_t*) malloc(sizeof(cpu_context_t));

    c->device_id = device_id;
      
    c->err_occurred = pycbc_err_occurred;
    c->err_message  = pycbc_err_message;
    c->set_error    = pycbc_set_error;
    
    // c->set_error(1);
    
    return c;
}

void delete_cpu_context_t( cpu_context_t* p )
{
    free( p );
}

