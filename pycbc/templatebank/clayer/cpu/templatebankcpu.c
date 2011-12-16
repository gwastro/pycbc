// Copyright (C) 2011 Duncan Brown, Karsten Wiesner
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

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <complex.h>
#include "datavectorcpu.h"
#include "pycbccpu.h"
#include "templatebankcpu_private.h"
#include "templatebankcpu.h"


/* 
 * functions
 */


template_bank_cpu_t* new_template_bank_cpu_t(cpu_context_t* context)
{
    
    template_bank_cpu_t* c;
    c = (template_bank_cpu_t*) malloc( sizeof(template_bank_cpu_t) );
    
    return c;
}


void delete_template_bank_cpu_t( template_bank_cpu_t* p )
{
    free( p );
}


