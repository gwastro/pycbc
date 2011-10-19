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
// pycbccpu declarations that are not going to be swig wrapped 
// thus they are private property of the clayer

#ifndef PYCBCCPU_PRIVATE_H
#define PYCBCCPU_PRIVATE_H

#include <stdlib.h>

#define ERR_STRING_LEN 512 

extern unsigned pycbccpu_err_stash;
extern char pycbccpu_err_map[][ERR_STRING_LEN];

int pycbc_cpu_check_err_occurred(void);
char* pycbc_cpu_get_err_message(void);
void pycbc_cpu_set_error(unsigned);
void pycbc_cpu_clear_error(void);

// prototypes of all methodes that will extend pure c typedefs
cpu_context_t* new_cpu_context_t(unsigned);
void delete_cpu_context_t( cpu_context_t* );


#endif /* PYCBCCPU_PRIVATE_H */
