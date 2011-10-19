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
// matchedfiltercpu declarations that are not going to be swig wrapped 
// thus they are private property of the clayer

#ifndef MATCHEDFILTERCPU_PRIVATE_H
#define MATCHEDFILTERCPU_PRIVATE_H

#include <stdlib.h>
#include "pycbccpu.h"
#include "datavectorcpu.h"
#include "matchedfiltercpu_private.h"
#include "matchedfiltercpu.h"


// prototypes of all methodes that will extend pure c typedefs
matched_filter_cpu_t* new_matched_filter_cpu_t(cpu_context_t*);
void delete_matched_filter_cpu_t( matched_filter_cpu_t* );


void correlate_complex_freq_vectors( complex_vector_single_cpu_t* out,
				     complex_vector_single_cpu_t* x, 
				     complex_vector_single_cpu_t* y, 
				     double f_min);

#endif /* MATCHEDFILTERCPU_PROTOTYPES_H */
