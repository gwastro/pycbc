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
// waveformgeneratorcpu declarations that are not going to be swig wrapped 
// thus they are private property of the clayer

#ifndef WAVEFORMGENERATORCPU_PRIVATE_H
#define WAVEFORMGENERATORCPU_PRIVATE_H

#include <stdlib.h>
#include "pycbccpu.h"
#include "datavectorcpu.h"
#include "waveformgeneratorcpu.h"

// prototypes of all methodes that will extend pure c typedefs
waveform_generator_cpu_t* new_waveform_generator_cpu_t(cpu_context_t*);
void delete_waveform_generator_cpu_t( waveform_generator_cpu_t* );

#endif /* WAVEFORMGENERATORCPU_PRIVATE_H */
