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
// waveformgeneratorcpu extension functions


#include <stdio.h>
#include "datavectorcpu.h"
#include "pycbccpu.h"
#include "waveformgeneratorcpu_private.h"
#include "waveformgeneratorcpu.h"

waveform_generator_cpu_t* new_waveform_generator_cpu_t(cpu_context_t* context)
{
    
    waveform_generator_cpu_t* c;
    c = (waveform_generator_cpu_t*) malloc( sizeof(waveform_generator_cpu_t) );
    
    return c;
}


void delete_waveform_generator_cpu_t( waveform_generator_cpu_t* p )
{
    free( p );
}


// TaylorF2 section ------------------------------------------------------------

int gen_precon_vector_TaylorF2(real_vector_single_cpu_t* precon_vec)
{    

    int i;
    
    //printf("called clayer: gen_precon_vector_TaylorF2");
    
    for (i=0; i < precon_vec->meta_data.vector_length; i++) {
        precon_vec->data[i]=i;
    }

    
    return 1;  // Yes TaylorF2 _has_ precondition capability and data
               // otherwise we would return 0 here to inform the framework
               // to not do preconditioning the data
}

void gen_waveform_filter_TaylorF2( complex_vector_single_cpu_t* waveform_filter,
                                   float mass1, float mass2 )
{    
    
    int i;
    
    //printf("called clayer: gen_waveform_filter_TaylorF2");
    
    // prototyping
    
    for (i=0; i < waveform_filter->meta_data.vector_length; i++) {
        __real__ waveform_filter->data[i]=i*mass1;
        __imag__ waveform_filter->data[i]=i*mass2;
    }
    
    
}


