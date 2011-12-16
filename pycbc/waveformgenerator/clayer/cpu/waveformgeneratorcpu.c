// Copyright (C) 2011 Karsten Wiesner, Drew Keppel
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
// Interface to algorithms to compute waveform related things.
// FIXME: Need to document the functions here!!!


#include <stdio.h>
#include "datavectorcpu.h"
#include "pycbccpu.h"
#include "waveformgeneratorcpu_private.h"
#include "waveformgeneratorcpu.h"
//#include <lal/LALSimInspiral.h>

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

#include "taylorf2.c"

int gen_precon_vector_TaylorF2(
    real_vector_single_cpu_t* precon_vec
    )
{
    taylorf2_precondition_factor(precon_vec);
    
    return 1;  // Yes TaylorF2 _has_ precondition capability and data
               // otherwise we would return 0 here to inform the framework
               // to not do preconditioning the data
}

void gen_waveform_filter_TaylorF2(
    complex_vector_single_cpu_t* waveform_filter,
    double M,
    double eta,
    int order,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three
    )
{
    taylorf2_phasing(waveform_filter, M, eta, order, f_min, f_max, minus_one_by_three);
}

void gen_waveform_strain_TaylorF2(
    complex_vector_single_cpu_t* waveform_strain,
    double M,
    double eta,
    int order,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three,
    real_vector_single_cpu_t* precon_vec
    )
{

    int i;
    taylorf2_phasing(waveform_strain, M, eta, order, f_min, f_max, minus_one_by_three);

    for(i=0; i < waveform_strain->meta_data.vector_length; i++)
    {
        waveform_strain->data[i] *= precon_vec->data[i];
    }
}

// TaylorT1 section ------------------------------------------------------------

/**
 * Function to generate the precondition vector for TaylorT1
 *
 * This is actually a pass as there is no precondition vector for this approximant
 */
int gen_precon_vector_TaylorT1(
    real_vector_single_cpu_t* precon_vec /**< pointer to data for the precondition vector */
    )
{
    return 0; //TaylorT1 does not have precondition capability and data.
}

/**
 */
void gen_waveform_filter_TaylorT1(
    complex_vector_single_cpu_t* waveform_filter,
    double m1,
    double m2,
    int order,
    double f_min,
    double f_max
    )
{
/*
    int i;
    REAL8TimeSeries **hplus;
    REAL8TimeSeries **hcross;
    int ret = XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, 0., 0., waveform_filter->meta_data->delta_x, m1, m2, f_min, 1e6*LAL_PC_SI, 0., amplitudeOrder, phaseOrder);
    if ret == XLAL_FAIULRE
    {
        //FIXME: We need to handle this failure somehow
    }
    memset(waveform_filter->data, 0, waveform_filter->meta_data.vector_length * waveform_filter->meta_data.element_size_bytes));
    for(i = 0; i < hplus->data->length; i++)
    {
        if(i >= waveform_filter->meta_data.vector_length)
        {
            //FIXME: This should be an error. The generated waveform overfills the space we have allocated.
        }
        __real__ waveform_filter->data[i] = hplus->data->data[i];
        __imag__ waveform_filter->data[i] = hcross->data->data[i];
    }
    LALFree(hplus);
    LALFree(hcross);
 */
}

/** 
 */
void gen_waveform_strain_TaylorT1(
    complex_vector_single_cpu_t* waveform_strain,
    double M,
    double eta,
    int order,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three,
    real_vector_single_cpu_t* precon_vec
    )
{
/*
    int i;
    REAL8TimeSeries **hplus;
    REAL8TimeSeries **hcross;
    int ret = XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, 0., 0., waveform_filter->meta_data->delta_x, m1, m2, f_min, 1e6*LAL_PC_SI, 0., amplitudeOrder, phaseOrder);
    if ret == XLAL_FAIULRE
    {
        //FIXME: We need to handle this failure somehow
    }
    memset(waveform_filter->data, 0, waveform_filter->meta_data.vector_length * waveform_filter->meta_data.element_size_bytes));
    for(i = 0; i < hplus->data->length; i++)
    {
        if(i >= waveform_filter->meta_data.vector_length)
        {
            //FIXME: This should be an error. The generated waveform overfills the space we have allocated.
        }
        __real__ waveform_filter->data[i] = hplus->data->data[i];
        __imag__ waveform_filter->data[i] = hcross->data->data[i];
    }
    LALFree(hplus);
    LALFree(hcross);
 */
}
