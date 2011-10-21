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
// templatebankcpu declarations 


#ifndef TEMPLATEBANKCPU_H
#define TEMPLATEBANKCPU_H

#include <stdlib.h>

typedef struct
{
    int nothing_to_define_yet;
}
template_bank_cpu_t;

//prototypes of all functions that swig wraps to methods
void new_kfac_vec(
                  real_vector_single_cpu_t* vec,
                  const double kfac
                  );

void precondition_factor(
                         real_vector_single_cpu_t* vec
                         );

void compute_template_phasing(
                              complex_vector_single_cpu_t* exp_psi,
                              double M,
                              double eta,
                              int order,
                              double f_min,
                              double f_max,
                              real_vector_single_cpu_t* minus_one_by_three
                              );

#endif /* TEMPLATEBANKCPU_H */
