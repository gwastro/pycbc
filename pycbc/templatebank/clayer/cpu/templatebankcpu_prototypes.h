#ifndef TEMPLATEBANKCPU_PROTOTYPES_H
#define TEMPLATEBANKCPU_PROTOTYPES_H

#include <stdlib.h>
#include "../../../clayer/cpu/pycbccpu_types.h"

//prototypes of all functions that swig wraps to methods
void new_kfac_vec(
    real_vector_single_t* vec,
    const double kfac
    );

void precondition_factor(
    real_vector_single_t* vec
    );

void compute_template_phasing(
    complex_vector_single_t* exp_psi,
    double M,
    double eta,
    int order,
    double f_min,
    double f_max,
    real_vector_single_t* minus_one_by_three
    );

#endif /* TEMPLATEBANKCPU_PROTOTYPES_H */
