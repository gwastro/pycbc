#ifndef TEMPLATEBANKCPU_PROTOTYPES_H
#define TEMPLATEBANKCPU_PROTOTYPES_H

#include <stdlib.h>
#include "../../../clayer/cpu/pycbccpu_types.h"

//prototypes of all functions that swig wraps to methods
real_vector_single_t* new_kfac_vec(
    const unsigned length,
    const double kfac,
    const double deltax
    );

real_vector_single_t* precondition_factor(
    const unsigned length,
    const double deltat
    );

void compute_template_phasing(
    complex_vector_single_t* exp_psi,
    double m1,
    double m2,
    double beta,
    int order,
    double f_min,
    double f_max,
    int N,
    double dt,
    real_vector_single_t* minus_one_by_three
    );

#endif /* TEMPLATEBANKCPU_PROTOTYPES_H */
