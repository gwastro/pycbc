#ifndef DATAVECTORCPU_TYPES_H
#define DATAVECTORCPU_TYPES_H

#include <stdlib.h>
#include "../datavector_types.h"

typedef struct
{
    meta_data_t meta_data;
    float *data;
}
real_vector_single_t;

#endif /* DATAVECTORCPU_TYPES_H */
