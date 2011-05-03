#include <stdio.h>
#include "datavectorcpu_types.h"
#include "datavectorcpu_prototypes.h"

/* real vector manipulation */
real_vector_t* new_real_vector_t(int length, int memory_location )
{
    
    real_vector_t* c;
    c = (real_vector_t*) malloc( sizeof(real_vector_t) );
    
    c->t_start = 0;
    c->dx = 1;
    c->vector_length = length;
    c->memory_type = memory_location;
    
    c->element_size_bytes = sizeof(float);
    c->data = malloc( length * c->element_size_bytes );
    
    printf("created real_vector_t at %p\n", c );
    return c;
}

void delete_real_vector_t( real_vector_t* p )
{
    printf("deleting real_vector_t at %p\n", p );
    free( p->data );
    free( p );
}
