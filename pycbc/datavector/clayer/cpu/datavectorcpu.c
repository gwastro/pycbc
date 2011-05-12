#include <stdio.h>
#include "datavectorcpu_types.h"
#include "datavectorcpu_prototypes.h"

real_vector_single_t* new_real_vector_single_t(int length)
{
    
    real_vector_single_t* c;
    c = (real_vector_single_t*) malloc( sizeof(real_vector_single_t) );
    
    c->meta_data.start = 0;
    c->meta_data.dx = 1;
    c->meta_data.vector_length = length;
    
    c->meta_data.element_size_bytes = sizeof(float); // single precision specified here
    c->data = calloc( length , c->meta_data.element_size_bytes );

    printf("created real_vector_single_t at %p\n", c );
    return c;
}

void delete_real_vector_single_t( real_vector_single_t* p )
{
    printf("deleting real_vector_single_t at %p\n", p );
    free( p->data );
    free( p );
}
