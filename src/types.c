#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <lal/LALDatatypes.h>
#include "pycbc.h"

/* real vector manipulation */
real_vector_t* new_real_vector_t( 
    int length, enum cbc_memory_meta_types_t memory_location )
{
  real_vector_t* c = (real_vector_t*) malloc( sizeof(real_vector_t) );

  c->t_start = 0;
  c->dx = 1;
  c->vector_length = length;
  c->memory_type = memory_location;

  if ( memory_location == gpu_cuda_global_memory )
  {
    c->element_size_bytes = sizeof(cufftReal);
    cudaMalloc( (void**) &(c->data), length * c->element_size_bytes );
  }
  else if ( memory_location == cpu_generic_memory )
  {
    c->element_size_bytes = sizeof(REAL4);
    c->data = malloc( length * c->element_size_bytes );
  }

#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "created real_vector_t at %p\n", c );
#endif
  return c;
}
  
void delete_real_vector_t( real_vector_t* p )
{
#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "deleting real_vector_t at %p\n", p );
#endif
    if ( p->memory_type == gpu_cuda_global_memory )
    {
      cudaFree( p->data );
    }
    else if ( p->memory_type == cpu_generic_memory )
    {
      free( p->data );
    }
  free( p );
}

/* complex vector manipulation */
complex_vector_t* new_complex_vector_t( 
    int length, enum cbc_memory_meta_types_t memory_location )
{
  complex_vector_t* c = (complex_vector_t*) malloc( sizeof(complex_vector_t) );

  c->t_start = 0;
  c->dx = 1;
  c->vector_length = length;
  c->memory_type = memory_location;

  if ( memory_location == gpu_cuda_global_memory )
  {
    c->element_size_bytes = sizeof(cufftComplex);
    cudaMalloc( (void**) &(c->data), length * c->element_size_bytes );
  }
  else if ( memory_location == cpu_generic_memory )
  {
    c->element_size_bytes = sizeof(COMPLEX8);
    c->data = malloc( length * c->element_size_bytes );
  }

#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "created complex_vector_t at %p\n", c );
#endif
  return c;
}
  
void delete_complex_vector_t( complex_vector_t* p )
{
#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "deleting complex_vector_t at %p\n", p );
#endif
    if ( p->memory_type == gpu_cuda_global_memory )
    {
      cudaFree( p->data );
    }
    else if ( p->memory_type == cpu_generic_memory )
    {
      free( p->data );
    }
  free( p );
}
