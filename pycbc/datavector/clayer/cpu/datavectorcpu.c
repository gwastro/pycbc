// Copyright (C) 2011 Karsten Wiesner, Josh Willis
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
// datavector constructors and destructors implementation for pycbc


#include <stdio.h>
#include <stdlib.h> /* Needed for posix_memalign decl */
#include <pycbc/clayer/except.h>
/* Following needed for memset() */
#include <string.h>

#include "pycbccpu.h"
#include "datavectorcpu.h"
#include "datavectorcpu_private.h"

/*
  The constant below is used by posix_memalign() to ensure that the
  datavectors generated for the CPU architecture are properly aligned
  to allow the use of FFTW SIMD extensions (a capability assumed in the
  FFTW wrappings of PyCBC).

  It is defined here so that the CPU datavectors can be built even if
  FFTW is not being used as the FFT library (though the alignment may,
  in that case, be unnecessary).

  Also note that maintainers should keep up enough with the internals
  of FFTW to know the appropriate value to use here.  At present, LAL
  and the LSC use FFTW 3.2, for which 16-byte alignment is all that is
  needed.  However, for FFTW 3.3 (which supports AVX instructions), the
  documentation states that still only 16-byte alignment is -needed-
  (which is true) but internally it used 32-byte alignment because that
  can give a performance benefit (even though not needed to avoid
  segmentation faults).

 */
#define PYCBC_CPU_DATAVECTOR_ALIGNMENT 16

real_vector_single_cpu_t* new_real_vector_single_cpu_t(cpu_context_t* context,
                                                       unsigned long length,
                                                       double delta_x)
{
    CONSTRUCTOR_TEMPLATE(real_vector_single_cpu_t, float)

    if (posix_memalign((void **) &(c->data),PYCBC_CPU_DATAVECTOR_ALIGNMENT,
		       c->meta_data.vector_length*c->meta_data.element_size_bytes))
	{ /*
	    We failed, for some reason.  There should not be any memory allocated to
	    c->data, but we must clean up c itself.
	  */
	  free(c);
	  pycbc_throw_exception(PYCBC_MEMORY_ERROR,"real_vector_single_cpu_t allocation failed");
	  return NULL;
	}
    else
      {
	memset((void *) c->data,0,
	       c->meta_data.vector_length*c->meta_data.element_size_bytes);
	return c;
      }
}

void delete_real_vector_single_cpu_t( real_vector_single_cpu_t* p )
{
    free( p->data );
    free( p );
}

real_vector_double_cpu_t* new_real_vector_double_cpu_t(cpu_context_t* context,
                                                       unsigned long length,
                                                       double delta_x)
{
    CONSTRUCTOR_TEMPLATE(real_vector_double_cpu_t, double)

    if (posix_memalign((void **) &(c->data),PYCBC_CPU_DATAVECTOR_ALIGNMENT,
		       c->meta_data.vector_length*c->meta_data.element_size_bytes))
	{/*
	    We failed, for some reason.  There should not be any memory allocated to
	    c->data, but we must clean up c itself.
	  */
	  free(c);
	  pycbc_throw_exception(PYCBC_MEMORY_ERROR,"real_vector_double_cpu_t allocation failed");
	  return NULL;
	}
    else
      {
	memset((void *) c->data,0,
	       c->meta_data.vector_length*c->meta_data.element_size_bytes);
	return c;
      }
}

void delete_real_vector_double_cpu_t( real_vector_double_cpu_t* p )
{
    free( p->data );
    free( p );
}

complex_vector_single_cpu_t* new_complex_vector_single_cpu_t(cpu_context_t* context,
                                                             unsigned long length,
                                                             double delta_x)
{
    CONSTRUCTOR_TEMPLATE(complex_vector_single_cpu_t, complex float)

    if (posix_memalign((void **) &(c->data),PYCBC_CPU_DATAVECTOR_ALIGNMENT,
		       c->meta_data.vector_length*c->meta_data.element_size_bytes))
	{/*
	    We failed, for some reason.  There should not be any memory allocated to
	    c->data, but we must clean up c itself.
	  */
	  free(c);
	  pycbc_throw_exception(PYCBC_MEMORY_ERROR,"complex_vector_single_cpu_t allocation failed");
	  return NULL;
	}
    else
      {
	memset((void *) c->data,0,
	       c->meta_data.vector_length*c->meta_data.element_size_bytes);
	return c;
      }
}

void delete_complex_vector_single_cpu_t( complex_vector_single_cpu_t* p )
{

    free( p->data );
    free( p );
}

complex_vector_double_cpu_t* new_complex_vector_double_cpu_t(cpu_context_t* context,
                                                             unsigned long length,
                                                             double delta_x)
{

    CONSTRUCTOR_TEMPLATE(complex_vector_double_cpu_t, complex double)

    if (posix_memalign((void **) &(c->data),PYCBC_CPU_DATAVECTOR_ALIGNMENT,
		       c->meta_data.vector_length*c->meta_data.element_size_bytes))
	{/*
	    We failed, for some reason.  There should not be any memory allocated to
	    c->data, but we must clean up c itself.
	  */
	  free(c);
	  pycbc_throw_exception(PYCBC_MEMORY_ERROR,"complex_vector_double_cpu_t allocation failed");
	  return NULL;
	}
    else
      {
	memset((void *) c->data,0,
	       c->meta_data.vector_length*c->meta_data.element_size_bytes);
	return c;
      }
}

void delete_complex_vector_double_cpu_t( complex_vector_double_cpu_t* p )
{
    free( p->data );
    free( p );
}

