/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2006, 2011 Duncan Brown <dabrown@physics.syr.edu>."
%enddef
%module(docstring=DOCSTRING) pycbc
%feature("autodoc", "1");
%{
#include "pycbc.h"
%}

%pythoncode %{
import types
import warnings
%}

/* Parse the header file to generate wrappers */
%include "pycbc.h"

%extend real_vector_t {
  real_vector_t(int vector_length, int on_gpu);
  ~real_vector_t();

  char *__str__() {
    static char a[1024];
    if ( self->memory_type == gpu_cuda_global_memory )
      snprintf( a, sizeof(a)/sizeof(*a), 
          "<real_vector_t in cuda global memory, length %d, step %e>", 
          self->vector_length, self->dx );
    else if ( self->memory_type == cpu_generic_memory )
      snprintf( a, sizeof(a)/sizeof(*a), 
          "<real_vector_t, in cpu memory, length %d, step %e>", 
          self->vector_length, self->dx );
    return a;
  }

  int __len__() {
    return self->vector_length;
  }

  double __getitem__(int i) {
    if ( self->memory_type == cpu_generic_memory )
    {
      REAL4* data = (REAL4*) self->data; 
      return (double) data[i];
    }
    else return 0;
  }

  void set_t_start( unsigned long int t_start ) {
    self->t_start = t_start;
  }

  double get_t_start( void ) {
    return self->t_start;
  }

  void set_dx( double dx ) {
    self->dx = dx;
  }

  double get_dx( void ) {
    return self->dx;
  }
}

%extend complex_vector_t {
  complex_vector_t(int n);
  ~complex_vector_t();

  char *__str__() {
    static char a[1024];
    if ( self->memory_type == gpu_cuda_global_memory )
      snprintf( a, sizeof(a)/sizeof(*a), 
          "<complex_vector_t in cuda global memory, length %d, step %e>", 
          self->vector_length, self->dx );
    else if ( self->memory_type == cpu_generic_memory )
      snprintf( a, sizeof(a)/sizeof(*a), 
          "<complex_vector_t, in cpu memory, length %d, step %e>", 
          self->vector_length, self->dx );
  }

  int __len__() {
    return self->vector_length;
  }

  complex __getitem__(int i) {
    if ( self->memory_type == cpu_generic_memory )
    {
      COMPLEX8* data = (COMPLEX8*) self->data; 
      return (complex) data[i];
    }
    else return 0;
  }

  void set_t_start( unsigned long int t_start ) {
    self->t_start = t_start;
  }

  double get_t_start( void ) {
    return self->t_start;
  }

  void set_dx( double dx ) {
    self->dx = dx;
  }

  double get_dx( void ) {
    return self->dx;
  }
}
