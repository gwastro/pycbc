/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef
%module(docstring=DOCSTRING) datavectorcpu
%feature("autodoc", "1");
%{
#include "datavectorcpu.h"
    %}

%pythoncode %{
    import types
    import warnings
    %}

/* Parse the header file to generate wrappers */
%include "datavectorcpu.h"

%extend real_vector_t {
    real_vector_t(int vector_length, int on_gpu);
    ~real_vector_t();
    
    char *__str__() {
        static char a[1024];
        snprintf( a, sizeof(a)/sizeof(*a), 
                     "<real_vector_t, in cpu memory, length %d, step %e>", 
                     self->vector_length, self->dx );
        return a;
    }
    
    int __len__() {
        return self->vector_length;
    }
    
    double __getitem__(int i) {
        float* data = (float*) self->data; 
        return (float) data[i];
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
