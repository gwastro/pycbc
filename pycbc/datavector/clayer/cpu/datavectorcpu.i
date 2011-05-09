/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) datavectorcpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "datavectorcpu_types.h"
#include "datavectorcpu_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

// To extend the c-types by methodes they have to be defined here
// but to declare function prototypes as well would raise a 
// "is multiply defined error". That is the reason for splitting 
// the headerfiles
%include "datavectorcpu_types.h"

%extend real_vector_single_t {
    real_vector_single_t(int vector_length, int on_gpu);
    ~real_vector_single_t();
    
    char *__str__() {
        static char a[1024];
        snprintf( a, sizeof(a)/sizeof(*a), 
                     "<real_vector_single_t, in cpu memory, length %d, data ptr %p>", 
                     self->meta_data.vector_length, self->data );
        return a;
    }
    
    int __len__() {
        return self->meta_data.vector_length;
    }
    
    
    // 
    // TBD a zero init function which is called automatically in the constructor
    // but could also be called from outside
    //
    
    double __getitem__(int i) {
        float* data = (float*) self->data; 
        return (float) data[i];
    }
    
    void set_t_start( unsigned long int t_start ) {
        self->meta_data.t_start = t_start;
    }
    
    double get_t_start( void ) {
        return self->meta_data.t_start;
    }
    
    void set_dx( double dx ) {
        self->meta_data.dx = dx;
    }
    
    double get_dx( void ) {
        return self->meta_data.dx;
    }
}
