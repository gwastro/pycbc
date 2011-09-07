/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) matchedfiltercpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"
#include "../../../clayer/cpu/pycbccpu_types.h"
#include "matchedfiltercpu_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

void gen_snr_cpu(cpu_context_t* context,
		 real_vector_single_cpu_t* snr,
		 complex_vector_single_cpu_t* stilde, 
		 complex_vector_single_cpu_t* htilde,
		 complex_vector_single_cpu_t* q,
		 complex_vector_single_cpu_t* qtilde,
		 /*complex_fft_plan_t* plan,*/
		 double f_min,
		 double sigma_sq);

