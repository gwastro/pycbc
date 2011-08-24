/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) templatebankcpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"
#include "../../../clayer/cpu/pycbccpu_types.h"
#include "templatebankcpu_prototypes.h"
%}


