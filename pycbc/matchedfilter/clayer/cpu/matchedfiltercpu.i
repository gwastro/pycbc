// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:syntax=c

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING,module="pycbc.matchedfilter.clayer") cpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "datavectorcpu.h"
#include "pycbccpu.h"    
#include "matchedfiltercpu.h"
#include "matchedfiltercpu_private.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

%include "matchedfiltercpu.h"

