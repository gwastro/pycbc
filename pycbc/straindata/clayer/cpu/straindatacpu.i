// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:syntax=c

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>, Josh Willis <joshua.willis@ligo.org>"
%enddef

%module(docstring=DOCSTRING,module="pycbc.straindata.clayer") straindatacpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "datavectorcpu.h"
#include "pycbccpu.h" 
#include "straindatacpu.h"
#include "straindatacpu_private.h"
#include <except.i>
%}

%include <except.i>

%pythoncode %{
    import types
    import warnings
    %}

%include "straindatacpu.h"
