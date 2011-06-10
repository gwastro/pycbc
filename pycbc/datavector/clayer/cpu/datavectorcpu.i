// Copyright (C) 2011 Karsten Wiesner
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
// datavector swig file for pycbc

/* vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab */

%define DOCSTRING
"Copyright 2011, 2011 Karsten Wiesner <karsten.wiesner@ligo.org>."
%enddef

%module(docstring=DOCSTRING) datavectorcpu

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
// wrap code needs to have typedefs and function prototypes!
%{
#include "../except.h"
#include "datavectorcpu_types.h"
#include "datavectorcpu_prototypes.h"
%}

%pythoncode %{
    import types
    import warnings
    %}

%include "exception.i"

// To extend the c-types by methodes they have to be defined here
// but to declare function prototypes as well would raise a 
// "is multiply defined error". That is the reason for splitting 
// the headerfiles into _types and _prototypes
%include "datavectorcpu_types.h"
%include "../datavector_types.i"

%extend real_vector_single_t {
    TYPE_INTERFACE_TEMPLATE(real_vector_single_t,float)
}

%extend real_vector_double_t {
    TYPE_INTERFACE_TEMPLATE(real_vector_double_t,double)
}

%extend complex_vector_single_t {
    TYPE_INTERFACE_TEMPLATE(complex_vector_single_t,complex_float_t)
}

%extend complex_vector_double_t {
    TYPE_INTERFACE_TEMPLATE(complex_vector_double_t,complex_double_t)
}



