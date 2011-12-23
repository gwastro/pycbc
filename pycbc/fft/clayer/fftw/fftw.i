// Copyright (C) 2011 Josh Willis
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
// FFTW swig file for pycbc


%define DOCSTRING
"Copyright 2011, Josh Willis <joshua.willis@ligo.org>."
%enddef

%module(docstring=DOCSTRING) fftw

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
%{
#include "fftw.h"
#include "fftw_private.h"
%}

%include <pycbc/clayer/except.i>
%include "fftw.h"

%extend fftw_real_single_plan {
	fftw_real_single_plan(unsigned long size, int fwdflag, int measurelvl);
	~fftw_real_single_plan();
}

%extend fftw_complex_single_plan {
	fftw_complex_single_plan(unsigned long size, int fwdflag, int measurelvl);
	~fftw_complex_single_plan();
}

%extend fftw_real_double_plan {
	fftw_real_double_plan(unsigned long size, int fwdflag, int measurelvl);
	~fftw_real_double_plan();
}

%extend fftw_complex_double_plan {
	fftw_complex_double_plan(unsigned long size, int fwdflag, int measurelvl);
	~fftw_complex_double_plan();
}

