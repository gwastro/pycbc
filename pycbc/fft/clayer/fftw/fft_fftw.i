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
// fft swig file for pycbc


%define DOCSTRING
"Copyright 2011, Josh Willis <joshua.willis@ligo.org>."
%enddef

%module(docstring=DOCSTRING) fft_fftw

%feature("autodoc", "1");

// This goes directly to the wrap-code (no swig preprocess)
%{
#include "fft_fftw.h"
#include "fft_fftw_private.h"
%}

%include "fft_fftw.h"

%extend fft_real_single_plan_fftw {
	fft_real_single_plan_fftw(unsigned long size, int fwdflag, int measurelvl);
	~fft_real_single_plan_fftw();
}

%extend fft_complex_single_plan_fftw {
	fft_complex_single_plan_fftw(unsigned long size, int fwdflag, int measurelvl);
	~fft_complex_single_plan_fftw();
}

%extend fft_real_double_plan_fftw {
	fft_real_double_plan_fftw(unsigned long size, int fwdflag, int measurelvl);
	~fft_real_double_plan_fftw();
}

%extend fft_complex_double_plan_fftw {
	fft_complex_double_plan_fftw(unsigned long size, int fwdflag, int measurelvl);
	~fft_complex_double_plan_fftw();
}

