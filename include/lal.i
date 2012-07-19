/* **********************************************************************
Copyright (C) 2012  Josh Willis

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*************************************************************************/


%module lal

%{
#include <lal/VectorOps.h>
#include <lal/ComplexFFT.h>
%}

%ignore XLALCOMPLEX8VectorFFT;
%ignore XLALSSVectorMultiply;
%ignore XLALCreateREAL4Vector;

%include "laltypemaps.i"

// Some tests:
%rename("%s") XLALSSVectorMultiply;
%apply REAL4Vector *NONEOUT_REAL4V { REAL4Vector *XLALSSVectorMultiply };
%apply REAL4Vector *INPUT_REAL4V {REAL4Vector *};
extern REAL4Vector *XLALSSVectorMultiply(REAL4Vector *out, REAL4Vector *in1, REAL4Vector *in2);

%rename("%s") XLALCreateREAL4Vector;
%newobject XLALCreateReal4Vector;
%apply REAL4Vector *NEWOUT_REAL4V {REAL4Vector *XLALCreateREAL4Vector};
extern REAL4Vector *XLALCreateREAL4Vector(UINT4 length);


%rename("%s") XLALCOMPLEX8VectorFFT;
%apply COMPLEX8Vector *INPUT_COMPLEX8V {COMPLEX8Vector *};
extern int XLALCOMPLEX8VectorFFT(COMPLEX8Vector *output, COMPLEX8Vector *input, const COMPLEX8FFTPlan *plan );

