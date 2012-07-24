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


%module lalwrap

%{
#include <lal/VectorOps.h>
#include <lal/ComplexFFT.h>
%}

%include "pycbc_laltypes.i"

// Some tests:

%unignore(XLALSSVectorMultiply);
%apply REAL4Vector *NONEOUT_REAL4V { REAL4Vector *XLALSSVectorMultiply };
%apply REAL4Vector *INPUT_REAL4V {REAL4Vector *};
extern REAL4Vector *XLALSSVectorMultiply(REAL4Vector *out, REAL4Vector *in1, REAL4Vector *in2);

%unignore(XLALCOMPLEX8VectorFFT);
%apply COMPLEX8Vector *INPUT_COMPLEX8V {COMPLEX8Vector *};
extern int XLALCOMPLEX8VectorFFT(COMPLEX8Vector *output, COMPLEX8Vector *input, const COMPLEX8FFTPlan *plan );

%unignore(XLALCreateREAL4Vector);
%newobject XLALCreateREAL4Vector;
%apply REAL4Vector *NEWOUT_REAL4V {REAL4Vector *XLALCreateREAL4Vector};
extern REAL4Vector *XLALCreateREAL4Vector(UINT4 length);

%unignore(XLALCreateREAL8Vector);
%newobject XLALCreateREAL8Vector;
%apply REAL8Vector *NEWOUT_REAL8V {REAL8Vector *XLALCreateREAL8Vector};
extern REAL8Vector *XLALCreateREAL8Vector(UINT4 length);

%unignore(XLALCreateCOMPLEX8Vector);
%newobject XLALCreateCOMPLEX8Vector;
%apply COMPLEX8Vector *NEWOUT_COMPLEX8V {COMPLEX8Vector *XLALCreateCOMPLEX8Vector};
extern COMPLEX8Vector *XLALCreateCOMPLEX8Vector(UINT4 length);

%unignore(XLALCreateCOMPLEX16Vector);
%newobject XLALCreateCOMPLEX16Vector;
%apply COMPLEX16Vector *NEWOUT_COMPLEX16V {COMPLEX16Vector *XLALCreateCOMPLEX16Vector};
extern COMPLEX16Vector *XLALCreateCOMPLEX16Vector(UINT4 length);
