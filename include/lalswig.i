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


%module lalswig

%{
#include <lal/VectorOps.h>
#include <lal/ComplexFFT.h>
%}

%include "laltypemaps.i"

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

%unignore(XLALCreateREAL4TimeSeries);
%newobject XLALCreateREAL4TimeSeries;
%apply REAL4TimeSeries *NEWOUT_REAL4TS {REAL4TimeSeries *XLALCreateREAL4TimeSeries};
extern REAL4TimeSeries *XLALCreateREAL4TimeSeries(UINT4 length);

%unignore(XLALCreateREAL8TimeSeries);
%newobject XLALCreateREAL8TimeSeries;
%apply REAL8TimeSeries *NEWOUT_REAL8TS {REAL8TimeSeries *XLALCreateREAL8TimeSeries};
extern REAL8TimeSeries *XLALCreateREAL8TimeSeries(UINT4 length);

%unignore(XLALCreateCOMPLEX8TimeSeries);
%newobject XLALCreateCOMPLEX8TimeSeries;
%apply COMPLEX8TimeSeries *NEWOUT_COMPLEX8TS {COMPLEX8TimeSeries *XLALCreateCOMPLEX8TimeSeries};
extern COMPLEX8TimeSeries *XLALCreateCOMPLEX8TimeSeries(UINT4 length);

%unignore(XLALCreateCOMPLEX16TimeSeries);
%newobject XLALCreateCOMPLEX16TimeSeries;
%apply COMPLEX16TimeSeries *NEWOUT_COMPLEX16TS {COMPLEX16TimeSeries *XLALCreateCOMPLEX16TimeSeries};
extern COMPLEX16TimeSeries *XLALCreateCOMPLEX16TimeSeries(UINT4 length);

%unignore(XLALCreateREAL4FrequencySeries);
%newobject XLALCreateREAL4FrequencySeries;
%apply REAL4FrequencySeries *NEWOUT_REAL4FS {REAL4FrequencySeries *XLALCreateREAL4FrequencySeries};
extern REAL4FrequencySeries *XLALCreateREAL4FrequencySeries(UINT4 length);

%unignore(XLALCreateREAL8FrequencySeries);
%newobject XLALCreateREAL8FrequencySeries;
%apply REAL8FrequencySeries *NEWOUT_REAL8FS {REAL8FrequencySeries *XLALCreateREAL8FrequencySeries};
extern REAL8FrequencySeries *XLALCreateREAL8FrequencySeries(UINT4 length);

%unignore(XLALCreateCOMPLEX8FrequencySeries);
%newobject XLALCreateCOMPLEX8FrequencySeries;
%apply COMPLEX8FrequencySeries *NEWOUT_COMPLEX8FS {COMPLEX8FrequencySeries *XLALCreateCOMPLEX8FrequencySeries};
extern COMPLEX8FrequencySeries *XLALCreateCOMPLEX8FrequencySeries(UINT4 length);

%unignore(XLALCreateCOMPLEX16FrequencySeries);
%newobject XLALCreateCOMPLEX16FrequencySeries;
%apply COMPLEX16FrequencySeries *NEWOUT_COMPLEX16FS {COMPLEX16FrequencySeries *XLALCreateCOMPLEX16FrequencySeries};
extern COMPLEX16FrequencySeries *XLALCreateCOMPLEX16FrequencySeries(UINT4 length);


