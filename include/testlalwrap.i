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

/*

This module creates several functions to be used in the unittesting of the
PyCBC custom LAL wrapping, as documented in include/lalwrap.i and
include/pycbc_laltypes.i.

*/


%module testlalwrap

%include "pycbc_laltypes.i"

%unignore(invec1);
%unignore(invec2);
%unignore(ints1);
%unignore(ints2);
%unignore(infs1);
%unignore(infs2);

%unignore(invec);
%unignore(outvec);
%unignore(ints);
%unignore(outts);
%unignore(infs);
%unignore(outfs);

%unignore(length);
%unignore(value);
%unignore(epoch);
%unignore(deltaT);
%unignore(deltaF);

// First, our input tests.  They take a single input
// and double everything. For time and frequency series
// they also check that all members not shadowed in
// pycbc.types are initially zeroed out.
%define MAKE_TEST_INPUTV(CTYPE,SWTYPE)
%inline %{
int TestInput ## CTYPE (CTYPE *invec1, CTYPE *invec2) {
  UINT4 i;

  if (invec1 == NULL) return 1;
  if (invec1->data == NULL) return 1;
  if (invec1->length == 0) return 1;

  if (invec2 == NULL) return 1;
  if (invec2->data == NULL) return 1;
  if (invec2->length == 0) return 1;

  for (i=0; i < invec1->length; i++){
    invec1->data[i] = 2.0*invec1->data[i];
  }
  return 0;
}
%}
%unignore(TestInput ## CTYPE);
%apply CTYPE *INPUT_ ## SWTYPE {CTYPE *};
int TestInput ## CTYPE(CTYPE *invec1, CTYPE *invec2);
%enddef

%define MAKE_TEST_INPUTTS(CTYPE,SWTYPE)
%inline %{
int TestInput ## CTYPE (CTYPE *ints1, CTYPE *ints2) {
  UINT4 i;

  if (ints1 == NULL) return 1;
  if (ints1->data == NULL) return 1;
  if (ints1->data->data == NULL) return 1;
  if (ints1->data->length == 0) return 1;

  for (i=0; i < LALNameLength; i++){
    if (ints1->name[i] != '\0') return 1;
  }

  if (ints1->f0 != 0.0) return 1;

  if ((ints1->sampleUnits).powerOfTen != 0) return 1;
  for (i=0; i < LALNumUnits; i++){
    if ( (ints1->sampleUnits).unitNumerator[i] != 0) return 1;
    if ( (ints1->sampleUnits).unitDenominatorMinusOne[i] != 0) return 1;
  }

  if (ints2 == NULL) return 1;
  if (ints2->data == NULL) return 1;
  if (ints2->data->data == NULL) return 1;
  if (ints2->data->length == 0) return 1;

  for (i=0; i < LALNameLength; i++){
    if (ints2->name[i] != '\0') return 1;
  }

  if (ints2->f0 != 0.0) return 1;

  if ((ints2->sampleUnits).powerOfTen != 0) return 1;
  for (i=0; i < LALNumUnits; i++){
    if ( (ints2->sampleUnits).unitNumerator[i] != 0) return 1;
    if ( (ints2->sampleUnits).unitDenominatorMinusOne[i] != 0) return 1;
  }

  for (i=0; i < ints1->data->length; i++){
    ints1->data->data[i] = 2.0 * ints1->data->data[i];
  }

  ints1->deltaT = 2.0 * ints1->deltaT;
  (ints1->epoch).gpsSeconds = 2 * (ints1->epoch).gpsSeconds;
  (ints1->epoch).gpsNanoSeconds = 2 * (ints1->epoch).gpsNanoSeconds;

  return 0;
}
%}
%unignore(TestInput ## CTYPE);
%apply CTYPE *INPUT_ ## SWTYPE {CTYPE *};
int TestInput ## CTYPE(CTYPE *ints1, CTYPE *ints2);
%enddef

%define MAKE_TEST_INPUTFS(CTYPE,SWTYPE)
%inline %{
int TestInput ## CTYPE (CTYPE *infs1, CTYPE *infs2) {
  UINT4 i;

  if (infs1 == NULL) return 1;
  if (infs1->data == NULL) return 1;
  if (infs1->data->data == NULL) return 1;
  if (infs1->data->length == 0) return 1;

  for (i=0; i < LALNameLength; i++){
    if (infs1->name[i] != '\0') return 1;
  }

  if (infs1->f0 != 0.0) return 1;

  if ((infs1->sampleUnits).powerOfTen != 0) return 1;
  for (i=0; i < LALNumUnits; i++){
    if ( (infs1->sampleUnits).unitNumerator[i] != 0) return 1;
    if ( (infs1->sampleUnits).unitDenominatorMinusOne[i] != 0) return 1;
  }

  if (infs2 == NULL) return 1;
  if (infs2->data == NULL) return 1;
  if (infs2->data->data == NULL) return 1;
  if (infs2->data->length == 0) return 1;

  for (i=0; i < LALNameLength; i++){
    if (infs2->name[i] != '\0') return 1;
  }

  if (infs2->f0 != 0.0) return 1;

  if ((infs2->sampleUnits).powerOfTen != 0) return 1;
  for (i=0; i < LALNumUnits; i++){
    if ( (infs2->sampleUnits).unitNumerator[i] != 0) return 1;
    if ( (infs2->sampleUnits).unitDenominatorMinusOne[i] != 0) return 1;
  }

  for (i=0; i < infs1->data->length; i++){
    infs1->data->data[i] = 2.0 * infs1->data->data[i];
  }

  infs1->deltaF = 2.0 * infs1->deltaF;
  (infs1->epoch).gpsSeconds = 2 * (infs1->epoch).gpsSeconds;
  (infs1->epoch).gpsNanoSeconds = 2 * (infs1->epoch).gpsNanoSeconds;

  return 0;
}
%}
%unignore(TestInput ## CTYPE);
%apply CTYPE *INPUT_ ## SWTYPE {CTYPE *};
int TestInput ## CTYPE(CTYPE *infs1, CTYPE *infs2);
%enddef

// Next, NONEOUT tests.  They take two inputs: on return the
// second is the doubled version of the first.  But they also
// return the second (which should be mapped to 'None' by the
// typemap).

%define MAKE_TEST_NONEOUTV(CTYPE,SWTYPE)
%inline %{
  CTYPE *TestNoneout ## CTYPE (CTYPE *invec, CTYPE *outvec) {
  UINT4 i;

  if (invec == NULL) return NULL;
  if (invec->data == NULL) return NULL;

  if (outvec == NULL) return NULL;
  if (outvec->data == NULL) return NULL;

  if (invec->length != outvec->length) return NULL;
  if (invec->length == 0) return NULL;

  for (i=0; i < invec->length; i++){
    outvec->data[i] = 2.0*invec->data[i];
  }
  return outvec;
}
%}
%unignore(TestNoneout ## CTYPE);
%apply CTYPE *INPUT_ ## SWTYPE {CTYPE *};
%apply CTYPE *NONEOUT_ ## SWTYPE {CTYPE *TestNoneout ## CTYPE};
CTYPE *TestNoneout ## CTYPE(CTYPE *invec, CTYPE *outvec);
%enddef

%define MAKE_TEST_NONEOUTTS(CTYPE,SWTYPE)
%inline %{
  CTYPE *TestNoneout ## CTYPE (CTYPE *ints, CTYPE *outts) {
  UINT4 i;

  if (ints == NULL) return NULL;
  if (ints->data == NULL) return NULL;
  if (ints->data->data == NULL) return NULL;
  if (ints->data->length == 0) return NULL;

  if (outts == NULL) return NULL;
  if (outts->data == NULL) return NULL;
  if (outts->data->data == NULL) return NULL;
  if (outts->data->length == 0) return NULL;

  if (outts->data->length != ints->data->length) return NULL;

  for (i=0; i < LALNameLength; i++){
    if (ints->name[i] != '\0') return NULL;
    if (outts->name[i] != '\0') return NULL;
  }

  if (ints->f0 != 0.0) return NULL;
  if (outts->f0 != 0.0) return NULL;

  if ((ints->sampleUnits).powerOfTen != 0) return NULL;
  if ((outts->sampleUnits).powerOfTen != 0) return NULL;
  for (i=0; i < LALNumUnits; i++){
    if ( (ints->sampleUnits).unitNumerator[i] != 0) return NULL;
    if ( (ints->sampleUnits).unitDenominatorMinusOne[i] != 0) return NULL;
    if ( (outts->sampleUnits).unitNumerator[i] != 0) return NULL;
    if ( (outts->sampleUnits).unitDenominatorMinusOne[i] != 0) return NULL;
  }

  for (i=0; i < ints->data->length; i++){
    outts->data->data[i] = 2.0 * ints->data->data[i];
  }

  outts->deltaT = 2.0 * ints->deltaT;
  (outts->epoch).gpsSeconds = 2 * (ints->epoch).gpsSeconds;
  (outts->epoch).gpsNanoSeconds = 2 * (ints->epoch).gpsNanoSeconds;

  return outts;
}
%}
%unignore(TestNoneout ## CTYPE);
%apply CTYPE *INPUT_ ## SWTYPE {CTYPE *};
%apply CTYPE *NONEOUT_ ## SWTYPE {CTYPE *TestNoneout ## CTYPE};
CTYPE *TestNoneout ## CTYPE(CTYPE *ints, CTYPE *outts);
%enddef

%define MAKE_TEST_NONEOUTFS(CTYPE,SWTYPE)
%inline %{
  CTYPE *TestNoneout ## CTYPE (CTYPE *infs, CTYPE *outfs) {
  UINT4 i;

  if (infs == NULL) return NULL;
  if (infs->data == NULL) return NULL;
  if (infs->data->data == NULL) return NULL;
  if (infs->data->length == 0) return NULL;

  if (outfs == NULL) return NULL;
  if (outfs->data == NULL) return NULL;
  if (outfs->data->data == NULL) return NULL;
  if (outfs->data->length == 0) return NULL;

  if (outfs->data->length != infs->data->length) return NULL;

  for (i=0; i < LALNameLength; i++){
    if (infs->name[i] != '\0') return NULL;
    if (outfs->name[i] != '\0') return NULL;
  }

  if (infs->f0 != 0.0) return NULL;
  if (outfs->f0 != 0.0) return NULL;

  if ((infs->sampleUnits).powerOfTen != 0) return NULL;
  if ((outfs->sampleUnits).powerOfTen != 0) return NULL;
  for (i=0; i < LALNumUnits; i++){
    if ( (infs->sampleUnits).unitNumerator[i] != 0) return NULL;
    if ( (infs->sampleUnits).unitDenominatorMinusOne[i] != 0) return NULL;
    if ( (outfs->sampleUnits).unitNumerator[i] != 0) return NULL;
    if ( (outfs->sampleUnits).unitDenominatorMinusOne[i] != 0) return NULL;
  }

  for (i=0; i < infs->data->length; i++){
    outfs->data->data[i] = 2.0 * infs->data->data[i];
  }

  outfs->deltaF = 2.0 * infs->deltaF;
  (outfs->epoch).gpsSeconds = 2 * (infs->epoch).gpsSeconds;
  (outfs->epoch).gpsNanoSeconds = 2 * (infs->epoch).gpsNanoSeconds;

  return outfs;
}
%}
%unignore(TestNoneout ## CTYPE);
%apply CTYPE *INPUT_ ## SWTYPE {CTYPE *};
%apply CTYPE *NONEOUT_ ## SWTYPE {CTYPE *TestNoneout ## CTYPE};
CTYPE *TestNoneout ## CTYPE(CTYPE *infs, CTYPE *outfs);
%enddef

// Next, NEWOUT tests.  They take inputs specifying
// what should be constructed.

%define MAKE_TEST_NEWOUTV(CTYPE,SWTYPE,DTYPE)
%inline %{
CTYPE *TestNewout ## CTYPE (UINT4 length, DTYPE value) {
  UINT4 i;
  CTYPE *retvec;

  if (length == 0) return NULL;
  retvec = XLALCreate ## CTYPE (length);

  if (retvec == NULL) return NULL;
  for (i=0; i < retvec->length; i++){
    retvec->data[i] = value;
  }
  return retvec;
}
%}
%unignore(TestNewout ## CTYPE);
%apply CTYPE *NEWOUT_ ## SWTYPE {CTYPE *TestNewout ## CTYPE};
CTYPE *TestNewout ## CTYPE(UINT4 length, DTYPE value);
%enddef

%define MAKE_TEST_NEWOUTTS(CTYPE,SWTYPE,VTYPE,DTYPE)
%inline %{
CTYPE *TestNewout ## CTYPE (UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaT) {
  UINT4 i;
  CTYPE *retts;

  if (length == 0) return NULL;
  retts = calloc(1,sizeof(CTYPE));
  if (retts == NULL) return NULL;
  retts->data = XLALCreate ## VTYPE (length);

  if (retts->data == NULL) return NULL;
  for (i=0; i < retts->data->length; i++){
    retts->data->data[i] = value;
  }

  retts->deltaT = deltaT;
  (retts->epoch).gpsSeconds = epoch.gpsSeconds;
  (retts->epoch).gpsNanoSeconds = epoch.gpsNanoSeconds;

  return retts;
}
%}
%unignore(TestNewout ## CTYPE);
%apply CTYPE *NEWOUT_ ## SWTYPE {CTYPE *TestNewout ## CTYPE};
CTYPE *TestNewout ## CTYPE(UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaT);
%enddef

%define MAKE_TEST_NEWOUTFS(CTYPE,SWTYPE,VTYPE,DTYPE)
%inline %{
CTYPE *TestNewout ## CTYPE (UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaF) {
  UINT4 i;
  CTYPE *retfs;

  if (length == 0) return NULL;
  retfs = calloc(1,sizeof(CTYPE));
  if (retfs == NULL) return NULL;
  retfs->data = XLALCreate ## VTYPE (length);

  if (retfs->data == NULL) return NULL;
  for (i=0; i < retfs->data->length; i++){
    retfs->data->data[i] = value;
  }

  retfs->deltaF = deltaF;
  (retfs->epoch).gpsSeconds = epoch.gpsSeconds;
  (retfs->epoch).gpsNanoSeconds = epoch.gpsNanoSeconds;

  return retfs;
}
%}
%unignore(TestNewout ## CTYPE);
%apply CTYPE *NEWOUT_ ## SWTYPE {CTYPE *TestNewout ## CTYPE};
CTYPE *TestNewout ## CTYPE(UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaF);
%enddef

// Finally, ARGOUT tests.  They also take inputs specifying
// what should be constructed.

%define MAKE_TEST_ARGOUTV(CTYPE,SWTYPE,DTYPE)
%inline %{
int TestArgout ## CTYPE (CTYPE **outvec, UINT4 length, DTYPE value) {
  UINT4 i;

  if (*outvec != NULL) return 1;
  if (length == 0) return 1;
  *outvec = XLALCreate ## CTYPE (length);

  if (*outvec == NULL) return 1;
  for (i=0; i < (*outvec)->length; i++){
    (*outvec)->data[i] = value;
  }
  return 0;
}
%}
%unignore(TestArgout ## CTYPE);
%apply CTYPE **ARGOUT_ ## SWTYPE {CTYPE **};
int TestArgout ## CTYPE(CTYPE **outvec, UINT4 length, DTYPE value);
%enddef

%define MAKE_TEST_ARGOUTTS(CTYPE,SWTYPE,VTYPE,DTYPE)
%inline %{
int TestArgout ## CTYPE (CTYPE **outts, UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaT) {
  UINT4 i;

  if (*outts != NULL) return 1;
  if (length == 0) return 1;
  *outts = calloc(1,sizeof(CTYPE));
  if (*outts == NULL) return 1;
  (*outts)->data = XLALCreate ## VTYPE (length);

  if ((*outts)->data == NULL) return 1;
  for (i=0; i < (*outts)->data->length; i++){
    (*outts)->data->data[i] = value;
  }

  (*outts)->deltaT = deltaT;
  ((*outts)->epoch).gpsSeconds = epoch.gpsSeconds;
  ((*outts)->epoch).gpsNanoSeconds = epoch.gpsNanoSeconds;

  return 0;
}
%}
%unignore(TestArgout ## CTYPE);
%apply CTYPE **ARGOUT_ ## SWTYPE {CTYPE **};
int TestArgout ## CTYPE(CTYPE **outts, UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaT);
%enddef

%define MAKE_TEST_ARGOUTFS(CTYPE,SWTYPE,VTYPE,DTYPE)
%inline %{
int TestArgout ## CTYPE (CTYPE **outfs, UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaF) {
  UINT4 i;

  if (*outfs != NULL) return 1;
  if (length == 0) return 1;
  *outfs = calloc(1,sizeof(CTYPE));
  if (*outfs == NULL) return 1;
  (*outfs)->data = XLALCreate ## VTYPE (length);

  if ((*outfs)->data == NULL) return 1;
  for (i=0; i < (*outfs)->data->length; i++){
    (*outfs)->data->data[i] = value;
  }

  (*outfs)->deltaF = deltaF;
  ((*outfs)->epoch).gpsSeconds = epoch.gpsSeconds;
  ((*outfs)->epoch).gpsNanoSeconds = epoch.gpsNanoSeconds;

  return 0;
}
%}
%unignore(TestArgout ## CTYPE);
%apply CTYPE **ARGOUT_ ## SWTYPE {CTYPE **};
int TestArgout ## CTYPE(CTYPE **outfs, UINT4 length, DTYPE value, LIGOTimeGPS epoch, REAL8 deltaF);
%enddef

// Now actually create the wrapped functions

// Vector test functions
MAKE_TEST_INPUTV(REAL4Vector,REAL4V)
MAKE_TEST_NONEOUTV(REAL4Vector,REAL4V)
MAKE_TEST_NEWOUTV(REAL4Vector,REAL4V,REAL4)
MAKE_TEST_ARGOUTV(REAL4Vector,REAL4V,REAL4)

MAKE_TEST_INPUTV(REAL8Vector,REAL8V)
MAKE_TEST_NONEOUTV(REAL8Vector,REAL8V)
MAKE_TEST_NEWOUTV(REAL8Vector,REAL8V,REAL8)
MAKE_TEST_ARGOUTV(REAL8Vector,REAL8V,REAL8)

MAKE_TEST_INPUTV(COMPLEX8Vector,COMPLEX8V)
MAKE_TEST_NONEOUTV(COMPLEX8Vector,COMPLEX8V)
MAKE_TEST_NEWOUTV(COMPLEX8Vector,COMPLEX8V,COMPLEX8)
MAKE_TEST_ARGOUTV(COMPLEX8Vector,COMPLEX8V,COMPLEX8)

MAKE_TEST_INPUTV(COMPLEX16Vector,COMPLEX16V)
MAKE_TEST_NONEOUTV(COMPLEX16Vector,COMPLEX16V)
MAKE_TEST_NEWOUTV(COMPLEX16Vector,COMPLEX16V,COMPLEX16)
MAKE_TEST_ARGOUTV(COMPLEX16Vector,COMPLEX16V,COMPLEX16)

// Time series test functions
MAKE_TEST_INPUTTS(REAL4TimeSeries,REAL4TS)
MAKE_TEST_NONEOUTTS(REAL4TimeSeries,REAL4TS)
MAKE_TEST_NEWOUTTS(REAL4TimeSeries,REAL4TS,REAL4Vector,REAL4)
MAKE_TEST_ARGOUTTS(REAL4TimeSeries,REAL4TS,REAL4Vector,REAL4)

MAKE_TEST_INPUTTS(REAL8TimeSeries,REAL8TS)
MAKE_TEST_NONEOUTTS(REAL8TimeSeries,REAL8TS)
MAKE_TEST_NEWOUTTS(REAL8TimeSeries,REAL8TS,REAL8Vector,REAL8)
MAKE_TEST_ARGOUTTS(REAL8TimeSeries,REAL8TS,REAL8Vector,REAL8)

MAKE_TEST_INPUTTS(COMPLEX8TimeSeries,COMPLEX8TS)
MAKE_TEST_NONEOUTTS(COMPLEX8TimeSeries,COMPLEX8TS)
MAKE_TEST_NEWOUTTS(COMPLEX8TimeSeries,COMPLEX8TS,COMPLEX8Vector,COMPLEX8)
MAKE_TEST_ARGOUTTS(COMPLEX8TimeSeries,COMPLEX8TS,COMPLEX8Vector,COMPLEX8)

MAKE_TEST_INPUTTS(COMPLEX16TimeSeries,COMPLEX16TS)
MAKE_TEST_NONEOUTTS(COMPLEX16TimeSeries,COMPLEX16TS)
MAKE_TEST_NEWOUTTS(COMPLEX16TimeSeries,COMPLEX16TS,COMPLEX16Vector,COMPLEX16)
MAKE_TEST_ARGOUTTS(COMPLEX16TimeSeries,COMPLEX16TS,COMPLEX16Vector,COMPLEX16)

// Frequency series test functions
MAKE_TEST_INPUTFS(REAL4FrequencySeries,REAL4FS)
MAKE_TEST_NONEOUTFS(REAL4FrequencySeries,REAL4FS)
MAKE_TEST_NEWOUTFS(REAL4FrequencySeries,REAL4FS,REAL4Vector,REAL4)
MAKE_TEST_ARGOUTFS(REAL4FrequencySeries,REAL4FS,REAL4Vector,REAL4)

MAKE_TEST_INPUTFS(REAL8FrequencySeries,REAL8FS)
MAKE_TEST_NONEOUTFS(REAL8FrequencySeries,REAL8FS)
MAKE_TEST_NEWOUTFS(REAL8FrequencySeries,REAL8FS,REAL8Vector,REAL8)
MAKE_TEST_ARGOUTFS(REAL8FrequencySeries,REAL8FS,REAL8Vector,REAL8)

MAKE_TEST_INPUTFS(COMPLEX8FrequencySeries,COMPLEX8FS)
MAKE_TEST_NONEOUTFS(COMPLEX8FrequencySeries,COMPLEX8FS)
MAKE_TEST_NEWOUTFS(COMPLEX8FrequencySeries,COMPLEX8FS,COMPLEX8Vector,COMPLEX8)
MAKE_TEST_ARGOUTFS(COMPLEX8FrequencySeries,COMPLEX8FS,COMPLEX8Vector,COMPLEX8)

MAKE_TEST_INPUTFS(COMPLEX16FrequencySeries,COMPLEX16FS)
MAKE_TEST_NONEOUTFS(COMPLEX16FrequencySeries,COMPLEX16FS)
MAKE_TEST_NEWOUTFS(COMPLEX16FrequencySeries,COMPLEX16FS,COMPLEX16Vector,COMPLEX16)
MAKE_TEST_ARGOUTFS(COMPLEX16FrequencySeries,COMPLEX16FS,COMPLEX16Vector,COMPLEX16)
