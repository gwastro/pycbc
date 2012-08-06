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

PyCBC custom LAL wrappings

Use this file to wrap time-critical functions in LAL that should be accessed
in PyCBC, and which accept or return a LAL Vector, TimeSeries, or
FrequencySeries.  To add a function do the following:

(0) Ensure that it is a function that may be safely wrapped with these
    typemaps (see below for more details).  In particular, none of the
    LAL Vectors, TimeSeries, or FrequencySeries that it takes as input may have
    their data pointers changed.  Any function which reallocs or resizes its
    inputs will cause such a problem.
(1) Add the header file in which it is defined to the section beneath the
    comment "Add LAL header files here".  It must be between the %{ and %}.
    The basic headers <lal/LALAtomicDatatypes.h> and <lal/LALDatatypes.h>
    have already been included.
(2) In the section beneath the comment "Add wrappings here" add the
    following:
    (a) %unignore(<name of function>);
    (b) %newobject directives for each variable that will have a
         NEWOUT or ARGOUT typemap applied to it.  Unlike the typemaps
	 themselves, these much each be specified by the name of the
	 variable; it is not possible (for example) to apply this
	 directive to everything declared as **COMPLEX8Vector.
    (c) %apply directives for any typemaps FOR ARRAY, TIMESERIES, OR
         FREQUENCYSERIES OBJECTS ONLY that the function accepts as input,
	 returns as output, or returns as an "output argument". Also ensure
	 that if the return value of a function is one of those types but it
	 should be ignored (generally because it returns one of its arguments,
	 but that argument is also modified in place) then you apply the
	 appropriate NONEOUT typemap to the function.
    (d) Finally, include the function declaration, as "extern"

In the resulting wrapped function---which will appear in the module
pycbc.lalwrap---all places where the LAL function accepts as input a LAL Vector,
TimeSeries, or FrequencySeries (of one of the four supported floating-point types),
the wrapped function will expect directly the appropriate object from pycbc.types.
If the function takes a double pointer to such an object because it will allocate
and return the object (this is common in functions that need to return more than
one such object) then the wrapped function will NOT accept an input for that
argument, but will instead return the corresponding object as part of a list or
tuple of return values.

So, for example, if we have the following declaration in <lal/MyStuff.h>:

  COMPLEX16FrequencySeries *XLALSomeFunc(REAL4TimeSeries **hplus,
                                         REAL4TimeSeries **hcross,
					 COMPLEX16FrequencySeries *input,
					 LIGOTimeGPS MyStart);

then to wrap it we would in all likelihood do the following:

(0) Check, as per above.  In this case making sure "input" is not resized
    by the function, either directly or through any functions it calls on its
    input.  We also determine whether the return value is the same as input,
    or a newly allocated frequency series.  For the sake of this example, assume
    the latter.
(1) We add the line:
       #include <lal/MyStuff.h>
    to the header section of the file.
(2) We add the following lines to the declaration section:
       %unignore(XLALSomeFunc);
       %newobject hplus;
       %newobject hcross;
       %newobject XLALSomeFunc;
       %apply REAL4TimeSeries **ARGOUT_REAL4TS {REAL4TimeSeries **};
       %apply COMPLEX16FrequencySeries *NEWOUT_COMPLEX16FS {COMPLEX16FrequencySeries *XLALSomeFunc};
       %apply COMPLEX16FrequencySeries *INPUT_COMPLEX16FS {COMPLEX16FrequencySeries *};
       extern COMPLEX16FrequencySeries *XLALSomeFunc(REAL4TimeSeries **hplus,
              REAL4TimeSeries **hcross,COMPLEX16FrequencySeries *input,
	      LIGOTimeGPS MyStart);

Note that we do not give a variable after the INPUT_COMPLEX16FS typemap; if there
were several such arguments this would cause them all to be treated the same
(likewise with the {REAL4TimeSeries **} declaration, which affects both hplus and
hcross).  Since we did not want that same typemap also applied to the output of the
function, we there specified the name.  The order does not matter; SWIG will apply
typemaps according to a "most specific match wins" rule.  But such typemap
declarations remain in effect until cleared, so all following pointers to
COMPLEX16FrequencySeries will be treated with the INPUT_COMPLEX16FS typemap, unless
the argument happens to be named "XLALSomeFunc" (which is unlikely). So if your
wrapping is behaving oddly, it could be due to a typemap left in effect from a
declaration above yours.  Note also that we do not specify any typemap for the
MyStart variable, since it is not a PyCBC basic type.

The resulting wrapped function will be pycbc.lalwrap.XLALSomeFunc, and it will
take exactly two arguments: an instance of pycbc.types.FrequencySeries (with
dtype of 'complex128') and an instance of lal.LIGOTimeGPS (i.e. the swiglal
wrapped struct).  It will return either a list or a tuple with, in order, the
return value of the function, then hplus, then hcross.  So you might call it in
Python as:

NewFS, hplus, hcross = pycbc.lalwrap.XLALSomeFunc(inputFS,MyEpoch)

assuming the inputs were of the appropriate type.

The available typemaps, for each of the four floating-point LAL array-like types,
are the following:

(1) INPUT_<TYPE>{V,TS,FS}    *<TYPE>{Vector,TimeSeries,FrequencySeries}
(2) NEWOUT_<TYPE>{V,TS,FS}   *<TYPE>{Vector,TimeSeries,FrequencySeries}
(3) NONEOUT_<TYPE>{V,TS,FS}  *<TYPE>{Vector,TimeSeries,FrequencySeries}
(4) ARGOUT_<TYPE>{V,TS,FS}  **<TYPE>{Vector,TimeSeries,FrequencySeries}

where <TYPE> is one of REAL4, REAL8, COMPLEX8, or COMPLEX16.

The INPUT typemap should be applied to any argument of that array-like type to a LAL function,
and will cause the corresponding Python wrapped function to expect (directly) an instance of the
appropriate type from pycbc.types.  The _data property of this input argument should be a
C-contiguous, one-dimensional Numpy array of the dtype appropriate to that vector. Note that
this object may in fact be treated by the XLAL function as "output", with its contents
modified in place. The elements of _data could have been modified, or also any of the object's
other properties that correspond to members of the corresponding LAL struct.

The NEWOUT typemaps are for when a function returns a newly-allocated vector
of that type. This function will be wrapped into one which returns a newly allocated
instance of pycbc.types.{Array,TimeSeries,FrequencySeries}, with its '_data' element set to the
appropriate numpy type, and avoiding any memory copies. Other properties will be set as appropriate
to that specific array-like type (e.g. epoch, delta_t, etc).

The NONEOUT typemap is for functions that return an object type, but that object is in
fact the same as one of the function's input arguments (which has been modified in
place).  In this case the wrapped function will always return "None" on successful
completion, and the effect of the function should be accessed through its effect on the
appropriate input object.

The ARGOUT typemap is for arguments that are double-pointers to an array-like type. These will
be wrapped into Python functions returning (possibly in a list, if there is more than one
return value) the new instances of the appropriate pycbc.type object created by those functions.
There will be no corresponding input argument to the wrapped function.

** IMPORTANT NOTE **

Not all XLAL functions will leave the location of their internal memory unchanged.
In particular any function that calls realloc() on an object's data (including
functions that "cut" or "resize" objects) can change the location of the data
pointer. Certainly anything than can change the size of 'self._data' is suspect, though
problematic functions are not necessarily limited to those. There is no way to accomodate
this behavior into Numpy's memory model, so you should ensure that you never SWIG-wrap
(with the typemaps of this file) any function that may change the location in memory of the
associated struct's data. The best you might expect is a difficult-to-reproduce segfault.

If you find you need to call such a function and it is not time-critical, one alternative
is to call the usual swiglal wrapped function and copy the output into a pycbc type for
further processing.

YOU HAVE BEEN WARNED!

*/


%module lalwrap

%{
// Add LAL header files here:
#include <lal/VectorOps.h>
#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>
%}

// DO NOT! change the following:
%include "pycbc_laltypes.i"

// Add wrappings here:

// An example of the usage of "NONEOUT": the vector returned by
// the function is the same as the function argument "out"
%unignore(XLALSSVectorMultiply);
%apply REAL4Vector *NONEOUT_REAL4V { REAL4Vector *XLALSSVectorMultiply };
%apply REAL4Vector *INPUT_REAL4V {REAL4Vector *};
extern REAL4Vector *XLALSSVectorMultiply(REAL4Vector *out, REAL4Vector *in1, REAL4Vector *in2);

// Functions to perform FFTs:
%unignore(XLALCOMPLEX8VectorFFT);
%apply COMPLEX8Vector *INPUT_COMPLEX8V {COMPLEX8Vector *};
extern int XLALCOMPLEX8VectorFFT(COMPLEX8Vector *output, COMPLEX8Vector *input, const COMPLEX8FFTPlan *plan );

%unignore(XLALCOMPLEX16VectorFFT);
%apply COMPLEX16Vector *INPUT_COMPLEX16V {COMPLEX16Vector *};
extern int XLALCOMPLEX16VectorFFT(COMPLEX16Vector *output, COMPLEX16Vector *input, const COMPLEX16FFTPlan *plan );

%unignore(XLALREAL4ForwardFFT);
%apply COMPLEX8Vector *INPUT_COMPLEX8V {COMPLEX8Vector *};
%apply REAL4Vector *INPUT_REAL4V {REAL4Vector *};
extern int XLALREAL4ForwardFFT(COMPLEX8Vector *output, REAL4Vector *input, REAL4FFTPlan *plan );

%unignore(XLALREAL4ReverseFFT);
%apply COMPLEX8Vector *INPUT_COMPLEX8V {COMPLEX8Vector *};
%apply REAL4Vector *INPUT_REAL4V {REAL4Vector *};
extern int XLALREAL4ReverseFFT(REAL4Vector *output, COMPLEX8Vector *input, REAL4FFTPlan *plan );

%unignore(XLALREAL4VectorFFT);
%apply REAL4Vector *INPUT_REAL4V {REAL4Vector *};
extern int XLALREAL4VectorFFT(REAL4Vector *output, REAL4Vector *input, REAL4FFTPlan *plan );

%unignore(XLALREAL8ForwardFFT);
%apply COMPLEX16Vector *INPUT_COMPLEX16V {COMPLEX16Vector *};
%apply REAL8Vector *INPUT_REAL8V {REAL8Vector *};
extern int XLALREAL8ForwardFFT(COMPLEX16Vector *output, REAL8Vector *input, REAL8FFTPlan *plan );

%unignore(XLALREAL8ReverseFFT);
%apply COMPLEX16Vector *INPUT_COMPLEX16V {COMPLEX16Vector *};
%apply REAL8Vector *INPUT_REAL8V {REAL8Vector *};
extern int XLALREAL8ReverseFFT(REAL8Vector *output, COMPLEX16Vector *input, REAL8FFTPlan *plan );

%unignore(XLALREAL8VectorFFT);
%apply REAL8Vector *INPUT_REAL8V {REAL8Vector *};
extern int XLALREAL8VectorFFT(REAL8Vector *output, REAL8Vector *input, REAL8FFTPlan *plan );

// Constructors for all four basic Vector types:
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
