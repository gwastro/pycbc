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

This is the pycbc_laltypes.i file.

It provides the typemaps used by PyCBC when time-critical code needs to
access LAL functions.  See the (PyCBC) file lalwrap.i for more instructions on
how to use the typemaps in this file to provide such a wrapping.  Even PyCBC
developers are not expected to routinely need to modify this file.

 */

%{
#define SWIG_FILE_WITH_INIT
#include <string.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/XLALError.h>
#include <numpy/arrayobject.h>
%}

// Ignore all of the names from the swiglal wrappings, so that they
// do not override what we will wrap.

%ignore "";

// Seems we must manually include <exception.i> before importing
// <lal/lalswig.i>.  Might be a bug in lalswig, or might just be
// unavoidable given that we %import rather than %include the
// lal wrappings.
%include <exception.i>
%import <lal/lalswig.i>

%begin %{
#define SWIG_TYPE_TABLE swiglal
%}

// We define global variables and the functions that will initialize
// them in the header portion of our interface file; the initialization
// functions are actually called in the %init section.

%header %{
// The next three global variables hold handles to the three basic
// PyCBC data types that can be used by C-API functions to create
// new instances of these objects when needed as return types by
// various typemaps below.  They can also be used by typemaps to check
// that input variables are in fact instances of the intended input
// types. These three are initialized by the function import_pycbc().
PyObject *CBC_Arr = NULL;
PyObject *CBC_TS = NULL;
PyObject *CBC_FS = NULL;
// The next global variables is simply an empty tuple that it is
// useful to have around for the various PyCBC constructors. The 'arg'
// parameter cannot be NULL, but for safety we always call with only
// keyword arguments.  Hence the empty tuple to be reused as the 'arg'
// parameter of constructors.  This is initialized in init_etuple().
PyObject *EmptyTuple = NULL;

// This code imports from pycbc.types the three basic
// datatypes of Array, TimeSeries, and FrequencySeries. They are
// imported as CBC_Arr, CBC_TS, and CBC_FS.
void import_pycbc(void){
  PyObject *CBC_ArrName, *CBC_TSName, *CBC_FSName;
  PyObject *CBC_TypesModule, *CBC_FromList, *CBC_Globals;

  CBC_ArrName = NULL;
  CBC_TSName = NULL;
  CBC_FSName = NULL;
  CBC_FromList = NULL;
  CBC_TypesModule = NULL;

  CBC_ArrName = PyString_FromString("Array");
  CBC_TSName = PyString_FromString("TimeSeries");
  CBC_FSName = PyString_FromString("FrequencySeries");
  if (!CBC_ArrName || !CBC_TSName || !CBC_FSName) goto fail;
  CBC_FromList = PyList_New(0);
  if (!CBC_FromList) goto fail;
  if (PyList_Append(CBC_FromList,CBC_ArrName)) goto fail;
  if (PyList_Append(CBC_FromList,CBC_TSName)) goto fail;
  if (PyList_Append(CBC_FromList,CBC_FSName)) goto fail;


  // Now we're ready to import the actual pycbc.types module

  // We only attempt to access the globals(), and just use NULL
  // where locals() would go in the module import, because
  // __import__ is documented as ignoring the locals() argument
  // and using globals only to determine package context. For
  // the same reason we explicitly create our global variables and
  // load them, rather than relying on typemaps to look them up in
  // either globals() or locals() (which would be slower, as well).

  CBC_Globals = PyEval_GetGlobals();
  Py_XINCREF(CBC_Globals); // Because we've just borrowed a reference
  CBC_TypesModule = PyImport_ImportModuleEx("pycbc.types",CBC_Globals,
					    NULL,CBC_FromList);
  Py_XDECREF(CBC_Globals);
  if (!CBC_TypesModule) goto fail;

  CBC_Arr = PyObject_GetAttr(CBC_TypesModule,CBC_ArrName);
  CBC_TS = PyObject_GetAttr(CBC_TypesModule,CBC_TSName);
  CBC_FS = PyObject_GetAttr(CBC_TypesModule,CBC_FSName);

  if (!CBC_Arr || !CBC_TS || !CBC_FS) goto fail;

  Py_DECREF(CBC_FromList);
  Py_DECREF(CBC_ArrName);
  Py_DECREF(CBC_TSName);
  Py_DECREF(CBC_FSName);
  Py_DECREF(CBC_TypesModule);

  return;

 fail:
  Py_XDECREF(CBC_Arr);
  Py_XDECREF(CBC_TS);
  Py_XDECREF(CBC_FS);
  Py_XDECREF(CBC_FromList);
  Py_XDECREF(CBC_ArrName);
  Py_XDECREF(CBC_TSName);
  Py_XDECREF(CBC_FSName);
  Py_XDECREF(CBC_TypesModule);
  PyErr_SetString(PyExc_ImportError,"Error importing 'pycbc.types'");
  PyErr_Print();
  return;
}

void init_etuple(void) {
  EmptyTuple = PyTuple_New(0);
  if (!EmptyTuple) {
    PyErr_SetString(PyExc_RuntimeError,"Error creating empty tuple for internal use");
    PyErr_Print();
    return;
  }
  return;
}

%}

%init {
  import_array();
  import_pycbc();
  init_etuple();
}

/*

Fragments

The next section contains several SWIG typemap
fragments, that are reused in several of the actual
typemaps.  They are essentially all either struct or
function definitions that are used or called when
sanity-checking the results of a typemap conversion
between C and Python.

*/

// Force inclusion of SWIG_From_int
// This is possibly a bug in the swiglal wrappings?

%fragment("SWIG_From_int");


%fragment("GenericVector","header") {
// The following struct is just so we can
// sanity check some output vectors without
// knowing exactly what datatype they are.
// This struct should not be SWIG-wrapped,
// nor used outside of the functions defined
// elsewhere in this file.
typedef struct {
  UINT4 length;
  void *data;
} GenericVector;

}

%fragment("GenericTS","header",fragment="GenericVector") {
// The following struct is just so we can
// sanity check some output time series without
// knowing exactly what datatype they are.
// This struct should not be SWIG-wrapped,
// nor used outside of the functions defined
// elsewhere in this file.
typedef struct {
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  LALUnit       sampleUnits;
  GenericVector *data;
} GenericTS;

}

%fragment("GenericFS","header",fragment="GenericVector") {
// The following struct is just so we can
// sanity check some output freq. series without
// knowing exactly what datatype they are.
// This struct should not be SWIG-wrapped,
// nor used outside of the functions defined
// elsewhere in this file.
typedef struct {
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         f0;
  REAL8         deltaF;
  LALUnit       sampleUnits;
  GenericVector *data;
} GenericFS;

}

// SWIG does not insert temporary object cleanup code into SWIG_Fail
// when called with "newout" code.  Thus all of our functions that
// marshall a newly allocated PyCBC type must destroy the temporary
// LAL struct created to hold the return; that cannot be done just
// in the typemap.  So we create common functions to do that here:

%fragment("ClearVect","header",fragment="GenericVector"){
  void ClearVect(GenericVector *self) {
    if (self) {
      free(self);
    }
    return;
  }
}

%fragment("ClearTS","header",fragment="GenericTS",fragment="ClearVect"){
  void ClearTS(GenericTS *self){
    if (self) {
      ClearVect(self->data);
      free(self);
    }
    return;
  }
}

%fragment("ClearFS","header",fragment="GenericFS",fragment="ClearVect"){
  void ClearFS(GenericFS *self){
    if (self) {
      ClearVect(self->data);
      free(self);
    }
    return;
  }
}

// The following fragment is used by any of the "ARGOUT" types to build
// a (possible) tuple of return values, since only in this case can
// there be more than one (that we must worry about).  It doesn't care
// at all what kind of Python objects its two inputs represent; the
// first is whatever output has been built up to now, and the second is
// the new value to add.

// This function directly calls a SWIG runtime function, and it seems
// impossible to tell whether this API is considered stable (as with most
// things in SWIG internals, it's completely undocumented).  Also, the
// SWIG function we call appears to perform no error checking on the Python
// C-API functions it calls.  However, its behavior depends on compile-time
// options of SWIG and this seems (at least for now) the easiest way to handle
// this.

%fragment("BuildReturnFromValue","header") {
  PyObject *BuildReturnFromValue(PyObject *CurrReturn, PyObject *value){
    // Since this function is often called directly before the second
    // argument can be checked in the caller, return failure if second
    // argument is NULL.  It's unclear what, if any, Python cleanup (i.e.
    // appropriate Py_DECREF's) is done in SWIG_Python_AppendOutput(). But
    // at least we do no worse than that function...
    if (!(value)) return NULL;
    return SWIG_Python_AppendOutput(CurrReturn,value);
  }
}

// Since all three of Array, TimeSeries, and FrequencySeries have a "_data" property
// that must be managed in the same way, we separate that handling into the following
// function, called by all "MarshallInput" functions.

%fragment("VectFromPyCBCType","header",fragment="GenericVector"){
  GenericVector *VectFromPyCBCType(PyObject *obj, const int numpy_type, const char *objname){
    GenericVector *returnptr;
    PyObject *tmpobj;

    tmpobj = PyObject_GetAttrString(obj,"_data");
    if (!tmpobj){
      PyErr_Format(PyExc_TypeError,
		   "Could not get _data property of argument '%s'", objname);
      return NULL;
    }
    if (!PyArray_Check(tmpobj)){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._data' must be a numpy array", objname);
      Py_DECREF(tmpobj);
      return NULL;
    }
    if (!(PyArray_ISCARRAY((PyArrayObject *) tmpobj)
	  || PyArray_ISCARRAY_RO((PyArrayObject *) tmpobj)) ){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._data' is not C-order contiguous",objname);
      Py_DECREF(tmpobj);
      return NULL;
    }
    if ( PyArray_NDIM((PyArrayObject *) tmpobj) != 1) {
      PyErr_Format(PyExc_ValueError,
		   "Argument '%s._data' is not one-dimensional",objname);
      Py_DECREF(tmpobj);
      return NULL;
    }
    if (PyArray_TYPE((PyArrayObject *) tmpobj) != numpy_type) {
      PyErr_Format(PyExc_ValueError,
		   "Argument '%s._data' has wrong dtype for corresponding LAL vector",objname);
      Py_DECREF(tmpobj);
      return NULL;
    }

    // Assemble our return object, which is a GenericVector.  If LAL should ever change
    // its definitions of Vectors so that different types change by more than underlying datatypes,
    // this would all have to be redone into a case statement based on the dtype.

    returnptr = (GenericVector *) calloc(1,sizeof(GenericVector));
    if (!returnptr) {
      PyErr_Format(PyExc_MemoryError,
		   "Could not allocate temporary Vector for argument '%s'",objname);
      Py_DECREF(tmpobj);
      return NULL;
    }
    returnptr->data = PyArray_DATA(tmpobj);
    returnptr->length = (UINT4) PyArray_DIM(tmpobj,0);

    // We should now release the reference count acquired through 'GetAttrString':

    Py_DECREF(tmpobj);

    return returnptr;
  }
}

// Similar to above, but for handling output.  This function returns (if
// successful) a PyDict with constructor entries for data, dtype, and copy;
// for pycbc.types.Array those are all that is needed.  Time and frequency
// series will further append to this dict before calling the constructor.

%fragment("ArrayDictFromVect","header",fragment="GenericVector",fragment="ClearVect"){
  PyObject *ArrayDictFromVect(GenericVector *vect, int numpy_type){
    PyObject *tmpobj, *constrdict;
    npy_intp dimensions[1];

    if (!(vect)) {
      PyErr_SetString(PyExc_ValueError,"Unexpected null vector returned from function");
      return NULL;
    }
    if ( (vect->length) &&  !(vect->data) ) {
      PyErr_SetString(PyExc_ValueError,"Null data pointer returned for non-zero length object");
      ClearVect(vect);
      return NULL;
    }

    constrdict = PyDict_New();
    if (!constrdict) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create dictionary for return value constructor");
      ClearVect(vect);
      return NULL;
    }

    dimensions[0] = (npy_intp) vect->length;
    tmpobj = PyArray_SimpleNewFromData(1,dimensions,numpy_type,(void *) vect->data);
    if (!tmpobj) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output data object");
      Py_DECREF(constrdict); // Dict still empty, so just delete
      ClearVect(vect);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"initial_array",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add data object to cosntructor dict");
      Py_DECREF(constrdict); // Dict still empty, so just delete
      Py_DECREF(tmpobj);
      ClearVect(vect);
      return NULL;
    }
    Py_DECREF(tmpobj);

    tmpobj = (PyObject *) PyArray_DescrFromType(numpy_type);
    if (!tmpobj){
      PyErr_SetString(PyExc_RuntimeError,"Could not create output dtype object");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      ClearVect(vect);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"dtype",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add dtype object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(tmpobj);
      ClearVect(vect);
      return NULL;
    }
    Py_DECREF(tmpobj);

    Py_INCREF(Py_False);
    tmpobj = Py_False;
    if (PyDict_SetItemString(constrdict,"copy",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add copy object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(tmpobj);
      ClearVect(vect);
      return NULL;
    }
    Py_DECREF(tmpobj);

    return constrdict;
  }
 }

// Common code for converting an instance of LIGOTimeGPS into a SWIG-wrapped
// Python object, with a dynamically allocated copy of the original object
// (which can therefore be garbage-collected by Python).

// Note: it is important that the PyObject *self be passed in, and have that
// name, for the SWIG runtime function SWIG_NewPointerObj() to work properly
// on the builtin LIGOTimeGPS wrapper that swiglal uses.  This is why that
// object pointer is passed not only to this function, but also from the actual
// typemaps into the marshalling functions that call the following function.

%fragment("EpochFromLAL","header") {
  PyObject *EpochFromLAL(PyObject *self,LIGOTimeGPS theepoch) {
    LIGOTimeGPS *epoch_ptr;
    PyObject *returnobj;

    epoch_ptr = calloc(1,sizeof(LIGOTimeGPS));
    if (!epoch_ptr) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output epoch object");
      return NULL;
    }
    epoch_ptr->gpsSeconds = theepoch.gpsSeconds;
    epoch_ptr->gpsNanoSeconds = theepoch.gpsNanoSeconds;
    returnobj = SWIG_NewPointerObj((void *) epoch_ptr,SWIG_TypeQuery("LIGOTimeGPS *"),SWIG_POINTER_OWN);
    if (!returnobj) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output epoch object");
      free(epoch_ptr);
    }
    return returnobj;
  }
}

%fragment("MarshallInputVector","header",fragment="VectFromPyCBCType") {
  GenericVector *MarshallInputVector(PyObject *obj, const int numpy_type, const char *objname) {
    GenericVector *returnptr;
    PyObject *tmpobj;

    tmpobj = PyObject_GetAttrString(obj,"_swighelper");

    // We explicitly access the '_swighelper' attribute of the argument, to force it onto
    // the CPU (if it was on the GPU and the current scheme is CPU) or to raise an
    // exception (if the current scheme is GPU).

    // We should have a '_swighelper' attribute, and it should point back to our argument, or
    // there's a problem.

    if (tmpobj != obj) {
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' has no '_swighelper' attribute---it is not an instance of pycbc.types.Array",
		   objname);
      Py_XDECREF(tmpobj);
      return NULL;
    }

    // If we get here, it means that the _swighelper property did behave as expected,
    // so to avoid an ever-increasing refcount, we must now decrement it:

    Py_DECREF(tmpobj);

    if (PyObject_IsInstance(obj,CBC_Arr) !=1){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' must be an instance of pycbc.types.Array or subclass", objname);
      return NULL;
    }

    returnptr = VectFromPyCBCType(obj,numpy_type,objname);

    return returnptr;
  }
}

%fragment("MarshallOutputVector","header",fragment="ArrayDictFromVect") {
  PyObject *MarshallOutputVector(GenericVector *vect, const int numpy_type) {
    PyObject *result, *constrdict;

    constrdict = ArrayDictFromVect(vect,numpy_type);
    if (!constrdict) {
      return NULL;
    }

    result = PyObject_Call(CBC_Arr,EmptyTuple,constrdict);
    if (!result) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create new instance of pycbc.types.Array");
      ClearVect(vect);
    }

    PyDict_Clear(constrdict);
    Py_DECREF(constrdict);

    return result;
  }
 }

%fragment("BuildArgoutVector","header",
	  fragment="BuildReturnFromValue",fragment="MarshallOutputVector") {};

%fragment("MarshallInputTS","header",fragment="GenericTS",fragment="VectFromPyCBCType") {
  GenericTS *MarshallInputTS(PyObject *obj, const int numpy_type, const char *objname) {
    GenericTS *returnptr;
    PyObject *tmpobj;
    LIGOTimeGPS *epoch_ptr;

    tmpobj = PyObject_GetAttrString(obj,"_swighelper");

    // We explicitly access the '_swighelper' attribute of the argument, to force it onto
    // the CPU (if it was on the GPU and the current scheme is CPU) or to raise an
    // exception (if the current scheme is GPU).

    // We should have a '_swighelper' attribute, and it should point back to our argument, or
    // there's a problem.

    if (tmpobj != obj) {
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' has wrong or missing '_swighelper' attribute---it is not an instance of pycbc.types.TimeSeries",
		   objname);
      Py_XDECREF(tmpobj);
      return NULL;
    }

    // If we get here, it means that the lal property did behave as expected, so to avoid
    // an ever-increasing refcount, we must now decrement it:

    Py_DECREF(tmpobj);

    if (PyObject_IsInstance(obj,CBC_TS) !=1){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' must be an instance of pycbc.types.TimeSeries or subclass", objname);
      return NULL;
    }

    // Start assembling our return object, which is a GenericTS.  If LAL should ever change
    // its definitions of TS so that different types change by more than underlying datatypes,
    // this would all have to be redone into a case statement based on the dtype.

    // Another important point: we use *calloc* not *malloc*.  This means that things we do not
    // read from our input (name, f0, and sampleUnits) are set to zero.  As this corresponds to
    // a blank string, 0.0, and lalDimensionlessUnit, respectively, this is OK.  Hence there is
    // no code below to set those members of our temporary struct that we pass to the wrapped
    // function.

    returnptr = (GenericTS *) calloc(1,sizeof(GenericTS));
    if (!returnptr) {
      PyErr_Format(PyExc_MemoryError,
		   "Could not allocate temporary TimeSeries for argument '%s'",objname);
      return NULL;
    }

    // First, marshall everything for the numpy array that is "self._data"
    returnptr->data = VectFromPyCBCType(obj,numpy_type,objname);
    if (!(returnptr->data)) {
      // VectFromPyCBCType has already set the error string, so just return
      return NULL;
    }

    // Next, marshall all of the other pieces of a TimeSeries that we do want to
    // get from our input.

    // First, epoch:

    tmpobj = PyObject_GetAttrString(obj,"_epoch");
    if (!tmpobj ||
	(SWIG_ConvertPtr(tmpobj,(void **) &epoch_ptr,SWIG_TypeQuery("LIGOTimeGPS *"),SWIG_POINTER_EXCEPTION) == -1)){
      Py_XDECREF(tmpobj);
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._epoch' does not exist or is not an instance of LIGOTimeGPS",objname);
      return NULL;
    }
    (returnptr->epoch).gpsSeconds = epoch_ptr->gpsSeconds;
    (returnptr->epoch).gpsNanoSeconds = epoch_ptr->gpsNanoSeconds;
    Py_DECREF(tmpobj);

    // Next, delta_t:

    tmpobj = PyObject_GetAttrString(obj,"_delta_t");
    if (!tmpobj || !PyFloat_Check(tmpobj)){
      Py_XDECREF(tmpobj);
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._delta_t' does not exist or is not a valid Python double",objname);
      return NULL;
    }
    // We use the macro form below---which doesn't check for errors--because
    // we've already checked that we have a float.
    returnptr->deltaT = (REAL8) PyFloat_AS_DOUBLE(tmpobj);
    Py_DECREF(tmpobj);

    // We're done, so return something non-NULL:

    return returnptr;
  }
}

%fragment("MarshallOutputTS","header",fragment="GenericTS",fragment="ClearTS",
	  fragment="ArrayDictFromVect",fragment="EpochFromLAL") {
  PyObject *MarshallOutputTS(PyObject *parentself, GenericTS *ts, const int numpy_type) {
    PyObject *result, *tmpobj, *constrdict;

    if (!(ts)) {
      PyErr_SetString(PyExc_ValueError,"Unexpected null time-series returned from function");
      return NULL;
    }

    constrdict = ArrayDictFromVect(ts->data,numpy_type);
    if (!constrdict) {
      // We have already cleared ts->data
      ts->data = NULL;
      ClearTS(ts);
      return NULL;
    }

    tmpobj = EpochFromLAL(parentself,ts->epoch);
    if (!tmpobj) { // Error message already set
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      ClearTS(ts);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"epoch",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add epoch object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(tmpobj);
      ClearTS(ts);
      return NULL;
    }
    Py_DECREF(tmpobj);

    tmpobj = PyFloat_FromDouble(ts->deltaT);
    if (!tmpobj) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output delta_t object");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      ClearTS(ts);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"delta_t",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add delta_t object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(tmpobj);
      ClearTS(ts);
      return NULL;
    }
    Py_DECREF(tmpobj);

    result = PyObject_Call(CBC_TS,EmptyTuple,constrdict);
    if (!result) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create new instance of pycbc.types.TimeSeries");
      ClearTS(ts);
    }
    // We don't need to do anything else for that last failure, as we'll be returning NULL
    // anyway and we have to do the same cleanup
    PyDict_Clear(constrdict);
    Py_DECREF(constrdict);

    return result;
  }
 }

%fragment("MarshallArgoutTS","header",fragment="GenericTS",fragment="ClearTS",
	  fragment="EpochFromLAL") {
  int MarshallArgoutTS(PyObject *parentself, PyObject *argument,
		       GenericTS *ts, const char *objname) {
    PyObject *tmpobj;

    if (!(ts)) {
      PyErr_Format(PyExc_ValueError,
		   "Unexpected null time-series for argument '%s' after return",objname);
      return 1;
    }
    if ( !(ts->data) ) {
      PyErr_Format(PyExc_ValueError,
		   "Time series argument '%s' had empty data vector after return",objname);
      ClearTS(ts);
      return 1;
    }
    if ( !(ts->data->data) &&  !(ts->data->length) ) {
      PyErr_Format(PyExc_ValueError,
		   "Time series argument '%s' had null data pointer returned for non-zero length",
		   objname);
      ClearTS(ts);
      return 1;
    }

    // The _data attribute should have automatically been modified in place, as it was
    // wrapped from a numpy array.  For all other elements of the returned TimeSeries,
    // we modify the existing attribute in place, since there's no easy way to know whether
    // or not the wrapped LAL function changed it.

    // The _epoch attribute:

    tmpobj = EpochFromLAL(parentself,ts->epoch);
    if (!tmpobj) {
      ClearTS(ts);
      return 1;
    }
    if (PyObject_SetAttrString(argument,"_epoch",tmpobj) == -1) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not modify '%s._epoch'",objname);
      Py_DECREF(tmpobj);
      ClearTS(ts);
      return 1;
    }
    Py_DECREF(tmpobj);

    // The _delta_t attribute:
    tmpobj = PyFloat_FromDouble((double) ts->deltaT);
    if (!tmpobj) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not create output delta_t object for argument '%s._delta_t'",objname);
      ClearTS(ts);
      return 1;
    }
    if (PyObject_SetAttrString(argument,"_delta_t",tmpobj) == -1) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not modify '%s._delta_t'",objname);
      Py_DECREF(tmpobj);
      ClearTS(ts);
      return 1;
    }
    Py_DECREF(tmpobj);

    return 0;
  }
}

%fragment("BuildArgoutTS","header",
	  fragment="BuildReturnFromValue",fragment="MarshallOutputTS") {};


%fragment("MarshallInputFS","header",fragment="GenericFS",fragment="VectFromPyCBCType") {
  GenericFS *MarshallInputFS(PyObject *obj, const int numpy_type, const char *objname) {
    GenericFS *returnptr;
    PyObject *tmpobj;
    LIGOTimeGPS *epoch_ptr;

    tmpobj = PyObject_GetAttrString(obj,"_swighelper");

    // We explicitly access the '_swighelper' attribute of the argument, to force it onto
    // the CPU (if it was on the GPU and the current scheme is CPU) or to raise an
    // exception (if the current scheme is GPU).

    // We should have a '_swighelper' attribute, and it should point back to our argument, or
    // there's a problem.

    if (tmpobj != obj) {
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' has wrong or missing '_swighelper' attribute---it is not an instance of pycbc.types.FrequencySeries",
		   objname);
      Py_XDECREF(tmpobj);
      return NULL;
    }
    Py_DECREF(tmpobj);

    if (PyObject_IsInstance(obj,CBC_FS) !=1){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' must be an instance of pycbc.types.FrequencySeries or subclass", objname);
      return NULL;
    }

    // Start assembling our return object, which is a GenericFS.  If LAL should ever change
    // its definitions of FS so that different types change by more than underlying datatypes,
    // this would all have to be redone into a case statement based on the dtype.

    // Another important point: we use *calloc* not *malloc*.  This means that things we do not
    // read from our input (name, f0, and sampleUnits) are set to zero.  As this corresponds to
    // a blank string, 0.0, and lalDimensionlessUnit, respectively, this is OK.  Hence there is
    // no code below to set those members of our temporary struct that we pass to the wrapped
    // function.

    returnptr = (GenericFS *) calloc(1,sizeof(GenericFS));
    if (!returnptr) {
      PyErr_Format(PyExc_MemoryError,
		   "Could not allocate temporary FrequencySeries for argument '%s'",objname);
      return NULL;
    }

    // First, marshall everything for the numpy array that is "self._data"
    returnptr->data = VectFromPyCBCType(obj,numpy_type,objname);
    if (!(returnptr->data)) {
      // VectFromPyCBCType has already set the error string, so just return
      return NULL;
    }

    // Next, marshall all of the other pieces of a FrequencySeries that we do want to
    // get from our input.

    // First, epoch:

    tmpobj = PyObject_GetAttrString(obj,"_epoch");
    if (!tmpobj ||
	(SWIG_ConvertPtr(tmpobj,(void **) &epoch_ptr,SWIG_TypeQuery("LIGOTimeGPS *"),SWIG_POINTER_EXCEPTION) == -1)){
      Py_XDECREF(tmpobj);
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._epoch' does not exist or is not an instance of LIGOTimeGPS",objname);
      return NULL;
    }
    (returnptr->epoch).gpsSeconds = epoch_ptr->gpsSeconds;
    (returnptr->epoch).gpsNanoSeconds = epoch_ptr->gpsNanoSeconds;
    Py_DECREF(tmpobj);


    /****************************************************************************************
       The code below is commented out and not used, but retained here in case it is decided
       to shadow the "f0" member in LAL FrequencySeries structs in the PyCBC type.

    // Next, f0:

    tmpobj = PyObject_GetAttrString(obj,"_f0");
    if (!tmpobj || !PyFloat_Check(tmpobj)){
      Py_XDECREF(tmpobj);
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._f0' does not exist or is not a valid Python double",objname);
      return NULL;
    }
    // We use the macro form below---which doesn't check for errors--because
    // we've already checked that we have a float.
    returnptr->f0 = (REAL8) PyFloat_AS_DOUBLE(tmpobj);
    Py_DECREF(tmpobj);
    ****************************************************************************************/

    // Finally, delta_f:

    tmpobj = PyObject_GetAttrString(obj,"_delta_f");
    if (!tmpobj || !PyFloat_Check(tmpobj)){
      Py_XDECREF(tmpobj);
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._delta_f' does not exist or is not a valid Python double",objname);
      return NULL;
    }
    // We use the macro form below---which doesn't check for errors--because
    // we've already checked that we have a float.
    returnptr->deltaF = (REAL8) PyFloat_AS_DOUBLE(tmpobj);
    Py_DECREF(tmpobj);

    // We're done, so return something non-NULL:

    return returnptr;
  }
}

%fragment("MarshallOutputFS","header",fragment="GenericFS",fragment="ClearFS",
	  fragment="ArrayDictFromVect",fragment="EpochFromLAL") {
  PyObject *MarshallOutputFS(PyObject *parentself, GenericFS *fs, const int numpy_type) {
    PyObject *result, *tmpobj, *constrdict;

    if (!(fs)) {
      PyErr_SetString(PyExc_ValueError,"Unexpected null frequency-series returned from function");
      return NULL;
    }

    constrdict = ArrayDictFromVect(fs->data,numpy_type);
    if (!constrdict) {
      fs->data = NULL;
      ClearFS(fs);
      return NULL;
    }

    tmpobj = EpochFromLAL(parentself,fs->epoch);
    if (!tmpobj) { // Error message already set
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      ClearFS(fs);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"epoch",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add epoch object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(tmpobj);
      ClearFS(fs);
      return NULL;
    }
    Py_DECREF(tmpobj);

    tmpobj = PyFloat_FromDouble(fs->deltaF);
    if (!tmpobj) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output delta_f object");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      ClearFS(fs);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"delta_f",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add delta_f object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(tmpobj);
      ClearFS(fs);
      return NULL;
    }
    Py_DECREF(tmpobj);

    /************************************************************************************

    The code below is commented out and not used, but retained here in case it is decided
    to shadow the "f0" member in LAL FrequencySeries structs in the PyCBC type.

    tmpobj = PyFloat_FromDouble(fs->f0);
    if (!tmpobj) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output f0 object");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      ClearFS(fs);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"f0",tmpobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add f0 object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(tmpobj);
      ClearFS(fs);
      return NULL;
    }
    Py_DECREF(tmpobj);

    *************************************************************************************/

    result = PyObject_Call(CBC_FS,EmptyTuple,constrdict);
    if (!result) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create new instance of pycbc.types.FrequencySeries");
      ClearFS(fs);
    }
    // We don't need to do anything else for that last failure, as we'll be returning NULL
    // anyway and we have to do the same cleanup
    PyDict_Clear(constrdict);
    Py_DECREF(constrdict);

    return result;
  }
 }

%fragment("MarshallArgoutFS","header",fragment="GenericFS",fragment="ClearTS",
	  fragment="EpochFromLAL") {
  int MarshallArgoutFS(PyObject *parentself, PyObject *argument,
		       GenericFS *fs, const char *objname) {
    PyObject *tmpobj;

    if (!(fs)) {
      PyErr_Format(PyExc_ValueError,
		   "Unexpected null frequency-series for argument '%s' after return",objname);
      return 1;
    }
    if ( !(fs->data) ) {
      PyErr_Format(PyExc_ValueError,
		   "Frequency series argument '%s' had empty data vector after return",objname);
      ClearFS(fs);
      return 1;
    }
    if ( !(fs->data->data) &&  !(fs->data->length) ) {
      PyErr_Format(PyExc_ValueError,
		   "Frequency series argument '%s' had null data pointer returned for non-zero length",
		   objname);
      ClearFS(fs);
      return 1;
    }

    // The _data attribute should have automatically been modified in place, as it was
    // wrapped from a numpy array.  For all other elements of the returned FrequencySeries,
    // we modify the existing attribute in place, since there's no easy way to know whether
    // or not the wrapped LAL function changed it.

    // The _epoch attribute:
    tmpobj = EpochFromLAL(parentself,fs->epoch);
    if (!tmpobj) {
      ClearFS(fs);
      return 1;
    }
    if (PyObject_SetAttrString(argument,"_epoch",tmpobj) == -1) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not modify '%s._epoch'",objname);
      Py_DECREF(tmpobj);
      ClearFS(fs);
      return 1;
    }
    Py_DECREF(tmpobj);


    // The _delta_f attribute:
    tmpobj = PyFloat_FromDouble((double) fs->deltaF);
    if (!tmpobj) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not create output delta_f object for argument '%s._delta_f'",objname);
      ClearFS(fs);
      return 1;
    }
    if (PyObject_SetAttrString(argument,"_delta_f",tmpobj) == -1) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not modify '%s._delta_f'",objname);
      Py_DECREF(tmpobj);
      ClearFS(fs);
      return 1;
    }
    Py_DECREF(tmpobj);

    /***********************************************************************************

    The code below is commented out and not used, but retained here in case it is decided
    to shadow the "f0" member in LAL FrequencySeries structs in the PyCBC type.

    // The _f0 attribute:
    tmpobj = PyFloat_FromDouble((double) fs->f0);
    if (!tmpobj) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not create output f0 object for argument '%s._f0'",objname);
      ClearFS(fs);
      return 1;
    }
    if (PyObject_SetAttrString(argument,"_f0",tmpobj) == -1) {
      PyErr_Format(PyExc_RuntimeError,
		   "Could not modify '%s._f0'",objname);
      Py_DECREF(tmpobj);
      ClearFS(fs);
      return 1;
    }
    Py_DECREF(tmpobj);

    *************************************************************************************/

    return 0;
  }
}

%fragment("BuildArgoutFS","header",
	  fragment="BuildReturnFromValue",fragment="MarshallOutputFS") {};

// Force inclusion of all of our marshalling functions, to help in debugging

%fragment("MarshallInputVector");
%fragment("MarshallOutputVector");
%fragment("MarshallInputTS");
%fragment("MarshallOutputTS");
%fragment("MarshallArgoutTS");
%fragment("MarshallInputFS");
%fragment("MarshallOutputFS");
%fragment("MarshallArgoutFS");

/*

pycbc.types LAL typemaps


Below we have typemaps for each of the four floating-point LAL array-like types.
For each type there are four sets of typemaps:

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
other properties that correspond to members of the corresponding XLAL struct.

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

All of the typemaps except the NONEOUT typemaps have additional features which
will deallocate the temporary struct needed to pass the necessary memory between
Python and LAL.  These typemaps ("freearg" and "newfree") should not be modified
or disabled or memory leaks will result.

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


// Typemaps for REAL4Vectors:

%typemap(in, fragment="MarshallInputVector") REAL4Vector *INPUT_REAL4V {
  $1 =(REAL4Vector *) MarshallInputVector($input,NPY_FLOAT32,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(freearg) REAL4Vector *INPUT_REAL4V {
  ClearVect((GenericVector *) $1);
}

%typemap(out, fragment="MarshallOutputVector") REAL4Vector *NEWOUT_REAL4V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_FLOAT32);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL4Vector *NEWOUT_REAL4V{
  ClearVect((GenericVector *) $1);
}

%typemap(out) REAL4Vector *NONEOUT_REAL4V{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) REAL4Vector **ARGOUT_REAL4V (REAL4Vector *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutVector") REAL4Vector **ARGOUT_REAL4V {
  $result = BuildReturnFromValue($result,
				 MarshallOutputVector((GenericVector *) *($1),NPY_FLOAT32));
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL4Vector **ARGOUT_REAL4V {
  ClearVect((GenericVector *) *($1));
}


// Typemaps for REAL8Vectors:

%typemap(in, fragment="MarshallInputVector") REAL8Vector *INPUT_REAL8V {
  $1 =(REAL8Vector *) MarshallInputVector($input,NPY_FLOAT64,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(freearg) REAL8Vector *INPUT_REAL8V {
  ClearVect((GenericVector *) $1);
}

%typemap(out, fragment="MarshallOutputVector") REAL8Vector *NEWOUT_REAL8V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_FLOAT64);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL8Vector *NEWOUT_REAL8V{
  ClearVect((GenericVector *) $1);
}

%typemap(out) REAL8Vector *NONEOUT_REAL8V{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) REAL8Vector **ARGOUT_REAL8V (REAL8Vector *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutVector") REAL8Vector **ARGOUT_REAL8V {
  $result = BuildReturnFromValue($result,
				 MarshallOutputVector((GenericVector *) *($1),NPY_FLOAT64));
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL8Vector **ARGOUT_REAL8V {
  ClearVect((GenericVector *) *($1));
}

// Typemaps for COMPLEX8Vectors:

%typemap(in, fragment="MarshallInputVector") COMPLEX8Vector *INPUT_COMPLEX8V {
  $1 =(COMPLEX8Vector *) MarshallInputVector($input,NPY_COMPLEX64,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(freearg) COMPLEX8Vector *INPUT_COMPLEX8V {
  ClearVect((GenericVector *) $1);
}

%typemap(out, fragment="MarshallOutputVector") COMPLEX8Vector *NEWOUT_COMPLEX8V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_COMPLEX64);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) COMPLEX8Vector *NEWOUT_COMPLEX8V{
  ClearVect((GenericVector *) $1);
}

%typemap(out) COMPLEX8Vector *NONEOUT_COMPLEX8V{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) COMPLEX8Vector **ARGOUT_COMPLEX8V (COMPLEX8Vector *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutVector") COMPLEX8Vector **ARGOUT_COMPLEX8V {
  $result = BuildReturnFromValue($result,
				 MarshallOutputVector((GenericVector *) *($1),NPY_COMPLEX64));
  if (!($result)) SWIG_fail;
}

%typemap(newfree) COMPLEX8Vector **ARGOUT_COMPLEX8V {
  ClearVect((GenericVector *) *($1));
}

// Typemaps for COMPLEX16Vectors:

%typemap(in, fragment="MarshallInputVector") COMPLEX16Vector *INPUT_COMPLEX16V {
  $1 =(COMPLEX16Vector *) MarshallInputVector($input,NPY_COMPLEX128,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(freearg) COMPLEX16Vector *INPUT_COMPLEX16V {
  ClearVect((GenericVector *) $1);
}

%typemap(out, fragment="MarshallOutputVector") COMPLEX16Vector *NEWOUT_COMPLEX16V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_COMPLEX128);
  if (!($result)) SWIG_fail;
 }

%typemap(newfree) COMPLEX16Vector *NEWOUT_COMPLEX16V{
  ClearVect((GenericVector *) $1);
}

%typemap(out) COMPLEX16Vector *NONEOUT_COMPLEX16V{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) COMPLEX16Vector **ARGOUT_COMPLEX16V (COMPLEX16Vector *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutVector") COMPLEX16Vector **ARGOUT_COMPLEX16V {
  $result = BuildReturnFromValue($result,
				 MarshallOutputVector((GenericVector *) *($1),NPY_COMPLEX128));
  if (!($result)) SWIG_fail;
}

%typemap(newfree) COMPLEX16Vector **ARGOUT_COMPLEX16V {
  ClearVect((GenericVector *) *($1));
}

// Typemaps for time and frequency series.  Note that those calling the functions
// MarshallOutput{T,F}S and MarshallArgout{T,F}S have as their first argument "self".
// This variable is defined in the wrapped function, and hence will exist even though
// it is not a $special SWIG variable.


// Typemaps for REAL4 Time Series

%typemap(in, fragment="MarshallInputTS") REAL4TimeSeries *INPUT_REAL4TS {
  $1 =(REAL4TimeSeries *) MarshallInputTS($input,NPY_FLOAT32,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutTS") REAL4TimeSeries *INPUT_REAL4TS {
  if (MarshallArgoutTS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) REAL4TimeSeries *INPUT_REAL4TS {
  ClearTS($1);
}

%typemap(out, fragment="MarshallOutputTS") REAL4TimeSeries *NEWOUT_REAL4TS{
  $result = MarshallOutputTS(self,(GenericTS *) $1,NPY_FLOAT32);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL4TimeSeries *NEWOUT_REAL4TS{
  ClearTS($1);
}

%typemap(out) REAL4TimeSeries *NONEOUT_REAL4TS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) REAL4TimeSeries **ARGOUT_REAL4TS (REAL4TimeSeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutTS") REAL4TimeSeries **ARGOUT_REAL4TS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputTS(self,(GenericTS *) *($1),NPY_FLOAT32));
  if (!($result)) SWIG_fail;
}

// Typemaps for REAL8 Time Series

%typemap(in, fragment="MarshallInputTS") REAL8TimeSeries *INPUT_REAL8TS {
  $1 =(REAL8TimeSeries *) MarshallInputTS($input,NPY_FLOAT64,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutTS") REAL8TimeSeries *INPUT_REAL8TS {
  if (MarshallArgoutTS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) REAL8TimeSeries *INPUT_REAL8TS {
  ClearTS($1);
}

%typemap(out, fragment="MarshallOutputTS") REAL8TimeSeries *NEWOUT_REAL8TS{
  $result = MarshallOutputTS(self,(GenericTS *) $1,NPY_FLOAT64);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL8TimeSeries *NEWOUT_REAL8TS{
  ClearTS($1);
}

%typemap(out) REAL8TimeSeries *NONEOUT_REAL8TS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) REAL8TimeSeries **ARGOUT_REAL8TS (REAL8TimeSeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutTS") REAL8TimeSeries **ARGOUT_REAL8TS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputTS(self,(GenericTS *) *($1),NPY_FLOAT64));
  if (!($result)) SWIG_fail;
}


// Typemaps for COMPLEX8 Time Series

%typemap(in, fragment="MarshallInputTS") COMPLEX8TimeSeries *INPUT_COMPLEX8TS {
  $1 =(COMPLEX8TimeSeries *) MarshallInputTS($input,NPY_COMPLEX64,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutTS") COMPLEX8TimeSeries *INPUT_COMPLEX8TS {
  if (MarshallArgoutTS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) COMPLEX8TimeSeries *INPUT_COMPLEX8TS {
  ClearTS($1);
}

%typemap(out, fragment="MarshallOutputTS") COMPLEX8TimeSeries *NEWOUT_COMPLEX8TS{
  $result = MarshallOutputTS(self,(GenericTS *) $1,NPY_COMPLEX64);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) COMPLEX8TimeSeries *NEWOUT_COMPLEX8TS{
  ClearTS($1);
}

%typemap(out) COMPLEX8TimeSeries *NONEOUT_COMPLEX8TS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) COMPLEX8TimeSeries **ARGOUT_COMPLEX8TS (COMPLEX8TimeSeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutTS") COMPLEX8TimeSeries **ARGOUT_COMPLEX8TS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputTS(self,(GenericTS *) *($1),NPY_COMPLEX64));
  if (!($result)) SWIG_fail;
}

// Typemaps for COMPLEX16 Time Series

%typemap(in, fragment="MarshallInputTS") COMPLEX16TimeSeries *INPUT_COMPLEX16TS {
  $1 =(COMPLEX16TimeSeries *) MarshallInputTS($input,NPY_COMPLEX128,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutTS") COMPLEX16TimeSeries *INPUT_COMPLEX16TS {
  if (MarshallArgoutTS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) COMPLEX16TimeSeries *INPUT_COMPLEX16TS {
  ClearTS($1);
}

%typemap(out, fragment="MarshallOutputTS") COMPLEX16TimeSeries *NEWOUT_COMPLEX16TS{
  $result = MarshallOutputTS(self,(GenericTS *) $1,NPY_COMPLEX128);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) COMPLEX16TimeSeries *NEWOUT_COMPLEX16TS{
  ClearTS($1);
}

%typemap(out) COMPLEX16TimeSeries *NONEOUT_COMPLEX16TS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) COMPLEX16TimeSeries **ARGOUT_COMPLEX16TS (COMPLEX16TimeSeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutTS") COMPLEX16TimeSeries **ARGOUT_COMPLEX16TS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputTS(self,(GenericTS *) *($1),NPY_COMPLEX128));
  if (!($result)) SWIG_fail;
}

// Typemaps for REAL4 Frequency Series

%typemap(in, fragment="MarshallInputFS") REAL4FrequencySeries *INPUT_REAL4FS {
  $1 =(REAL4FrequencySeries *) MarshallInputFS($input,NPY_FLOAT32,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutFS") REAL4FrequencySeries *INPUT_REAL4FS {
  if (MarshallArgoutFS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) REAL4FrequencySeries *INPUT_REAL4FS {
  ClearFS($1);
}

%typemap(out, fragment="MarshallOutputFS") REAL4FrequencySeries *NEWOUT_REAL4FS{
  $result = MarshallOutputFS(self,(GenericFS *) $1,NPY_FLOAT32);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL4FrequencySeries *NEWOUT_REAL4FS{
  ClearFS($1);
}

%typemap(out) REAL4FrequencySeries *NONEOUT_REAL4FS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) REAL4FrequencySeries **ARGOUT_REAL4FS (REAL4FrequencySeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutFS") REAL4FrequencySeries **ARGOUT_REAL4FS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputFS(self,(GenericFS *) *($1),NPY_FLOAT32));
  if (!($result)) SWIG_fail;
}

// Typemaps for REAL8 Frequency Series

%typemap(in, fragment="MarshallInputFS") REAL8FrequencySeries *INPUT_REAL8FS {
  $1 =(REAL8FrequencySeries *) MarshallInputFS($input,NPY_FLOAT64,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutFS") REAL8FrequencySeries *INPUT_REAL8FS {
  if (MarshallArgoutFS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) REAL8FrequencySeries *INPUT_REAL8FS {
  ClearFS($1);
}

%typemap(out, fragment="MarshallOutputFS") REAL8FrequencySeries *NEWOUT_REAL8FS{
  $result = MarshallOutputFS(self,(GenericFS *) $1,NPY_FLOAT64);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL8FrequencySeries *NEWOUT_REAL8FS{
  ClearFS($1);
}

%typemap(out) REAL8FrequencySeries *NONEOUT_REAL8FS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) REAL8FrequencySeries **ARGOUT_REAL8FS (REAL8FrequencySeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutFS") REAL8FrequencySeries **ARGOUT_REAL8FS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputFS(self,(GenericFS *) *($1),NPY_FLOAT64));
  if (!($result)) SWIG_fail;
}


// Typemaps for COMPLEX8 Frequency Series

%typemap(in, fragment="MarshallInputFS") COMPLEX8FrequencySeries *INPUT_COMPLEX8FS {
  $1 =(COMPLEX8FrequencySeries *) MarshallInputFS($input,NPY_COMPLEX64,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutFS") COMPLEX8FrequencySeries *INPUT_COMPLEX8FS {
  if (MarshallArgoutFS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) COMPLEX8FrequencySeries *INPUT_COMPLEX8FS {
  ClearFS($1);
}

%typemap(out, fragment="MarshallOutputFS") COMPLEX8FrequencySeries *NEWOUT_COMPLEX8FS{
  $result = MarshallOutputFS(self,(GenericFS *) $1,NPY_COMPLEX64);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) COMPLEX8FrequencySeries *NEWOUT_COMPLEX8FS{
  ClearFS($1);
}

%typemap(out) COMPLEX8FrequencySeries *NONEOUT_COMPLEX8FS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) COMPLEX8FrequencySeries **ARGOUT_COMPLEX8FS (COMPLEX8FrequencySeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutFS") COMPLEX8FrequencySeries **ARGOUT_COMPLEX8FS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputFS(self,(GenericFS *) *($1),NPY_COMPLEX64));
  if (!($result)) SWIG_fail;
}

// Typemaps for COMPLEX16 Frequency Series

%typemap(in, fragment="MarshallInputFS") COMPLEX16FrequencySeries *INPUT_COMPLEX16FS {
  $1 =(COMPLEX16FrequencySeries *) MarshallInputFS($input,NPY_COMPLEX128,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(argout, fragment="MarshallArgoutFS") COMPLEX16FrequencySeries *INPUT_COMPLEX16FS {
  if (MarshallArgoutFS(self,$input,$1,"$1_name") ) SWIG_fail;
}

%typemap(freearg) COMPLEX16FrequencySeries *INPUT_COMPLEX16FS {
  ClearFS($1);
}

%typemap(out, fragment="MarshallOutputFS") COMPLEX16FrequencySeries *NEWOUT_COMPLEX16FS{
  $result = MarshallOutputFS(self,(GenericFS *) $1,NPY_COMPLEX128);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) COMPLEX16FrequencySeries *NEWOUT_COMPLEX16FS{
  ClearFS($1);
}

%typemap(out) COMPLEX16FrequencySeries *NONEOUT_COMPLEX16FS{
  Py_INCREF(Py_None);
  $result = Py_None;
}

%typemap(in,numinputs=0) COMPLEX16FrequencySeries **ARGOUT_COMPLEX16FS (COMPLEX16FrequencySeries *temp) {
  temp = NULL;
  $1 = &temp;
}

%typemap(argout, fragment="BuildArgoutFS") COMPLEX16FrequencySeries **ARGOUT_COMPLEX16FS {
  $result = BuildReturnFromValue($result,
				 MarshallOutputFS(self,(GenericFS *) *($1),NPY_COMPLEX128));
  if (!($result)) SWIG_fail;
}

// The following statment is more intuitive for developers to use than the direct
// use of the %rename("%s") syntax of SWIG.

%define %unignore(NAME)
%rename("%s") NAME;
%enddef
