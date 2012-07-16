%module testlal

%{
#define SWIG_FILE_WITH_INIT
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <numpy/arrayobject.h>
#include <lal/ComplexFFT.h>
#include <lal/XLALError.h>
%}

%ignore XLALCOMPLEX8VectorFFT;
%ignore XLALSSVectorMultiply;
%ignore XLALCreateREAL4Vector;

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

%fragment("MarshallInputVector","header") {
  void *MarshallInputVector(PyObject *obj, const int numpy_type, const char *objname) {
    void *returnptr;
    PyObject *dataobj, *tmpobj;

    tmpobj = PyObject_GetAttrString(obj,"lal");

    // We explicitly access the 'lal' attribute of the argument, to force it onto
    // the CPU (if it was on the GPU and the current scheme is CPU) or to raise an
    // exception (if the current scheme is GPU).

    // We should have a 'lal' attribute, and it should point back to our argument, or
    // there's a problem.

    if (tmpobj != obj) {
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' has no 'lal' attribute---it is not an instance of pycbc.types.Array",
		   objname);
      return NULL;
    }

    // If we get here, it means that the lal property did behave as expected, so to avoid
    // an ever-increasing refcount, we must now decrement it:

    Py_DECREF(tmpobj);

    if (PyObject_IsInstance(obj,CBC_Arr) !=1){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' must be an instance of pycbc.types.Array or subclass", objname);
      return NULL;
    }

    dataobj = PyObject_GetAttrString(obj,"_data");
    if (!dataobj){
      PyErr_Format(PyExc_TypeError,
		   "Could not get _data property of Argument '%s'", objname);
      return NULL;
    }
    if (!PyArray_Check(dataobj)){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._data' must be a numpy array", objname);
      return NULL;
    }
    if (!(PyArray_ISCARRAY((PyArrayObject *) dataobj)
	  || PyArray_ISCARRAY_RO((PyArrayObject *) dataobj)) ){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._data' is not C-order contiguous",objname);
      return NULL;
    }
    if ( PyArray_NDIM((PyArrayObject *) dataobj) != 1) {
      PyErr_Format(PyExc_ValueError,
		   "Argument '%s._data' is not one-dimensional",objname);
      return NULL;
    }
    if (PyArray_TYPE((PyArrayObject *) dataobj) != numpy_type) {
      PyErr_Format(PyExc_ValueError,
		   "Argument '%s._data' has wrong dtype for corresponding LAL vector",objname);
      return NULL;
    }

    switch(numpy_type){
    case NPY_FLOAT32:
      returnptr = calloc(1,sizeof(REAL4Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary REAL4Vector for argument '%s'",objname);
	return NULL;
      }
      ((REAL4Vector *) returnptr)->data = (REAL4 *) PyArray_DATA(dataobj);
      ((REAL4Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    case NPY_FLOAT64:
      returnptr = calloc(1,sizeof(REAL8Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary REAL8Vector for argument '%s'",objname);
	return NULL;
      }
      ((REAL8Vector *) returnptr)->data = (REAL8 *) PyArray_DATA(dataobj);
      ((REAL8Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    case NPY_COMPLEX64:
      returnptr = calloc(1,sizeof(COMPLEX8Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary COMPLEX8Vector for argument '%s'",objname);
	return NULL;
      }
      ((COMPLEX8Vector *) returnptr)->data = (COMPLEX8 *) PyArray_DATA(dataobj);
      ((COMPLEX8Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    case NPY_COMPLEX128:
      returnptr = calloc(1,sizeof(COMPLEX16Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary COMPLEX16Vector for argument '%s'",objname);
	return NULL;
      }
      ((COMPLEX16Vector *) returnptr)->data = (COMPLEX16 *) PyArray_DATA(dataobj);
      ((COMPLEX16Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    default:
      // If we get here, *we've* made a mistake in the wrapping---can't be user's fault.
      PyErr_Format(PyExc_TypeError,
		   "Argument array '%s._data' associated with wrong numpy type---this is a PyCBC bug",objname);
      return NULL;
    }
  return returnptr;
  }
}

%fragment("MarshallOutputVector","header",fragment="GenericVector") {
  PyObject *MarshallOutputVector(GenericVector *vect, const int numpy_type) {
    PyObject *result, *dataobj, *dtypeobj, *copybool, *constrdict;


    if (!(vect)) {
      PyErr_SetString(PyExc_ValueError,"Unexpected null vector returned from function");
      return NULL;
    }
    if ( (vect->length) &&  !(vect->data) ) {
      PyErr_SetString(PyExc_ValueError,"Null data pointer returned for non-zero length array");
      return NULL;
    }

    constrdict = PyDict_New();
    if (!constrdict) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create dictionary for return value constructor");
      return NULL;
    }

    npy_intp dimensions[1];
    dimensions[0] = (npy_intp) vect->length;
    dataobj = PyArray_SimpleNewFromData(1,dimensions,numpy_type,(void *) vect->data);
    if (!dataobj) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output data object");
      Py_DECREF(constrdict); // Dict still empty, so just delete
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"initial_array",dataobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add data object to cosntructor dict");
      Py_DECREF(constrdict); // Dict still empty, so just delete
      Py_DECREF(dataobj);
      return NULL;
    }

    dtypeobj = (PyObject *) PyArray_DescrNewFromType(numpy_type);
    if (!dtypeobj){
      PyErr_SetString(PyExc_RuntimeError,"Could not create output dtype object");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(dataobj);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"dtype",dtypeobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add dtype object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(dataobj);
      Py_DECREF(dtypeobj);
      return NULL;
    }

    Py_INCREF(Py_False);
    copybool = Py_False;
    if (PyDict_SetItemString(constrdict,"copy",copybool)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add copy object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(dataobj);
      Py_DECREF(dtypeobj);
      Py_DECREF(copybool);
      return NULL;
    }


    result = PyObject_Call(CBC_Arr,EmptyTuple,constrdict);
    if (!result) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create new instance of pycbc.types.Array");
    }
    // We don't need to do anything else for that last failure, as we'll be returning NULL
    // anyway and we have to do the same cleanup
    PyDict_Clear(constrdict);
    Py_DECREF(constrdict);
    Py_DECREF(dataobj);
    Py_DECREF(dtypeobj);
    Py_DECREF(copybool);

    return result;
  }
 }

%fragment("MarshallInputTS","header") {
  void *MarshallInputTS(PyObject *obj, const int numpy_type, const char *objname) {
    void *returnptr;
    PyObject *dataobj, *tmpobj;

    tmpobj = PyObject_GetAttrString(obj,"lal");

    // We explicitly access the 'lal' attribute of the argument, to force it onto
    // the CPU (if it was on the GPU and the current scheme is CPU) or to raise an
    // exception (if the current scheme is GPU).

    // We should have a 'lal' attribute, and it should point back to our argument, or
    // there's a problem.

    if (tmpobj != obj) {
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s' has no 'lal' attribute---it is not an instance of pycbc.types.TimeSeries",
		   objname);
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

    dataobj = PyObject_GetAttrString(obj,"_data");
    if (!dataobj){
      PyErr_Format(PyExc_TypeError,
		   "Could not get _data property of Argument '%s'", objname);
      return NULL;
    }
    if (!PyArray_Check(dataobj)){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._data' must be a numpy array", objname);
      return NULL;
    }
    if (!(PyArray_ISCARRAY((PyArrayObject *) dataobj)
	  || PyArray_ISCARRAY_RO((PyArrayObject *) dataobj)) ){
      PyErr_Format(PyExc_TypeError,
		   "Argument '%s._data' is not C-order contiguous",objname);
      return NULL;
    }
    if ( PyArray_NDIM((PyArrayObject *) dataobj) != 1) {
      PyErr_Format(PyExc_ValueError,
		   "Argument '%s._data' is not one-dimensional",objname);
      return NULL;
    }
    if (PyArray_TYPE((PyArrayObject *) dataobj) != numpy_type) {
      PyErr_Format(PyExc_ValueError,
		   "Argument '%s._data' has wrong dtype for corresponding LAL vector",objname);
      return NULL;
    }

    switch(numpy_type){
    case NPY_FLOAT32:
      returnptr = calloc(1,sizeof(REAL4Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary REAL4Vector for argument '%s'",objname);
	return NULL;
      }
      ((REAL4Vector *) returnptr)->data = (REAL4 *) PyArray_DATA(dataobj);
      ((REAL4Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    case NPY_FLOAT64:
      returnptr = calloc(1,sizeof(REAL8Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary REAL8Vector for argument '%s'",objname);
	return NULL;
      }
      ((REAL8Vector *) returnptr)->data = (REAL8 *) PyArray_DATA(dataobj);
      ((REAL8Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    case NPY_COMPLEX64:
      returnptr = calloc(1,sizeof(COMPLEX8Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary COMPLEX8Vector for argument '%s'",objname);
	return NULL;
      }
      ((COMPLEX8Vector *) returnptr)->data = (COMPLEX8 *) PyArray_DATA(dataobj);
      ((COMPLEX8Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    case NPY_COMPLEX128:
      returnptr = calloc(1,sizeof(COMPLEX16Vector));
      if (!returnptr) {
	PyErr_Format(PyExc_MemoryError,
		     "Could not allocate temporary COMPLEX16Vector for argument '%s'",objname);
	return NULL;
      }
      ((COMPLEX16Vector *) returnptr)->data = (COMPLEX16 *) PyArray_DATA(dataobj);
      ((COMPLEX16Vector *) returnptr)->length = (UINT4) PyArray_DIM(dataobj,0);
      break;
    default:
      // If we get here, *we've* made a mistake in the wrapping---can't be user's fault.
      PyErr_Format(PyExc_TypeError,
		   "Argument array '%s._data' associated with wrong numpy type---this is a PyCBC bug",objname);
      return NULL;
    }
  return returnptr;
  }
}

%fragment("MarshallOutputTS","header",fragment="GenericTS") {
  PyObject *MarshallOutputVector(GenericTS *vect, const int numpy_type) {
    PyObject *result, *dataobj, *dtypeobj, *copybool, *constrdict;


    if (!(vect)) {
      PyErr_SetString(PyExc_ValueError,"Unexpected null time-series returned from function");
      return NULL;
    }
    if ( (vect->length) &&  !(vect->data) ) {
      PyErr_SetString(PyExc_ValueError,"Null data pointer returned for non-zero length array");
      return NULL;
    }

    constrdict = PyDict_New();
    if (!constrdict) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create dictionary for return value constructor");
      return NULL;
    }

    npy_intp dimensions[1];
    dimensions[0] = (npy_intp) vect->length;
    dataobj = PyArray_SimpleNewFromData(1,dimensions,numpy_type,(void *) vect->data);
    if (!dataobj) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create output data object");
      Py_DECREF(constrdict); // Dict still empty, so just delete
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"initial_array",dataobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add data object to cosntructor dict");
      Py_DECREF(constrdict); // Dict still empty, so just delete
      Py_DECREF(dataobj);
      return NULL;
    }

    dtypeobj = (PyObject *) PyArray_DescrNewFromType(numpy_type);
    if (!dtypeobj){
      PyErr_SetString(PyExc_RuntimeError,"Could not create output dtype object");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(dataobj);
      return NULL;
    }
    if (PyDict_SetItemString(constrdict,"dtype",dtypeobj)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add dtype object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(dataobj);
      Py_DECREF(dtypeobj);
      return NULL;
    }

    Py_INCREF(Py_False);
    copybool = Py_False;
    if (PyDict_SetItemString(constrdict,"copy",copybool)) {
      PyErr_SetString(PyExc_RuntimeError,"Could not add copy object to constructor dict");
      PyDict_Clear(constrdict);
      Py_DECREF(constrdict);
      Py_DECREF(dataobj);
      Py_DECREF(dtypeobj);
      Py_DECREF(copybool);
      return NULL;
    }


    result = PyObject_Call(CBC_Arr,EmptyTuple,constrdict);
    if (!result) {
      PyErr_SetString(PyExc_RuntimeError,"Could not create new instance of pycbc.types.Array");
    }
    // We don't need to do anything else for that last failure, as we'll be returning NULL
    // anyway and we have to do the same cleanup
    PyDict_Clear(constrdict);
    Py_DECREF(constrdict);
    Py_DECREF(dataobj);
    Py_DECREF(dtypeobj);
    Py_DECREF(copybool);

    return result;
  }
 }



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
this vector may in fact be treated by the XLAL function as "output", with its contents
modified in place. The elements of _data could have been modified, or also any of its properties
that can be set by the constructor.

The NEWOUT typemaps are for when a function returns a newly-allocated vector
of that type. This function will be wrapped into one which returns a newly allocated
instance of pycbc.types.{Array,TimeSeries,FrequencySeries}, with its '_data' element set to the
appropriate numpy type, and avoiding any memory copies. Other properties will be set as appropriate
to that specific array-like type (e.g. epoch, delta_t, etc).

The NONEOUT typemap is for functions that return a vector, but that vector is in
fact the same as one of the function's input vectors (which has been modified in
place).  In this case the wrapped function will always return "None" on successful
completion, and the effect of the function should be accessed through the
appropriate input input numpy array.

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
In particular any function that calls realloc() on a vector's data (including
functions that "cut" or "resize" vectors) can change the location of the data
pointer.  There is no way to accomodate this into Numpy's memory model, so you
should ensure that you never SWIG-wrap (with the typemaps of this file) any
function that may change the location in memory of the associated vector struct's
data. The best you might expect is a difficult-to-reproduce segfault.

YOU HAVE BEEN WARNED!

*/


// Typemaps for REAL4Vectors:

%typemap(in, fragment="MarshallInputVector") REAL4Vector *INPUT_REAL4V {
  $1 =(REAL4Vector *) MarshallInputVector($input,NPY_FLOAT32,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(freearg) REAL4Vector *INPUT_REAL4V {
  if ($1) {
    free((REAL4Vector *) $1);
  }
 }

%typemap(out, fragment="MarshallOutputVector") REAL4Vector *NEWOUT_REAL4V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_FLOAT32);
  if (!($result)) SWIG_fail;
 }

%typemap(newfree) REAL4Vector *NEWOUT_REAL4V{
  free( (REAL4Vector *) $1);
 }

%typemap(out) REAL4Vector *IGNOUT_REAL4V{
  Py_INCREF(Py_None);
  $result = Py_None;
 }

// Typemaps for REAL8Vectors:

%typemap(in, fragment="MarshallInputVector") REAL8Vector *INPUT_REAL8V {
  $1 =(REAL8Vector *) MarshallInputVector($input,NPY_FLOAT64,"$1_name");
  if (!($1)) SWIG_fail;
 }

%typemap(freearg) REAL8Vector *INPUT_REAL8V {
  if ($1) {
    free((REAL8Vector *) $1);
  }
 }

%typemap(out, fragment="MarshallOutputVector") REAL8Vector *NEWOUT_REAL8V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_FLOAT64);
  if (!($result)) SWIG_fail;
}

%typemap(newfree) REAL8Vector *NEWOUT_REAL8V{
  free( (REAL8Vector *) $1);
 }

%typemap(out) REAL8Vector *IGNOUT_REAL8V{
  Py_INCREF(Py_None);
  $result = Py_None;
}

// Typemaps for COMPLEX8Vectors:

%typemap(in, fragment="MarshallInputVector") COMPLEX8Vector *INPUT_COMPLEX8V {
  $1 =(COMPLEX8Vector *) MarshallInputVector($input,NPY_COMPLEX64,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(freearg) COMPLEX8Vector *INPUT_COMPLEX8V {
  if ($1) {
    free((COMPLEX8Vector *) $1);
  }
}

%typemap(out, fragment="MarshallOutputVector") COMPLEX8Vector *NEWOUT_COMPLEX8V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_COMPLEX64);
  if (!($result)) SWIG_fail;
 }

%typemap(newfree) COMPLEX8Vector *NEWOUT_COMPLEX8V{
  free( (COMPLEX8Vector *) $1);
 }

%typemap(out) COMPLEX8Vector *IGNOUT_COMPLEX8V{
  Py_INCREF(Py_None);
  $result = Py_None;
 }

// Typemaps for COMPLEX16Vectors:

%typemap(in, fragment="MarshallInputVector") COMPLEX16Vector *INPUT_COMPLEX16V {
  $1 =(COMPLEX8Vector *) MarshallInputVector($input,NPY_COMPLEX128,"$1_name");
  if (!($1)) SWIG_fail;
}

%typemap(freearg) COMPLEX16Vector *INPUT_COMPLEX16V {
  if ($1) {
    free((COMPLEX16Vector *) $1);
  }
}

%typemap(out, fragment="MarshallOutputVector") COMPLEX16Vector *NEWOUT_COMPLEX16V{
  $result = MarshallOutputVector((GenericVector *) $1,NPY_COMPLEX128);
  if (!($result)) SWIG_fail;
 }

%typemap(newfree) COMPLEX16Vector *NEWOUT_COMPLEX16V{
  free( (COMPLEX16Vector *) $1);
 }

%typemap(out) COMPLEX16Vector *IGNOUT_COMPLEX16V{
  Py_INCREF(Py_None);
  $result = Py_None;
 }

// Typemaps for REAL4 Time Series



// Some tests:
%rename("%s") XLALSSVectorMultiply;
%apply REAL4Vector *IGNOUT_REAL4V { REAL4Vector *XLALSSVectorMultiply };
%apply REAL4Vector *INPUT_REAL4V {REAL4Vector *};
extern REAL4Vector *XLALSSVectorMultiply(REAL4Vector *out, REAL4Vector *in1, REAL4Vector *in2);

%rename("%s") XLALCreateREAL4Vector;
%apply REAL4Vector *NEWOUT_REAL4V {REAL4Vector *XLALCreateREAL4Vector};
extern REAL4Vector *XLALCreateREAL4Vector(UINT4 length);


%rename("%s") XLALCOMPLEX8VectorFFT;
%apply COMPLEX8Vector *INPUT_COMPLEX8V {COMPLEX8Vector *};
extern int XLALCOMPLEX8VectorFFT(COMPLEX8Vector *output, COMPLEX8Vector *input, const COMPLEX8FFTPlan *plan );

