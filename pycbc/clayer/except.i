%include exception.i 

%{
#include <pycbc/clayer/except.h>
%}

%exception {
  int num;
  char *err;

  pycbc_clear_exception();

  $action

  if ( (num = pycbc_check_exception()) )
  {
     err = pycbc_get_error_message();
     switch (num)
     {
       case PYCBC_ATTRIBUTE_ERROR:
         PyErr_SetString(PyExc_AttributeError, err);
         break;
       case PYCBC_EOF_ERROR:
         PyErr_SetString(PyExc_EOFError, err);
         break;
       case PYCBC_IO_ERROR:
         PyErr_SetString(PyExc_IOError, err);
         break;
       case PYCBC_INDEX_ERROR:
         PyErr_SetString(PyExc_IndexError, err);
         break;
       case PYCBC_TYPE_ERROR:
         PyErr_SetString(PyExc_TypeError, err);
         break;
       case PYCBC_VALUE_ERROR:
         PyErr_SetString(PyExc_ValueError, err);
         break;
       case PYCBC_MEMORY_ERROR:
         PyErr_SetString(PyExc_MemoryError, err);
         break;
       case PYCBC_NAME_ERROR:
         PyErr_SetString(PyExc_NameError, err);
         break;
       case PYCBC_OVERFLOW_ERROR:
         PyErr_SetString(PyExc_OverflowError, err);
         break;
       case PYCBC_ZERO_DIVISION_ERROR:
         PyErr_SetString(PyExc_ZeroDivisionError, err);
         break;
       case PYCBC_RUNTIME_ERROR:
         PyErr_SetString(PyExc_RuntimeError, err);
         break;
       case PYCBC_STOP_ITERATION:
         PyErr_SetNone(PyExc_StopIteration);
         break;
       default:
         PyErr_SetString(PyExc_RuntimeError, err);
     }
     return NULL;
  }
}
