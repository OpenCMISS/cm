/* opencmiss_wrapper.i */
%module opencmiss
%{
#include "stdlib.h"
#include "opencmiss.h"
%}

/* Typemaps for passing CMISS types to CMISS...Type initialise routines
   We don't need to pass an input value, we can create a NULL
   pointer within the C wrapper */
%typemap(in,numinputs=0) CMISSDummyInitialiseType *CMISSDummy($*1_ltype temp) {
  temp = ($*1_ltype)NULL;
  $1 = &temp;
}

/* Typemap to convert the output pointer to a SWIG pointer we can then
   pass it into other routines from Python */
%typemap(argout) CMISSDummyInitialiseType *CMISSDummy {
  PyObject *new_output_tuple;
  PyObject *previous_result;
  PyObject *output_pointer;

  output_pointer = SWIG_NewPointerObj(*$1, $*1_descriptor, 1);

  if ((!$result) || ($result == Py_None)) {
    /* If there isn't already a result, just return this */
    $result = output_pointer;
  } else {
    if (!PyTuple_Check($result)) {
      /* If the previous result isn't a tuple, create a tuple containing the result */
      PyObject *previous_result = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result,0,previous_result);
    }
    /* Add result from this argument to the tuple */
    new_output_tuple = PyTuple_New(1);
    PyTuple_SetItem(new_output_tuple,0,output_pointer);
    previous_result = $result; /* previous tuple result */
    $result = PySequence_Concat(previous_result,new_output_tuple);
    Py_DECREF(previous_result);
    Py_DECREF(new_output_tuple);
  }
}

%include "opencmiss.h"
