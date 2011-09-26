/* opencmiss_wrapper.i */
%module opencmiss
%{
#include "stdlib.h"
#include "opencmiss.h"
#define MAX_OUTPUT_STRING_SIZE 200
%}

/* TODO: Handle output of arrays */

/* Macros for returning a python object in a tuple of return values */
%define RESULT_VARS()
  PyObject *new_output_tuple;
  PyObject *previous_result;
%enddef
%define APPEND_TO_RESULT(new_result)
  if ((!$result) || ($result == Py_None)) {
    /* If there isn't already a result, just return this */
    $result = new_result;
  } else {
    if (!PyTuple_Check($result)) {
      /* If the previous result isn't a tuple, create a tuple containing the result */
      PyObject *previous_result = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result,0,previous_result);
    }
    /* Add result from this argument to the tuple */
    new_output_tuple = PyTuple_New(1);
    PyTuple_SetItem(new_output_tuple,0,new_result);
    previous_result = $result; /* previous tuple result */
    $result = PySequence_Concat(previous_result,new_output_tuple);
    Py_DECREF(previous_result);
    Py_DECREF(new_output_tuple);
  }
%enddef

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
  RESULT_VARS()
  PyObject *output_pointer;

  output_pointer = SWIG_NewPointerObj(*$1, $*1_descriptor, 0);

  APPEND_TO_RESULT(output_pointer)
}

/* String input */
%typemap(in,numinputs=1) (const int Size, const char *DummyInputString) {
  if (!PyString_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a string");
    return NULL;
  }
  $1 = PyString_Size($input);
  $2 = PyString_AsString($input);
}

/* String output arguments */
%typemap(in,numinputs=0) (const int Size, char *DummyOutputString)(char temp[MAX_OUTPUT_STRING_SIZE]) {
  $1 = MAX_OUTPUT_STRING_SIZE;
  $2 = &temp[0];
}
%typemap(argout) (const int Size, char *DummyOutputString) {
  RESULT_VARS
  PyObject *output_string;

  output_string = PyString_FromString($2);
  APPEND_TO_RESULT(output_string)
}

/* Integer output: */
%typemap(in,numinputs=0) (int *DummyOutputInt)(int temp) {
  $1 = &temp;
}
%typemap(argout) (int *DummyOutputInt) {
  RESULT_VARS
  PyObject *output_int;

  output_int = PyInt_FromLong((long) *$1);
  APPEND_TO_RESULT(output_int)
}

/* Float output: */
%typemap(in,numinputs=0) (double *DummyOutputDouble)(double temp) {
  $1 = &temp;
}
%typemap(argout) (double *DummyOutputDouble) {
  RESULT_VARS
  PyObject *output_double;

  output_double = PyFloat_FromDouble((double) *$1);
  APPEND_TO_RESULT(output_double)
}

%typemap(in,numinputs=0) (float *DummyOutputFloat)(float temp) {
  $1 = &temp;
}
%typemap(argout) (float *DummyOutputFloat) {
  RESULT_VARS
  PyObject *output_double;

  output_double = PyFloat_FromDouble((double) *$1);
  APPEND_TO_RESULT(output_double)
}

/* Boolean input: */
%typemap(in) (const int DummyInputBool) {
  $1 = PyObject_IsTrue($input);
}

/* Boolean output: */
%typemap(in,numinputs=0) (int *DummyOutputBool)(int temp) {
  $1 = &temp;
}
%typemap(argout) (int *DummyOutputBool) {
  RESULT_VARS
  PyObject *output_bool;

  output_bool = PyBool_FromLong((long) *$1);
  APPEND_TO_RESULT(output_bool)
}

%include "opencmiss.h"
