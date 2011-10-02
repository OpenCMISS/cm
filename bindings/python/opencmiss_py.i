/* OpenCMISS SWIG interface for Python */
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

/* Macro for array input */
%define ARRAY_INPUT(sequence_type, check_routine, convert_routine, readable_sequence_type)
%typemap(in,numinputs=1) (const int ArraySize, const sequence_type *DummyInputArray)(int len, int i) {
  PyObject *o;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a sequence");
    return NULL;
  }
  len = PyObject_Length($input);
  $2 = (sequence_type *) malloc(len * sizeof(sequence_type));
  if ($2 == NULL) {
    PyErr_SetString(PyExc_MemoryError,"Could not allocate memory for array");
    return NULL;
  } else {
    for (i=0; i < len; i++) {
      o = PySequence_GetItem($input,i);
      if (!check_routine(o)) {
        Py_XDECREF(o);
        PyErr_SetString(PyExc_ValueError,"Expected a sequence of readable_sequence_type");
        free($2);
        return NULL;
      }
      *($2 + i) = (sequence_type) convert_routine(o);
      Py_DECREF(o);
    }
  }
  $1 = len;
}
%typemap(freearg) (const int ArraySize, const sequence_type *DummyInputArray) {
    free($2);
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
  $1 = PyString_Size($input)+1;
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

/* Array input */
ARRAY_INPUT(int, PyInt_Check, PyInt_AsLong, integers)

ARRAY_INPUT(double, PyFloat_Check, PyFloat_AsDouble, floats)

ARRAY_INPUT(float, PyFloat_Check, PyFloat_AsDouble, floats)

/* Array of CMISS types */
%typemap(in,numinputs=1) (const int ArraySize, const CMISSDummyType *DummyTypes)(int len, int i) {
  PyObject *o;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a sequence");
    return NULL;
  }
  len = PyObject_Length($input);
  $2 = ($2_ltype) malloc(len * sizeof($*2_ltype));
  if ($2 == NULL) {
    PyErr_SetString(PyExc_MemoryError,"Could not allocate memory for array");
    return NULL;
  } else {
    for (i=0; i < len; i++) {
      o = PySequence_GetItem($input,i);
      if (SWIG_ConvertPtr(o, (void **) ($2+i), $*2_descriptor, SWIG_POINTER_EXCEPTION) == -1) {
        PyErr_SetString(PyExc_TypeError,"Expected a sequence of CMISS types.");
        free($2);
        return NULL;
      }
      Py_DECREF(o);
    }
  }
  $1 = len;
}
%typemap(freearg) (const int ArraySize, const CMISSDummyType *DummyTypes) {
    free($2);
}

/* Input array of strings */
%typemap(in,numinputs=1) (const int NumStrings, const int StringLength, const char *DummyStringList)(int len, int i, Py_ssize_t max_strlen) {
  PyObject *o;
  max_strlen = 0;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a sequence");
    return NULL;
  }
  len = PyObject_Length($input);
  for (i =0; i < len; i++) {
    o = PySequence_GetItem($input,i);
    if (!PyString_Check(o)) {
      Py_XDECREF(o);
      PyErr_SetString(PyExc_ValueError,"Expected a sequence of strings");
      return NULL;
    }
    if (PyString_Size(o) > max_strlen) {
      max_strlen = PyString_Size(o);
    }
    Py_DECREF(o);
  }
  max_strlen = max_strlen + 1; /* Null terminator */
  $3 = (char *) malloc(len * max_strlen * sizeof(char));
  if ($3 == NULL) {
    PyErr_SetString(PyExc_MemoryError,"Could not allocate memory for array");
    return NULL;
  } else {
    for (i=0; i < len; i++) {
      o = PySequence_GetItem($input,i);
      strncpy($3+i*max_strlen, PyString_AsString(o), PyString_Size(o)+1);
      Py_DECREF(o);
    }
  }
  $1 = len;
  $2 = (int) max_strlen;
}
%typemap(freearg) (const int NumStrings, const int StringLength, const char *DummyStringList) {
    free($3);
}

%include "opencmiss.h"
