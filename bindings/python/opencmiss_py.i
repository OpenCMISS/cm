/* Python specific typemaps for SWIG */
%module opencmiss
%{
#include "stdlib.h"
#include "opencmiss.h"
#define MAX_OUTPUT_STRING_SIZE 200
%}

/**** Macros ****/

/* Return a python object in a tuple of return values
 */
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

/* Array input
 */
%define ARRAY_INPUT(sequence_type, check_routine, convert_routine, readable_sequence_type)
%typemap(in,numinputs=1) (const int ArraySize, const sequence_type *DummyInputArray)(int len, int i, PyObject *o) {
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

/* 2D array input
 */
%define ARRAY_INPUT_2D(sequence_type, check_routine, convert_routine, readable_sequence_type)
%typemap(in,numinputs=1) (const int ArraySize1, const int ArraySize2, const sequence_type *DummyInputArray)(int i, int j, PyObject *oi, PyObject *oj) {
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a sequence");
    return NULL;
  }
  $1 = PyObject_Length($input);
  for (i=0; i < $1; i++) {
    oi = PySequence_GetItem($input,i);
    if (!PySequence_Check(oi)) {
      Py_DECREF(oi);
      PyErr_SetString(PyExc_TypeError,"Expected a sequence of sequences");
      return NULL;
    }
    if (i == 0) {
      $2 = PyObject_Length(oi);
    }
    else if (PyObject_Length(oi) != $2) {
      Py_DECREF(oi);
      PyErr_SetString(PyExc_TypeError,"Sub-lists must be of equal length");
      return NULL;
    }
    Py_DECREF(oi);
  }
  $3 = (sequence_type *) malloc($1 * $2 * sizeof(sequence_type));
  if ($3 == NULL) {
    PyErr_SetString(PyExc_MemoryError,"Could not allocate memory for array");
    return NULL;
  } else {
    for (i=0; i < $1; i++) {
      oi = PySequence_GetItem($input,i);
      for (j=0; j < $2; j++) {
        oj = PySequence_GetItem(oi,j);
        if (!check_routine(oj)) {
          Py_XDECREF(oj);
          PyErr_SetString(PyExc_ValueError,"Expected a 2D sequence of readable_sequence_type");
          free($3);
          return NULL;
        }
        *($3 + j*$1 + i) = (sequence_type) convert_routine(oj);
        Py_DECREF(oj);
      }
      Py_DECREF(oi);
    }
  }
}
%typemap(freearg) (const int ArraySize2, const int ArraySize2, const sequence_type *DummyInputArray) {
    free($3);
}
%enddef

/* Array output, with expected size passed as a parameter
 * Already allocated memory should be passed in, so require size parameter from Python
 */
%define ARRAY_OUTPUT(sequence_type,convert_call,convert_type)
%typemap(in,numinputs=1) (const int ArraySize, sequence_type *DummyOutputArray) {
  if (!PyInt_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected an integer array size");
    return NULL;
  }
  $1 = PyInt_AsLong($input);
  $2 = (sequence_type *) malloc($1 * sizeof(sequence_type));
}
%typemap(argout) (const int ArraySize, sequence_type *DummyOutputArray)(int i, PyObject *output_list, PyObject *list_item) {
  PyObject *new_output_tuple;
  PyObject *previous_result;

  output_list = PyList_New($1);
  for (i=0; i<$1 ;i++) {
    list_item = convert_call((convert_type) *($2 + i));
    PyList_SetItem(output_list, i, list_item);
  }

  if ((!$result) || ($result == Py_None)) {
    /* If there isn't already a result, just return this */
    $result = output_list;
  } else {
    if (!PyTuple_Check($result)) {
      /* If the previous result isn't a tuple, create a tuple containing the result */
      PyObject *previous_result = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result,0,previous_result);
    }
    /* Add result from this argument to the tuple */
    new_output_tuple = PyTuple_New(1);
    PyTuple_SetItem(new_output_tuple,0,output_list);
    previous_result = $result; /* previous tuple result */
    $result = PySequence_Concat(previous_result,new_output_tuple);
    Py_DECREF(previous_result);
    Py_DECREF(new_output_tuple);
  }
}
%typemap(freearg) (const int ArraySize, sequence_type *DummyOutputArray) {
  free($2);
}
%enddef

/* Array output of 2D arrays, with expected size passed as a parameter
 * SWIG only allows one Python input for each typemap, so we have to pass all
 * array dimensions through a tuple
 */
%define ARRAY_OUTPUT_2D(sequence_type,convert_call,convert_type)
%typemap(in,numinputs=1) (const int ArraySize1, const int ArraySize2, sequence_type *DummyOutputArray) {
  if (!PyTuple_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple of array dimensions");
    return NULL;
  }
  if (!PyArg_ParseTuple($input,"ii",&$1,&$2)) {
    PyErr_SetString(PyExc_TypeError,"Size tuple must contain two integers");
    return NULL;
  }
  $3 = (sequence_type *) malloc($1 * $2 * sizeof(sequence_type));
}
%typemap(argout) (const int ArraySize1, const int ArraySize2, sequence_type *DummyOutputArray)(
    int i, int j, PyObject *output_list, PyObject *sub_list, PyObject *list_item, PyObject *new_output_tuple, PyObject *previous_result) {

  output_list = PyList_New($1);
  for (i=0; i<$1; i++) {
    sub_list = PyList_New($2);
    PyList_SetItem(output_list, i, sub_list);
    for (j=0; j<$2; j++) {
      /* Fortran stores arrays by column */
      list_item = convert_call((convert_type) *($3 + j*$1 + i));
      PyList_SetItem(sub_list, j, list_item);
    }
  }

  if ((!$result) || ($result == Py_None)) {
    /* If there isn't already a result, just return this */
    $result = output_list;
  } else {
    if (!PyTuple_Check($result)) {
      /* If the previous result isn't a tuple, create a tuple containing the result */
      PyObject *previous_result = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result,0,previous_result);
    }
    /* Add result from this argument to the tuple */
    new_output_tuple = PyTuple_New(1);
    PyTuple_SetItem(new_output_tuple,0,output_list);
    previous_result = $result; /* previous tuple result */
    $result = PySequence_Concat(previous_result,new_output_tuple);
    Py_DECREF(previous_result);
    Py_DECREF(new_output_tuple);
  }
}
%typemap(freearg) (const int ArraySize1, const int ArraySize2, sequence_type *DummyOutputArray) {
  free($3);
}
%enddef

/* Array output of a known size, so the size isn't passed as a
 * parameter from Python
 */
%define ARRAY_OUTPUT_WITH_SIZE(sequence_type,convert_call,convert_type,array_size,array_name)
%typemap(in,numinputs=0) (sequence_type *array_name) {
  $1 = (sequence_type *) malloc(array_size * sizeof(sequence_type));
}
%typemap(argout) (sequence_type *array_name)(int i, PyObject *output_list, PyObject *list_item) {
  PyObject *new_output_tuple;
  PyObject *previous_result;

  output_list = PyList_New(array_size);
  for (i=0; i<array_size ;i++) {
    list_item = convert_call((convert_type) *($1 + i));
    PyList_SetItem(output_list, i, list_item);
  }

  if ((!$result) || ($result == Py_None)) {
    /* If there isn't already a result, just return this */
    $result = output_list;
  } else {
    if (!PyTuple_Check($result)) {
      /* If the previous result isn't a tuple, create a tuple containing the result */
      PyObject *previous_result = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result,0,previous_result);
    }
    /* Add result from this argument to the tuple */
    new_output_tuple = PyTuple_New(1);
    PyTuple_SetItem(new_output_tuple,0,output_list);
    previous_result = $result; /* previous tuple result */
    $result = PySequence_Concat(previous_result,new_output_tuple);
    Py_DECREF(previous_result);
    Py_DECREF(new_output_tuple);
  }
}
%typemap(freearg) (sequence_type *array_name) {
  free($1);
}
%enddef

/**** CMISS*Type typemaps ****/

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

/* CMISS*TypeFinalise routines. Convert SWIG pointer input. Can't modify the input pointer to nullify it though */
%typemap(in) CMISSDummyFinaliseType *CMISSDummy ($*1_ltype type_pointer) {
  if (SWIG_ConvertPtr($input, (void **) (&type_pointer), $*1_descriptor, SWIG_POINTER_EXCEPTION) == -1) {
    PyErr_SetString(PyExc_TypeError,"Input must be a SWIG pointer to the correct CMISS type.");
    return NULL;
  }
  $1 = &type_pointer;
}

/**** Strings ****/

/* String input */
%typemap(in,numinputs=1) (const int Size, const char *DummyInputString) {
  if (!PyString_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a string");
    return NULL;
  }
  $1 = PyString_Size($input)+1;
  $2 = PyString_AsString($input);
}

/* String output */
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

/**** Scalars ****/

/* Integer output */
%typemap(in,numinputs=0) (int *DummyOutputInt)(int temp) {
  $1 = &temp;
}
%typemap(argout) (int *DummyOutputInt) {
  RESULT_VARS
  PyObject *output_int;

  output_int = PyInt_FromLong((long) *$1);
  APPEND_TO_RESULT(output_int)
}

/* Float output */
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

/* Boolean input */
%typemap(in) (const int DummyInputBool) {
  $1 = PyObject_IsTrue($input);
}

/* Boolean output */
%typemap(in,numinputs=0) (int *DummyOutputBool)(int temp) {
  $1 = &temp;
}
%typemap(argout) (int *DummyOutputBool) {
  RESULT_VARS
  PyObject *output_bool;

  output_bool = PyBool_FromLong((long) *$1);
  APPEND_TO_RESULT(output_bool)
}

/**** Arrays ****/

/* Array input */
ARRAY_INPUT(int, PyInt_Check, PyInt_AsLong, integers)
ARRAY_INPUT(double, PyFloat_Check, PyFloat_AsDouble, floats)
ARRAY_INPUT(float, PyFloat_Check, PyFloat_AsDouble, floats)

/* 2D array input */
ARRAY_INPUT_2D(int, PyInt_Check, PyInt_AsLong, integers)
ARRAY_INPUT_2D(float, PyFloat_Check, PyFloat_AsDouble, floats)
ARRAY_INPUT_2D(double, PyFloat_Check, PyFloat_AsDouble, floats)

/* Array output */
ARRAY_OUTPUT(int, PyInt_FromLong, int)
ARRAY_OUTPUT(double, PyFloat_FromDouble, double)
ARRAY_OUTPUT(float, PyFloat_FromDouble, double)

/* Array output with known sizes, need to add more sizes as required */
ARRAY_OUTPUT_WITH_SIZE(double, PyFloat_FromDouble, double, 2, DummyOutputArraySize2)
ARRAY_OUTPUT_WITH_SIZE(double, PyFloat_FromDouble, double, 3, DummyOutputArraySize3)
ARRAY_OUTPUT_WITH_SIZE(double, PyFloat_FromDouble, double, 8, DummyOutputArraySize8)

/* Array output of 2d arrays */
ARRAY_OUTPUT_2D(int, PyInt_FromLong, int)
ARRAY_OUTPUT_2D(float, PyFloat_FromDouble, double)
ARRAY_OUTPUT_2D(double, PyFloat_FromDouble, double)

/* Array of CMISS types */
%typemap(in,numinputs=1) (const int ArraySize, const CMISSDummyType *DummyTypes)(int len, int i, PyObject *o) {
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
%typemap(in,numinputs=1) (const int NumStrings, const int StringLength, const char *DummyStringList)(int len, int i, Py_ssize_t max_strlen, PyObject *o) {
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

%include "opencmiss.i"
