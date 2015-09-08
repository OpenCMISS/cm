/* Python specific typemaps for SWIG */
%module iron_python
%{
#include "stdlib.h"
#include "opencmiss.h"
#define MAX_OUTPUT_STRING_SIZE 300
%}

/* Use numpy for input/output of arrays
   This means that arrays are contiguous
   memory chunks so we don't need to do
   expensive copying to/from Python lists
*/
%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%include "numpy_extra.i"
%init %{
  import_array();
%}

/* Apply numpy macros using const int for the dim type */
%numpy_typemaps(unsigned char , NPY_UBYTE , const int)
%numpy_typemaps(short , NPY_SHORT , const int)
%numpy_typemaps(unsigned short , NPY_USHORT , const int)
%numpy_typemaps(int , NPY_INT , const int)
%numpy_typemaps(unsigned int , NPY_UINT , const int)
%numpy_typemaps(long , NPY_LONG , const int)
%numpy_typemaps(unsigned long , NPY_ULONG , const int)
%numpy_typemaps(long long , NPY_LONGLONG , const int)
%numpy_typemaps(unsigned long long, NPY_ULONGLONG, const int)
%numpy_typemaps(float , NPY_FLOAT , const int)
%numpy_typemaps(double , NPY_DOUBLE , const int)
/* And extra typemaps defined for OpenCMISS */
%numpy_extra_typemaps(unsigned char , NPY_UBYTE , const int)
%numpy_extra_typemaps(short , NPY_SHORT , const int)
%numpy_extra_typemaps(unsigned short , NPY_USHORT , const int)
%numpy_extra_typemaps(int , NPY_INT , const int)
%numpy_extra_typemaps(unsigned int , NPY_UINT , const int)
%numpy_extra_typemaps(long , NPY_LONG , const int)
%numpy_extra_typemaps(unsigned long , NPY_ULONG , const int)
%numpy_extra_typemaps(long long , NPY_LONGLONG , const int)
%numpy_extra_typemaps(unsigned long long, NPY_ULONGLONG, const int)
%numpy_extra_typemaps(float , NPY_FLOAT , const int)
%numpy_extra_typemaps(double , NPY_DOUBLE , const int)

/**** Macros ****/

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
  PyObject *output_pointer;

  output_pointer = SWIG_NewPointerObj(*$1, $*1_descriptor, 0);

  $result = SWIG_Python_AppendOutput($result,output_pointer);
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
  PyObject *output_string;

  output_string = PyString_FromString($2);
  $result = SWIG_Python_AppendOutput($result,output_string);
}

/**** Scalars ****/

/* Integer output */
%typemap(in,numinputs=0) (int *DummyOutputScalar)(int temp) {
  $1 = &temp;
}
%typemap(argout) (int *DummyOutputScalar) {
  PyObject *output_int;

  output_int = PyInt_FromLong((long) *$1);
  $result = SWIG_Python_AppendOutput($result,output_int);
}

/* Float output */
%typemap(in,numinputs=0) (double *DummyOutputScalar)(double temp) {
  $1 = &temp;
}
%typemap(argout) (double *DummyOutputScalar) {
  PyObject *output_double;

  output_double = PyFloat_FromDouble((double) *$1);
  $result = SWIG_Python_AppendOutput($result,output_double);
}

%typemap(in,numinputs=0) (float *DummyOutputScalar)(float temp) {
  $1 = &temp;
}
%typemap(argout) (float *DummyOutputScalar) {
  PyObject *output_double;

  output_double = PyFloat_FromDouble((double) *$1);
  $result = SWIG_Python_AppendOutput($result,output_double);
}

/* Boolean input */
%typemap(in) (const CMISSBool DummyInputBool) {
  $1 = PyObject_IsTrue($input);
}

/* Boolean output */
%typemap(in,numinputs=0) (CMISSBool *DummyOutputScalar)(int temp) {
  $1 = &temp;
}
%typemap(argout) (CMISSBool *DummyOutputScalar) {
  PyObject *output_bool;

  output_bool = PyBool_FromLong((long) *$1);
  $result = SWIG_Python_AppendOutput($result,output_bool);
}

/**** Arrays ****/

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

%include "iron_generated.i"
