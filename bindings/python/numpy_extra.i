/* -*- C -*-  (not really, but good for syntax highlighting) */

/* Additional numpy typemaps required by OpenCMISS */

#ifdef SWIGPYTHON

%define %numpy_extra_typemaps(DATA_TYPE, DATA_TYPECODE, DIM_TYPE)
/* Typemap suite for (DATA_TYPE* ARGOUT_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
 */
%typemap(in,numinputs=1,
         fragment="NumPy_Fragments")
  (DATA_TYPE* ARGOUT_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
  (PyObject * array = NULL, PyObject * tupleItem, int tupleSize)
{
  if (!PyTuple_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple of array dimensions");
    SWIG_fail;
  }
  tupleSize = PyTuple_Size($input);
  if (tupleSize != 2) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple with length 2");
    SWIG_fail;
  }
  /* tupleItem is a borrowed reference, so don't call Py_DECREF */
  tupleItem = PyTuple_GetItem($input, 0);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $2 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  tupleItem = PyTuple_GetItem($input, 1);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $3 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  npy_intp dims[2] = { $2, $3 };
  array = PyArray_SimpleNew(2, dims, DATA_TYPECODE);
  if (!array) SWIG_fail;
  $1 = ($1_ltype) array_data(array);
}
%typemap(argout)
  (DATA_TYPE* ARGOUT_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
{
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* ARGOUT_ARRAY2)
 */
%typemap(in,numinputs=1,
         fragment="NumPy_Fragments")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* ARGOUT_ARRAY2)
  (PyObject * array = NULL, PyObject * tupleItem, int tupleSize)
{
  if (!PyTuple_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple of array dimensions");
    SWIG_fail;
  }
  tupleSize = PyTuple_Size($input);
  if (tupleSize != 2) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple with length 2");
    SWIG_fail;
  }
  /* tupleItem is a borrowed reference, so don't call Py_DECREF */
  tupleItem = PyTuple_GetItem($input, 0);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $1 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  tupleItem = PyTuple_GetItem($input, 1);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $2 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  npy_intp dims[2] = { $1, $2 };
  array = PyArray_SimpleNew(2, dims, DATA_TYPECODE);
  if (!array) SWIG_fail;
  $3 = ($3_ltype) array_data(array);
}
%typemap(argout)
  (DATA_TYPE* ARGOUT_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
{
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

/* Typemap suite for (DATA_TYPE* ARGOUT_FARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
 */
%typemap(in,numinputs=1,
         fragment="NumPy_Fragments")
  (DATA_TYPE* ARGOUT_FARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
  (PyObject * array = NULL, PyObject * tupleItem, int tupleSize)
{
  if (!PyTuple_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple of array dimensions");
    SWIG_fail;
  }
  tupleSize = PyTuple_Size($input);
  if (tupleSize != 2) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple with length 2");
    SWIG_fail;
  }
  /* tupleItem is a borrowed reference, so don't call Py_DECREF */
  tupleItem = PyTuple_GetItem($input, 0);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $2 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  tupleItem = PyTuple_GetItem($input, 1);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $3 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  npy_intp dims[2] = { $2, $3 };
  /* Last parameter means Fortran ordering */
  array = PyArray_EMPTY(2, dims, DATA_TYPECODE, 1);
  if (!array) SWIG_fail;
  $1 = ($1_ltype) array_data(array);
}
%typemap(argout)
  (DATA_TYPE* ARGOUT_FARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
  (PyArrayObject * arrayObj)
{
  arrayObj = obj_to_array_no_conversion(array$argnum, DATA_TYPECODE);
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* ARGOUT_FARRAY2)
 */
%typemap(in,numinputs=1,
         fragment="NumPy_Fragments")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* ARGOUT_FARRAY2)
  (PyObject * array = NULL, PyObject * tupleItem, int tupleSize)
{
  if (!PyTuple_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple of array dimensions");
    SWIG_fail;
  }
  tupleSize = PyTuple_Size($input);
  if (tupleSize != 2) {
    PyErr_SetString(PyExc_TypeError,"Expected a tuple with length 2");
    SWIG_fail;
  }
  /* tupleItem is a borrowed reference, so don't call Py_DECREF */
  tupleItem = PyTuple_GetItem($input, 0);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $1 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  tupleItem = PyTuple_GetItem($input, 1);
  if (!PyInt_Check(tupleItem))
  {
    const char* typestring = pytype_string(tupleItem);
    PyErr_Format(PyExc_TypeError,
                 "Int dimension expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  $2 = (DIM_TYPE) PyInt_AsLong(tupleItem);

  npy_intp dims[2] = { $1, $2 };
  /* Last parameter means Fortran ordering */
  array = PyArray_EMPTY(2, dims, DATA_TYPECODE, 1);
  if (!array) SWIG_fail;
  $3 = ($3_ltype) array_data(array);
}
%typemap(argout)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* ARGOUT_FARRAY2)
  (PyArrayObject * arrayObj)
{
  arrayObj = obj_to_array_no_conversion(array$argnum, DATA_TYPECODE);
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

/* Typemap suite for (DIM_TYPE* DIM1, DATA_TYPE** ARGINOUTVIEW_ARRAY1)
 *
 * This typemap is for cases where a pointer to an array is required
 * as an input and is also an output. This is used when the
 * input pointer points to memory allocated internally and is modified,
 * eg. for CMISSField_ParameterSetDataRestore where the pointer is null
 * on output.
 * Using it with data allocated from the Python side can result
 * in memory leaks.
 * The input typemap requires a NumPy array, it can't be any
 * Python iterable as we need to modify the data and dimensions
 * on return.
 */

%typemap(in,numinputs=1)
  (DIM_TYPE* DIM1, DATA_TYPE** ARGINOUTVIEW_ARRAY1)
  (DIM_TYPE dim1_temp, DATA_TYPE* data_temp, PyArrayObject * tempArray)
{
  /* Check input is an array */
  if (!is_array($input)) {
    const char* typestring = pytype_string($input);
    PyErr_Format(PyExc_TypeError,
                 "NumPy array expected.  '%s' given.",
                 typestring);
    SWIG_fail;
  }
  tempArray = (PyArrayObject*) $input;
  /* Check input array dimension and data type */
  int numdim = array_numdims(tempArray);
  if (numdim != 1) {
    PyErr_Format(PyExc_ValueError,
                 "1D array expected.  Array with %d dimensions given.",
                 numdim);
    SWIG_fail;
  }
  if (!PyArray_EquivTypenums(array_type(tempArray), DATA_TYPECODE)) {
    const char* desired_type = typecode_string(DATA_TYPECODE);
    const char* actual_type = typecode_string(array_type($input));
    PyErr_Format(PyExc_TypeError,
                 "Array with data type '%s' expected.  '%s' array given.",
                 desired_type, actual_type);
    SWIG_fail;
  }

  dim1_temp = (DIM_TYPE) array_size(tempArray, 0);
  data_temp = (DATA_TYPE *) array_data(tempArray);

  $1 = &dim1_temp;
  $2 = &data_temp;
}
%typemap(argout)
  (DIM_TYPE* DIM1, DATA_TYPE** ARGINOUTVIEW_ARRAY1)
{
  /* Update the input NumPy array.
   *
   * We only ever expect to have either the same array out again or a
   * null pointer. If the output is null then the input
   * array needs to know it can no longer access its data.
   * We can't alter the numpy array's pointer to the data
   * but by setting the size to zero it will not try to access the data.
   */
  PyArray_DIMS(tempArray$argnum)[0] = dim1_temp$argnum;
}

%enddef

/* Concrete instances of the %numpy_extra_typemaps() macro
 */
%numpy_extra_typemaps(unsigned char , NPY_UBYTE , int)
%numpy_extra_typemaps(short , NPY_SHORT , int)
%numpy_extra_typemaps(unsigned short , NPY_USHORT , int)
%numpy_extra_typemaps(int , NPY_INT , int)
%numpy_extra_typemaps(unsigned int , NPY_UINT , int)
%numpy_extra_typemaps(long , NPY_LONG , int)
%numpy_extra_typemaps(unsigned long , NPY_ULONG , int)
%numpy_extra_typemaps(long long , NPY_LONGLONG , int)
%numpy_extra_typemaps(unsigned long long, NPY_ULONGLONG, int)
%numpy_extra_typemaps(float , NPY_FLOAT , int)
%numpy_extra_typemaps(double , NPY_DOUBLE , int)

#endif /* SWIGPYTHON */
