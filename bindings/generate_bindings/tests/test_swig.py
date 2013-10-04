import unittest

import parse
from swig import *
from tests import mocks as m


def typemap_apply(typemap, to):
    return "%%apply (%s){(%s)};" % (typemap, to)


class SWIGTestClass(unittest.TestCase):
    """Test correct typemaps are applied to parameters"""

    def setUp(self):
        pass

    def test_scalar_parameters(self):
        """Test scalar parameters have correct SWIG typemaps"""

        # Input integer
        result = parameter_swig_lines(m.input_integer)
        expected = ""
        self.assertEqual(result[0], expected)

        # Output integer
        result = parameter_swig_lines(m.output_integer)
        expected = typemap_apply(
                "int *DummyOutputScalar",
                "int *test")
        self.assertEqual(result[0], expected)

        # Input real
        result = parameter_swig_lines(m.input_real)
        expected = ""
        self.assertEqual(result[0], expected)

        # Output real
        result = parameter_swig_lines(m.output_real)
        expected = typemap_apply(
                "double *DummyOutputScalar",
                "double *test")
        self.assertEqual(result[0], expected)

    def test_array_parameters(self):
        """Test array parameters have correct SWIG typemaps"""

        # Input 1D array
        result = parameter_swig_lines(m.input_array)
        expected = typemap_apply(
                "const int DIM1, int *IN_ARRAY1",
                "const int testSize, const int *test")
        self.assertEqual(result[0], expected)

        # Input 2D array
        result = parameter_swig_lines(m.input_array_2d)
        expected = typemap_apply(
                "const int DIM1, const int DIM2, "
                    "int *IN_FARRAY2",
                "const int testSize1, const int testSize2, const int *test")
        self.assertEqual(result[0], expected)

        # Input array of known size
        result = parameter_swig_lines(m.output_array_known_size)
        expected = typemap_apply(
                "double ARGOUT_ARRAY1[ANY]",
                "double test[2]")
        self.assertEqual(result[0], expected)

        # Output array
        result = parameter_swig_lines(m.output_array)
        expected = typemap_apply(
                "const int DIM1, int *ARGOUT_ARRAY1",
                "const int testSize, int *test")
        self.assertEqual(result[0], expected)

        # Output array pointer
        result = parameter_swig_lines(m.output_array_pointer)
        expected = typemap_apply(
                "int *DIM1, int **ARGOUTVIEW_ARRAY1",
                "int *testSize, int **test")
        self.assertEqual(result[0], expected)

    def test_string_parameters(self):
        """Test string parameters have correct SWIG typemaps"""

        # Input string
        result = parameter_swig_lines(m.input_string)
        expected = typemap_apply(
                "const int Size, const char *DummyInputString",
                "const int testSize, const char *test")
        self.assertEqual(result[0], expected)

        # Output string
        result = parameter_swig_lines(m.output_string)
        expected = typemap_apply(
                "const int Size, char *DummyOutputString",
                "const int testSize, char *test")
        self.assertEqual(result[0], expected)

        # Input array of strings
        result = parameter_swig_lines(m.input_string_array)
        expected = typemap_apply(
                "const int NumStrings, const int StringLength, "
                    "const char *DummyStringList",
                "const int testNumStrings, const int testStringLength, "
                    "const char *test")
        self.assertEqual(result[0], expected)

    def test_cmiss_parameters(self):
        """Test CMFE type parameters have correct SWIG typemaps"""

        # Input CMFE type
        result = parameter_swig_lines(m.input_cmiss_type)
        expected = ""
        self.assertEqual(result[0], expected)

        # Input array of CMFE types
        result = parameter_swig_lines(m.input_cmiss_type_array)
        expected = typemap_apply(
                "const int ArraySize, const cmfe_DummyType *DummyTypes",
                "const int testSize, const cmfe_TestType *test")
        self.assertEqual(result[0], expected)

if __name__ == '__main__':
    unittest.main()
