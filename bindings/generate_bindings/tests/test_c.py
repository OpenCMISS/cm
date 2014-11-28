import unittest

import parse
from c import *
from tests import mocks as m


class CTestClass(unittest.TestCase):
    def setUp(self):
        pass

    def test_routine_names(self):
        """Test how routines are renamed for the C bindings

        Only interface routines that take C strings and arrays rather
        than scalars are kept, so those indicators not needed in the names.
        """

        test_obj_names = [
            "CMISSRoutineObj",
            "CMISSRoutine1",
            "CMISSRoutineC1",
            "CMISSRoutineCObj"]
        expected_obj_names = ("CMISSRoutine", "CMISSRoutineC")

        test_num_names = [
            "CMISSRoutineNumber",
            "CMISSRoutineCNumber",
            "CMISSRoutineNumber0",
            "CMISSRoutineCNumber1",
            "CMISSRoutineCCNumber1"]
        expected_num_names = ("CMISSRoutineNum", "CMISSRoutineCNum")

        for routine_name in test_obj_names:
            routine = m.Mock(name=routine_name)
            self.assertEqual(
                subroutine_c_names(routine), expected_obj_names)
        for routine_name in test_num_names:
            routine = m.Mock(name=routine_name)
            self.assertEqual(
                subroutine_c_names(routine), expected_num_names)

    def test_scalar_parameters(self):
        """Test how scalar parameters are translated to C"""

        # Input integer
        c_result = parameter_to_c(m.input_integer)
        c_expected = ["const int test"]
        cf90_result = parameter_c_f90_declaration(m.input_integer)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: test"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Output integer
        c_result = parameter_to_c(m.output_integer)
        c_expected = ["int *test"]
        cf90_result = parameter_c_f90_declaration(m.output_integer)
        cf90_expected = [
                "INTEGER(C_INT), INTENT(OUT) :: test"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Input real
        c_result = parameter_to_c(m.input_real)
        c_expected = ["const double test"]
        cf90_result = parameter_c_f90_declaration(m.input_real)
        cf90_expected = [
                "REAL(C_DOUBLE), VALUE, INTENT(IN) :: test"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Output real
        c_result = parameter_to_c(m.output_real)
        c_expected = ["double *test"]
        cf90_result = parameter_c_f90_declaration(m.output_real)
        cf90_expected = [
                "REAL(C_DOUBLE), INTENT(OUT) :: test"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

    def test_array_parameters(self):
        """Test how array parameters are translated to C"""

        # Input 1D array
        c_result = parameter_to_c(m.input_array)
        c_expected = ["const int testSize", "const int *test"]
        cf90_result = parameter_c_f90_declaration(m.input_array)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testSize",
                "TYPE(C_PTR), VALUE, INTENT(IN) :: testPtr"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Input 2D array
        c_result = parameter_to_c(m.input_array_2d)
        c_expected = ["const int testSize1", "const int testSize2",
                "const int *test"]
        cf90_result = parameter_c_f90_declaration(m.input_array_2d)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testSize1",
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testSize2",
                "TYPE(C_PTR), VALUE, INTENT(IN) :: testPtr"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Output array
        c_result = parameter_to_c(m.output_array)
        c_expected = ["const int testSize", "int *test"]
        cf90_result = parameter_c_f90_declaration(m.output_array)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testSize",
                "TYPE(C_PTR), VALUE, INTENT(IN) :: testPtr"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Output array pointer
        c_result = parameter_to_c(m.output_array_pointer)
        c_expected = ["int *testSize", "int **test"]
        cf90_result = parameter_c_f90_declaration(m.output_array_pointer)
        cf90_expected = [
                "INTEGER(C_INT), INTENT(OUT) :: testSize",
                "TYPE(C_PTR), INTENT(OUT) :: testPtr"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

    def test_string_parameters(self):
        """Test how string parameters are translated to C"""

        # Input string
        c_result = parameter_to_c(m.input_string)
        c_expected = ["const int testSize", "const char *test"]
        cf90_result = parameter_c_f90_declaration(m.input_string)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testSize",
                "TYPE(C_PTR), VALUE, INTENT(IN) :: test"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Output string
        c_result = parameter_to_c(m.output_string)
        c_expected = ["const int testSize", "char *test"]
        cf90_result = parameter_c_f90_declaration(m.output_string)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testSize",
                "TYPE(C_PTR), VALUE, INTENT(IN) :: test"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Input array of strings
        c_result = parameter_to_c(m.input_string_array)
        c_expected = [
                "const int testNumStrings", "const int testStringLength",
                "const char *test"]
        cf90_result = parameter_c_f90_declaration(m.input_string_array)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testNumStrings",
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testStringLength",
                "TYPE(C_PTR), VALUE, INTENT(IN) :: test"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

    def test_cmiss_parameters(self):
        """Test how CMISS type parameters are translated to C"""

        # Input CMISS type
        c_result = parameter_to_c(m.input_cmiss_type)
        c_expected = ["const CMISSTestType test"]
        cf90_result = parameter_c_f90_declaration(m.input_cmiss_type)
        cf90_expected = [
                "TYPE(C_PTR), VALUE, INTENT(IN) :: testPtr"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Output CMISS type
        c_result = parameter_to_c(m.output_cmiss_type)
        c_expected = ["CMISSTestType *test"]
        cf90_result = parameter_c_f90_declaration(m.output_cmiss_type)
        cf90_expected = [
                "TYPE(C_PTR), INTENT(OUT) :: testPtr"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

        # Input array of CMISS types
        c_result = parameter_to_c(m.input_cmiss_type_array)
        c_expected = ["const int testSize", "const CMISSTestType *test"]
        cf90_result = parameter_c_f90_declaration(m.input_cmiss_type_array)
        cf90_expected = [
                "INTEGER(C_INT), VALUE, INTENT(IN) :: testSize",
                "TYPE(C_PTR), VALUE, INTENT(IN) :: testPtr"]
        self.assertEqual(c_result, c_expected)
        self.assertEqual(cf90_result, cf90_expected)

if __name__ == '__main__':
    unittest.main()
