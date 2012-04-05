import unittest

from parse import *


class ParseTestClass(unittest.TestCase):
    def setUp(self):
        # Load source code for test library
        self.library = LibrarySource('tests/example_library')
        self.routines = dict(
            (r.name, r) for r in self.library.public_subroutines)
        self.test_routine = (
            self.routines["CMISSExample_SomeInterfaceObj"])

    def test_string_routines(self):
        """Only routines that use C strings should be included

        This is based purely on the routine name
        """

        self.assertFalse("CMISSStringRoutineVSObj" in self.routines)
        self.assertFalse("CMISSStringRoutineVSNumber" in self.routines)
        self.assertTrue("CMISSStringRoutineCObj" in self.routines)
        self.assertTrue("CMISSStringRoutineCNumber" in self.routines)

    def test_public(self):
        """Check we don't get any non-public routines"""

        self.assertFalse("CMISSNonPublicRoutine" in self.routines)

    def test_array_routines(self):
        """Routine that takes an array rather than scalar should be used

        This is based purely on the routine name
        """

        self.assertFalse("CMISSArrayRoutine0" in self.routines)
        self.assertTrue("CMISSArrayRoutine1" in self.routines)

    def test_parameter_intent(self):
        """Check parameter intents are correct"""

        for param in self.test_routine.parameters:
            if param.name.startswith("Input"):
                assert param.intent == 'IN'
            elif param.name.endswith("Output"):
                assert param.intent == 'OUT'

    def test_parameter_arrays(self):
        """Check array dimensions are correct"""

        for param in self.test_routine.parameters:
            if param.name.endswith("Array2D"):
                self.assertEqual(param.array_dims, 2)
                self.assertEqual(param.array_spec, [':', ':'])
                self.assertEqual(param.required_sizes, 2)
            elif param.name.endswith("Array"):
                self.assertEqual(param.array_dims, 1)
                self.assertEqual(param.array_spec, [':'])
                self.assertEqual(param.required_sizes, 1)
            elif param.var_type == Parameter.CHARACTER:
                self.assertEqual(param.array_dims, 1)
                self.assertEqual(param.array_spec, [':'])
                self.assertEqual(param.required_sizes, 1)
            elif param.name == "ArrayWithSize":
                self.assertEqual(param.array_dims, 1)
                self.assertEqual(param.array_spec, ['2'])
                self.assertEqual(param.required_sizes, 0)
            else:
                self.assertEqual(param.array_dims, 0)
                self.assertEqual(param.array_spec, [])

    def test_parameter_comments(self):
        """Check doxygen comments are correct"""

        for param in self.test_routine.parameters:
            self.assertEqual(param.comment,
                "Comment for %s" % param.name)

    def test_parameters(self):
        """Check number of parameters and that Err isn't included"""

        self.assertEqual(len(self.test_routine.parameters), 10)
        self.assertFalse(
            "Err" in [p.name for p in self.test_routine.parameters])

    def test_enum(self):
        """Check enum has correct parameters"""

        (enums, ungrouped_constants) = self.library.group_constants()
        self.assertEqual(len(ungrouped_constants), 1)
        self.assertEqual(len(enums), 1)
        enum = enums[0]
        self.assertEqual(len(enum.constants), 3)
        self.assertEqual(enum.comment, "Example of an enum")
        self.assertEqual(enum.name, "ExampleEnum")
        self.assertEqual(enum.constants[0].value, 1)
        self.assertFalse("NON_PUBLIC_CONSTANT" in ungrouped_constants)

    def test_method(self):
        """Test that correct routines are set as methods of type"""

        type = self.library.lib_source.types["CMISSExampleType"]
        # Make sure we get Initialise and CreateStart
        self.assertEqual(len(type.methods), 3)


if __name__ == '__main__':
    unittest.main()
