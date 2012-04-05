import unittest

from python import *
from tests import mocks as m


class PythonTestClass(unittest.TestCase):
    def setUp(self):
        pass

    def test_properties(self):
        """Test that get and set methods are translated to properties"""

        type = m.Mock(
                name="CMISSTestType")
        routines = [
                "CMISSTestLabelGet",
                "CMISSTestLabelSet",
                "CMISSTestSetPropertySet",
                "CMISSTestGetPropertyGet",
                "CMISSTestDoSomething"]
        type.methods = [
                m.Mock(name=r, comment_lines=[], parameters=[None, None])
                for r in routines]
        type.methods.append(
                m.Mock(
                    name="CMISSTestTooManyParamsGet",
                    comment_lines=[],
                    parameters=[None, None, None]))
        properties = [p[0] for p in type_properties(type)]

        self.assertTrue("Label" in properties)
        # If there's only a Get or Set, these are still included but
        # as read or write only properties
        self.assertTrue("SetProperty" in properties)
        self.assertTrue("GetProperty" in properties)
        self.assertFalse("DoSomething" in properties)
        self.assertFalse("TooManyParameters" in properties)

    def test_property_docstring(self):
        """Test conversion of get/set routine comment to property docstring"""

        comments = [
                "Sets/changes the test property.",
                "Returns the test property.",
                "Gets the test property."]
        for comment in comments:
            property = m.Mock(comment_lines=[comment])
            result = property_docstring(property)
            self.assertEqual(result, "The test property.")

    def test_parameter_docstrings(self):
        """Test docstring comments for parameters"""

        parameters = [
                m.input_integer,
                m.input_real,
                m.input_array,
                m.input_string]
        result = parameters_docstring(parameters)
        expected = (":param test: Test comment\n"
                ":type test: int\n"
                ":param test: Test comment\n"
                ":type test: float\n"
                ":param test: Test comment\n"
                ":type test: Array of ints\n"
                ":param test: Test comment\n"
                ":type test: string\n"
                ":rtype: None")
        self.assertEqual(result, expected)

        #Test return values
        parameters = [
                m.output_integer]
        result = parameters_docstring(parameters)
        expected = (":returns: test. Test comment\n"
                ":rtype: int")
        self.assertEqual(result, expected)

        parameters = [
                m.output_integer,
                m.output_real]
        result = parameters_docstring(parameters)
        expected = (":returns: (Test comment, Test comment)\n"
                ":rtype: tuple. (int, float)")
        self.assertEqual(result, expected)

    def test_enum(self):
        """Test how an enum is written in Python"""

        names = ["ONE", "TWO", "THREE"]
        values = [1, 2, 3]
        constants = [
                m.Mock(
                    name="TEST_ENUM_%s_VALUE" % n,
                    value=v,
                    comment="Value comment")
                for (n, v) in zip(names, values)]
        enum = m.Mock(
                name="TestEnum",
                comment="Enum description",
                constants=constants)

        result = enum_to_py(enum)
        expected = ('class TestEnum(Enum):\n'
            '    """Enum description\n'
            '    """\n'
            '\n'
            '    ONE = 1  # Value comment\n'
            '    TWO = 2  # Value comment\n'
            '    THREE = 3  # Value comment')
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
