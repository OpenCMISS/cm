from unittest import TestLoader, TextTestRunner, TestSuite
from tests import test_parse, test_c, test_swig, test_python


if __name__ == "__main__":
    loader = TestLoader()
    suite = TestSuite((
        loader.loadTestsFromTestCase(test_parse.ParseTestClass),
        loader.loadTestsFromTestCase(test_c.CTestClass),
        loader.loadTestsFromTestCase(test_swig.SWIGTestClass),
        loader.loadTestsFromTestCase(test_python.PythonTestClass)))

    runner = TextTestRunner(verbosity=2)
    runner.run(suite)
