import numpy


try:
    print(numpy.get_include())
except AttributeError:
    print(numpy.get_numpy_include())
