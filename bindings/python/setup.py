#!/usr/bin/env python

from distutils.core import setup

#Don't compile the extension module here, just include it as an
#additional file after building it from the main OpenCMISS Makefile

setup(
    name='OpenCMISS-Iron',
    version='0.3',
    description=('Python bindings for the OpenCMISS-Iron computational '
            'modelling library.'),
    long_description=('Python bindings to OpenCMISS-Iron. '
            'OpenCMISS (Open Continuum Mechanics, Imaging, Signal processing '
            'and System identification) is a mathematical modelling '
            'environment that enables the application of finite element '
            'analysis techniques to a variety of complex bioengineering '
            'problems. OpenCMISS-Iron is the computational backend component '
            'of OpenCMISS.'),
    author='Adam Reeve',
    license='Mozilla Tri-license',
    author_email='aree035@aucklanduni.ac.nz',
    url='http://www.opencmiss.org/',
    packages=['opencmiss'],
    package_data={'opencmiss': ['_iron_swig.so']}
)
