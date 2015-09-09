#!/usr/bin/env python

import os
from sys import platform
from setuptools import setup

#Don't compile the extension module here, just include it as an
#additional file after building it from the main OpenCMISS Makefile

requires = ['numpy']
package_data = {'opencmiss.iron': ['_@IRON_PYTHON_MODULE@.so']}

try:
    if platform == 'darwin':
        os.symlink('@IRON_TARGET_FILE@', 'opencmiss/iron/libiron.dylib')
        package_data['opencmiss.iron'].append('libiron.dylib')

    setup(
        name='OpenCMISS-Iron',
        version='0.4.0',
        description=('Python bindings for the OpenCMISS computational '
                'modelling library Iron.'),
        long_description=('Python bindings to OpenCMISS-Iron. '
                'OpenCMISS (Open Continuum Mechanics, Imaging, Signal processing '
                'and System identification) is a mathematical modelling '
                'environment that enables the application of finite element '
                'analysis techniques to a variety of complex '
                'bioengineering problems.'),
        author='Adam Reeve',
        license='Mozilla Tri-license',
        author_email='aree035@aucklanduni.ac.nz',
        url='http://www.opencmiss.org/',
        install_requires=requires,
        packages=['opencmiss', 'opencmiss.iron'],
        package_data=package_data
    )
finally:
    if platform == 'darwin':
        os.unlink('opencmiss/iron/libiron.dylib')
