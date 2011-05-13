import sys
try:
    import ez_setup
    ez_setup.use_setuptools()
except ImportError:
    pass

from setuptools import setup

setup(
    name='html output plugin',
    version='0.1',
    description = 'nose html output plugin for OpenCMISS',
    license = 'GNU LGPL',
    py_modules = ['htmlplug'],
    entry_points = {
        'nose.plugins.0.10': [
            'htmlout = htmlplug:HtmlOutput'
            ]
        }

    )

