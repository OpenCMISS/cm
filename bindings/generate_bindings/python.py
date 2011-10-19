import os
from parse import *
from c import subroutine_c_name

MODULE_DOCSTRING = """OpenCMISS (Open Continuum Mechanics, Imaging, Signal processing and System identification)

A mathematical modelling environment that enables the application of finite
element analysis techniques to a variety of complex bioengineering problems.

This Python module wraps the underlying OpenCMISS Fortran library.

http://www.opencmiss.org
"""

INITIALISE = """WorldCoordinateSystem = CoordinateSystemTypeInitialise()
WorldRegion = RegionTypeInitialise()
Initialise(WorldCoordinateSystem, WorldRegion)
ErrorHandlingModeSet(ReturnErrorCode) #Don't output errors, we'll include trace in exception
"""

def generate(cm_path,args):
    """
    Generate the Python module that wraps the lower level C module created by SWIG
    """
    module = open(os.sep.join((cm_path,'bindings','python','opencmiss','CMISS.py')),'w')

    library = LibrarySource(cm_path)

    module.write('"""%s"""\n\n' % MODULE_DOCSTRING)
    module.write("import _opencmiss_swig\n")
    module.write("from _utils import CMISSError, CMISSType, wrap_cmiss_routine as _wrap_routine\n\n")

    types = sorted(library.lib_source.types.values(), key=attrgetter('name'))
    for type in types:
        module.write(type_to_py(type))

    for routine in library.unbound_routines:
        module.write(routine_to_py(routine))

    module.write('# OpenCMISS constants\n')
    constants = [c for c in library.public_objects.values() if isinstance(c,Constant)]
    constants=sorted(constants,key=attrgetter('name'))
    for c in constants:
        module.write("%s = %d #%s\n" % (c.name[5:], c.value, c.doxygen_comment))
    module.write('\n')

    module.write(INITIALISE)

    module.close()


def type_to_py(type):
    cmiss_type = type.name[len('CMISS'):-len('Type')]
    docstring = '\n    '.join(type.comment_lines)

    py_class = "class %s(CMISSType):\n" % cmiss_type
    py_class += '    """%s\n    """\n' % docstring
    py_class += '\n\n'

    return py_class


def routine_to_py(routine):
    c_name = subroutine_c_name(routine)[0]
    name = c_name[len('CMISS'):]
    routine.get_parameters()

    docstring = '\n    '.join(routine.comment_lines)
    docstring += '\n\n'
    for param in routine.parameters:
        if param.intent == 'OUT':
            docstring += '    :return %s: %s\n' % (param.name, param.doxygen)
        else:
            docstring += '    :param %s: %s\n' % (param.name, param.doxygen)
    docstring = docstring.strip()

    args = ', '.join([p.name for p in routine.parameters if p.intent != 'OUT'])

    py_routine = "def %s(%s):\n" % (name, args)
    py_routine += '    """%s\n    """\n\n' % docstring
    py_routine += '    return _wrap_routine(_opencmiss_swig.%s, (%s))\n' % (c_name, args)
    py_routine += '\n\n'

    return py_routine

