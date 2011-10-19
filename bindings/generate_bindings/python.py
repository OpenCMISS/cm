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
    return_values = []
    for param in routine.parameters:
        if param.intent == 'OUT':
            return_values.append(param)
        else:
            docstring += '    :param %s: %s\n' % (param.name, param.doxygen)
            docstring += '    :type %s: %s\n' % (param.name, param_type(param))
    return_comments = [return_comment(r.doxygen) for r in return_values]
    if len(return_values) == 0:
        docstring += '    :rtype: None\n'
    elif len(return_values) == 1:
        docstring += '    :rtype: %s, %s\n' % (param_type(return_values[0]), return_comments[0])
    else:
        docstring += '    :rtype: tuple (%s)\n' % (', '.join(return_comments))
    docstring = docstring.strip()

    args = ', '.join([p.name for p in routine.parameters if p.intent != 'OUT'])

    py_routine = "def %s(%s):\n" % (name, args)
    py_routine += '    """%s\n    """\n\n' % docstring
    py_routine += '    return _wrap_routine(_opencmiss_swig.%s, (%s))\n' % (c_name, args)
    py_routine += '\n\n'

    return py_routine


def return_comment(comment):
    """Fix comment describing return value
    """
    on_return = 'on return, '
    if comment.lower().startswith(on_return):
        comment = comment[len(on_return):]
    return comment

PARAMETER_TYPES = {
    Parameter.INTEGER: 'int',
    Parameter.FLOAT: 'float',
    Parameter.DOUBLE: 'float',
    Parameter.CHARACTER: 'string',
    Parameter.LOGICAL: 'bool',
    Parameter.CUSTOM_TYPE: None
}

def param_type(param):
    """Python type corresponding to Fortran type"""
    if param.var_type == Parameter.CUSTOM_TYPE:
        type = param.type_name[len('CMISS'):-len('Type')]
    else:
        type = PARAMETER_TYPES[param.var_type]
    if param.array_dims == 1:
        if param.var_type == Parameter.CUSTOM_TYPE:
            type = "Array of %s objects" % type
        else:
            type = "Array of %ss" % type
    elif param.array_dims >= 1:
        if param.var_type == Parameter.CUSTOM_TYPE:
            type = "%dd list of %s objects" % (param.array_dims, type)
        else:
            type = "%dd list of %ss" % (param.array_dims, type)
    return type
