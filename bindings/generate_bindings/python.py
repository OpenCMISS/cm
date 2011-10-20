import os
from parse import *
from c import subroutine_c_name

MODULE_DOCSTRING = """OpenCMISS (Open Continuum Mechanics, Imaging, Signal processing and System identification)

A mathematical modelling environment that enables the application of finite
element analysis techniques to a variety of complex bioengineering problems.

This Python module wraps the underlying OpenCMISS Fortran library.

http://www.opencmiss.org
"""

INITIALISE = """WorldCoordinateSystem = CoordinateSystem()
WorldRegion = Region()
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
    """Convert CMISS type to Python class"""

    cmiss_type = type.name[len('CMISS'):-len('Type')]
    docstring = '\n    '.join(type.comment_lines)

    py_class = "class %s(CMISSType):\n" % cmiss_type
    py_class += '    """%s\n    """\n' % docstring
    py_class += '\n\n'

    py_class += "    def __init__(self):\n"
    py_class += '        """Initialise a null %s"""\n\n' % type.name
    py_class += "        self.cmiss_type = _wrap_routine(_opencmiss_swig.%sInitialise, None)\n\n" % type.name

    for method in type.methods:
        if not method.name.endswith('TypeInitialise'):
            py_class += py_method(type, method)

    return py_class


def py_method(type, routine):
    """Write subroutine as method of Python class"""

    c_name = subroutine_c_name(routine)[0]
    name = c_name[len(type.name)-len('Type'):]

    if name.endswith('CreateStart'):
        parameters = routine.parameters[:-1]
    else:
        parameters = routine.parameters[1:]

    py_args = [p.name for p in parameters if p.intent != 'OUT']
    method_args = ', '.join(['self']+py_args)
    if name.endswith('CreateStart'):
        py_swig_args = ', '.join(py_args + ['self'])
    else:
        py_swig_args = ', '.join(['self'] + py_args)

    docstring = '\n        '.join(routine.comment_lines)
    docstring += '\n\n'
    docstring += ' '*8 + '\n        '.join(parameters_docstring(parameters).splitlines())
    docstring = docstring.strip()

    method = "    def %s(%s):\n" % (name, method_args)
    method += '        """%s\n        """\n\n' % docstring
    method += '        return _wrap_routine(_opencmiss_swig.%s, [%s])\n' % (c_name, py_swig_args)
    method += '\n'

    return method


def routine_to_py(routine):
    c_name = subroutine_c_name(routine)[0]
    name = c_name[len('CMISS'):]

    docstring = '\n    '.join(routine.comment_lines)
    docstring += '\n\n'
    docstring += ' '*4 +'\n    '.join(parameters_docstring(routine.parameters).splitlines())
    docstring = docstring.strip()

    args = ', '.join([p.name for p in routine.parameters if p.intent != 'OUT'])

    py_routine = "def %s(%s):\n" % (name, args)
    py_routine += '    """%s\n    """\n\n' % docstring
    py_routine += '    return _wrap_routine(_opencmiss_swig.%s, [%s])\n' % (c_name, args)
    py_routine += '\n\n'

    return py_routine


def parameters_docstring(parameters):
    """Create docstring section for parameters and return values"""

    return_values = []
    docstring = ""
    for param in parameters:
        if param.intent == 'OUT':
            return_values.append(param)
        else:
            docstring += ':param %s: %s\n' % (param.name, param.doxygen)
            docstring += ':type %s: %s\n' % (param.name, param_type_comment(param))
    return_comments = [return_comment(r.doxygen) for r in return_values]
    if len(return_values) == 0:
        docstring += ':rtype: None\n'
    elif len(return_values) == 1:
        docstring += ':returns: %s\n' % (return_comments[0])
        docstring += ':rtype: %s\n' % (param_type_comment(return_values[0]))
    else:
        docstring += ':returns: (%s)\n' % (', '.join([c.rstrip('.') for c in return_comments]))
        docstring += ':rtype: tuple\n'

    return docstring


def return_comment(comment):
    """Fix comment describing return value"""

    on_return = 'on return, '
    if comment.lower().startswith(on_return):
        comment = comment[len(on_return)].upper()+comment[len(on_return)+1:]
    if not comment.strip():
        return 'No description'
    return comment

PARAMETER_TYPES = {
    Parameter.INTEGER: 'int',
    Parameter.FLOAT: 'float',
    Parameter.DOUBLE: 'float',
    Parameter.CHARACTER: 'string',
    Parameter.LOGICAL: 'bool',
    Parameter.CUSTOM_TYPE: None
}

def param_type_comment(param):
    """Python type corresponding to Fortran type for use in docstrings"""

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
