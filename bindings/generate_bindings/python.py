import os
import re

from parse import *
from c import subroutine_c_names

MODULE_DOCSTRING = ("""OpenCMISS (Open Continuum Mechanics, Imaging, """
"""Signal processing and System identification)

A mathematical modelling environment that enables the application of finite
element analysis techniques to a variety of complex bioengineering problems.

This Python module wraps the underlying OpenCMISS Fortran library.

http://www.opencmiss.org
""")

INITIALISE = """WorldCoordinateSystem = CoordinateSystem()
WorldRegion = Region()
Initialise(WorldCoordinateSystem, WorldRegion)
# Don't output errors, we'll include trace in exception
ErrorHandlingModeSet(ErrorHandlingModes.ReturnErrorCode)

# Ignore SIGPIPE generated when closing the help pager when it isn't fully
# buffered, otherwise it gets caught by OpenCMISS and crashes the interpreter
signal.signal(signal.SIGPIPE, signal.SIG_IGN)
"""

PREFIX = 'CMISS'


def generate(cm_path, args):
    """Generate the OpenCMISS Python module

    This wraps the lower level extension module created by SWIG
    """

    module = open(os.sep.join((cm_path, 'bindings', 'python', 'opencmiss',
        'CMISS.py')), 'w')

    library = LibrarySource(cm_path)

    module.write('"""%s"""\n\n' % MODULE_DOCSTRING)
    module.write("import _opencmiss_swig\n")
    module.write("import signal\n")
    module.write("from _utils import (CMISSError, CMISSType, Enum,\n"
        "    wrap_cmiss_routine as _wrap_routine)\n\n")

    types = sorted(library.lib_source.types.values(), key=attrgetter('name'))
    for type in types:
        module.write(type_to_py(type))
        module.write('\n' * 3)

    for routine in library.unbound_routines:
        module.write(routine_to_py(routine))
        module.write('\n' * 3)

    (enums, ungrouped_constants) = library.group_constants()
    for e in enums:
        module.write(enum_to_py(e))
        module.write('\n' * 3)
    if ungrouped_constants:
        for c in ungrouped_constants:
            doxygen_comment = remove_doxygen_commands(c.doxygen_comment)
            if doxygen_comment.strip():
                module.write("%s = %d  # %s\n" % (c.name[5:], c.value,
                    doxygen_comment))
            else:
                module.write("%s = %d\n" % (c.name[5:], c.value))
        module.write('\n')

    module.write(INITIALISE)
    module.close()


def type_to_py(type):
    """Convert CMISS type to Python class"""

    cmiss_type = type.name[len(PREFIX):-len('Type')]
    docstring = remove_doxygen_commands('\n    '.join(type.comment_lines))

    py_class = ["class %s(CMISSType):" % cmiss_type]
    py_class.append('    """%s\n    """\n' % docstring)
    py_class.append("    def __init__(self):")
    py_class.append('        """Initialise a null %s"""\n' % type.name)
    py_class.append("        self.cmiss_type = "
        "_wrap_routine(_opencmiss_swig.%sInitialise, None)\n" % type.name)

    for method in type.methods:
        if not method.name.endswith('TypeInitialise'):
            py_class.append(py_method(type, method))
            py_class.append('')

    return '\n'.join(py_class).rstrip()


def py_method(type, routine):
    """Write subroutine as method of Python class"""

    c_name = subroutine_c_names(routine)[0]
    name = c_name[len(type.name) - len('Type'):]
    if name == 'TypeFinalise':
        name = 'Finalise'
    create_start_name = type.name[:-len('Type')] + 'CreateStart'

    if c_name.startswith(create_start_name):
        parameters = routine.parameters[:-1]
    else:
        parameters = routine.parameters[1:]

    parameters = add_size_parameters(parameters)

    py_args = [p.name for p in parameters if p.intent != 'OUT']
    method_args = ', '.join(['self'] + py_args)
    if c_name.startswith(create_start_name):
        py_swig_args = ', '.join(py_args + ['self'])
    else:
        py_swig_args = ', '.join(['self'] + py_args)

    docstring = [remove_doxygen_commands(
            '\n        '.join(routine.comment_lines))]
    docstring.append('\n\n')
    docstring.append(' ' * 8)
    docstring.append('\n        '.join(
        parameters_docstring(parameters).splitlines()))
    docstring = ''.join(docstring).strip()

    method = ["    def %s(%s):" % (name, method_args)]
    method.append('        """%s\n        """\n' % docstring)
    method.append('        return _wrap_routine(_opencmiss_swig.%s, [%s])' %
        (c_name, py_swig_args))

    return '\n'.join(method)


def routine_to_py(routine):
    c_name = subroutine_c_names(routine)[0]
    name = c_name[len(PREFIX):]

    parameters = routine.parameters[:]
    parameters = add_size_parameters(parameters)

    docstring = [remove_doxygen_commands('\n    '.join(routine.comment_lines))]
    docstring.append('')
    docstring.append(' ' * 4 + '\n    '.join(
            parameters_docstring(parameters).splitlines()))
    docstring = '\n'.join(docstring).strip()

    args = ', '.join([p.name for p in parameters if p.intent != 'OUT'])

    py_routine = ["def %s(%s):" % (name, args)]
    py_routine.append('    """%s\n    """\n' % docstring)
    py_routine.append('    return _wrap_routine(_opencmiss_swig.%s, [%s])' %
                      (c_name, args))

    return '\n'.join(py_routine)


def parameters_docstring(parameters):
    """Create docstring section for parameters and return values"""

    return_values = []
    docstring = []
    for param in parameters:
        if param.intent == 'OUT':
            return_values.append(param)
        else:
            docstring.append(':param %s: %s' %
                (param.name, replace_doxygen_commands(param)))
            docstring.append(':type %s: %s' %
                (param.name, param_type_comment(param)))
    return_comments = [return_comment(r) for r in return_values]
    if len(return_values) == 0:
        docstring.append(':rtype: None')
    elif len(return_values) == 1:
        docstring.append(':returns: %s' % (return_comments[0]))
        docstring.append(':rtype: %s' % (param_type_comment(return_values[0])))
    else:
        docstring.append(':returns: (%s)' %
            (', '.join([c.rstrip('.')for c in return_comments])))
        docstring.append(':rtype: tuple')

    return '\n'.join(docstring)


def return_comment(return_param):
    """Fix comment describing return value"""

    comment = replace_doxygen_commands(return_param)

    on_return = 'on return, '
    if comment.lower().startswith(on_return):
        comment = (comment[len(on_return)].upper() +
            comment[len(on_return) + 1:])
    if not comment.strip():
        return 'No description'
    return comment.strip()


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
        type = param.type_name[len(PREFIX):-len('Type')]
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


def remove_doxygen_commands(comment):
    see_re = r'\.?\s*\\see\s*[^\s]*'
    match = re.search(see_re, comment)
    if match:
        comment = comment[0:match.start(0)] + comment[match.end(0):]
    return comment.strip()


def replace_doxygen_commands(param):
    """Replace doxygen see command with a reference to the Python enum class"""

    comment = param.doxygen

    if param.var_type == Parameter.INTEGER:
        see_re = r'\.?\s*\\see\s*OPENCMISS_([^\s,\.]*)'
        match = re.search(see_re, comment)
        if match:
            enum = match.group(1)
            if enum is not None:
                if enum.startswith(PREFIX):
                    enum = enum[len(PREFIX):]
                comment = comment[0:match.start(0)]
                if param.intent == 'IN':
                    comment += '. Must be a value from the ' + enum + ' enum.'
                else:
                    comment += '. Will be a value from the ' + enum + ' enum.'

    return comment


def enum_to_py(enum):
    """Create a Python class to represent and enum"""

    output = []
    if enum.name.startswith(PREFIX):
        name = enum.name[len(PREFIX):]
    else:
        name = enum.name
    output.append("class %s(Enum):" % name)
    output.append('    """%s\n    """\n' % enum.comment)
    constant_names = remove_prefix_and_suffix(
            [c.name for c in enum.constants])
    for (constant, constant_name) in zip(enum.constants, constant_names):
        doxygen_comment = remove_doxygen_commands(constant.doxygen_comment)
        if doxygen_comment.strip():
            output.append("    %s = %d  # %s" %
                    (constant_name, constant.value, doxygen_comment))
        else:
            output.append("    %s = %d" % (constant_name, constant.value))
    return '\n'.join(output)


def remove_prefix_and_suffix(names):
    """Remove any common prefix and suffix from a list
    of enum names. These are redundant due to the enum
    class name"""

    if len(names) == 0:
        return names

    prefix_length = 0
    suffix_length = 0
    if len(names) == 1:
        # Special cases we have to specify
        if names[0] == 'CMISSControlLoopNode':
            prefix_length = len('CMISSControlLoop')
        elif names[0] == 'CMISSEquationsSetHelmholtzEquationTwoDim1':
            prefix_length = len('CMISSEquationsSetHelmholtzEquation')
        elif names[0] == 'CMISSEquationsSetPoiseuilleTwoDim1':
            prefix_length = len('CMISSEquationsSetPoiseuille')
        elif names[0] == 'CMISSEquationsSetFiniteElasticityCylinder':
            prefix_length = len('CMISSEquationsSetFiniteElasticity')
        else:
            sys.stderr.write("Warning: Found an unknown enum "
                    "group with only one name: %s.\n" % names[0])
    else:
        min_length = min([len(n) for n in names])

        for i in range(min_length):
            chars = [n[i] for n in names]
            if chars.count(chars[0]) == len(chars):
                prefix_length += 1
            else:
                break

        for i in range(min_length):
            chars = [n[-i - 1] for n in names]
            if chars.count(chars[0]) == len(chars):
                suffix_length += 1
            else:
                break

        # Make sure the suffix starts with uppercase.  So we get eg.
        # EquationsLumpingTypes.Unlumped and Lumped rather than Unl and L
        # Do the same for the prefix so that TwoDim and ThreeDim don't become
        # woDim and hreeDim.  This breaks with a CMISS or CMISSCellML prefix
        # for example though.
        #
        # Constants will change to capitals with underscores soon so this
        # won't be an issue then, we can just check the prefix ends with an
        # underscore
        if prefix_length > 0:
            prefix = names[0][0:prefix_length]
            if prefix == PREFIX:
                pass
            elif prefix == PREFIX + 'CellML':
                pass
            else:
                while names[0][prefix_length - 1].isupper():
                    prefix_length -= 1
        if suffix_length > 0:
            while names[0][-suffix_length].islower():
                suffix_length -= 1

    if suffix_length == 0:
        new_names = [name[prefix_length:] for name in names]
    else:
        new_names = [name[prefix_length:-suffix_length] for name in names]
    for (i, name) in enumerate(new_names):
        # Eg. NoOutputType should become None, not No
        if name == 'No':
            new_names[i] = 'NONE'
        elif name == 'None':
            # Can't assign to None
            new_names[i] = 'NONE'
        elif name[0].isdigit():
            new_names[i] = digit_to_word(name[0]) + name[1:]
        elif name.endswith('VariableType'):
            # The NumberOfVariableSubtypes in this enum stuffs everything up
            new_names[i] = name[:-len('VariableType')]

    return new_names


def digit_to_word(digit):
    words = {
        0: 'Zero',
        1: 'One',
        2: 'Two',
        3: 'Three',
        4: 'Four',
        5: 'Five',
        6: 'Six',
        7: 'Seven',
        8: 'Eight',
        9: 'Nine'
    }
    return words[int(digit)]


class SizeParameter(object):
    def __init__(self, name, doxygen):
        self.name = name
        self.doxygen = doxygen
        self.var_type = Parameter.INTEGER
        self.intent = 'IN'
        self.array_dims = 0


def add_size_parameters(parameters):
    """Returns a new list of parameters, inserting extra size parameters
    required when retrieving an array"""

    new_parameters = []

    for param in parameters:
        if param.intent == 'OUT':
            if param.array_dims == 0:
                pass
            elif param.array_dims == 1:
                new_parameters.append(SizeParameter(param.name + 'Size',
                    'Expected length of %s list' % param.name))
            else:
                new_parameters.extend([SizeParameter(param.name + 'Size%d' % i,
                        'Expected length of dimension %d for %s list' %
                        (i, param.name))
                        for i in range(1, param.array_dims + 1)])

        new_parameters.append(param)

    return new_parameters
