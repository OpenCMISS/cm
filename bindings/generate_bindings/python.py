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
        "    wrap_cmiss_routine as _wrap_routine)\n\n\n")

    types = sorted(library.lib_source.types.values(), key=attrgetter('name'))
    for type in types:
        module.write(type_to_py(type))
        module.write('\n' * 3)

    for routine in library.unbound_routines:
        try:
            module.write(routine_to_py(routine))
            module.write('\n' * 3)
        except UnsupportedParameterError:
            sys.stderr.write("Skipping routine with unsupported "
                "parameter: %s\n" % routine.name)

    (enums, ungrouped_constants) = library.group_constants()
    for e in enums:
        module.write(enum_to_py(e))
        module.write('\n' * 3)
    if ungrouped_constants:
        for c in ungrouped_constants:
            doxygen_comment = remove_doxygen_commands(c.comment)
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

    # Find initialise routine
    for method in type.methods:
        if method.name.endswith('_Initialise'):
            initialise_method = method.name
            break
        if method.name.endswith('TypeInitialise'):
            initialise_method = method.name
            break
    else:
        raise RuntimeError("Couldn't find initialise routine for %s" %
                type.name)

    py_class = ["class %s(CMISSType):" % cmiss_type]
    py_class.append('    """%s\n    """\n' % docstring)
    py_class.append("    def __init__(self):")
    py_class.append('        """Initialise a null %s"""\n' % type.name)
    py_class.append("        self.cmiss_type = "
        "_wrap_routine(_opencmiss_swig.%s, None)\n" % initialise_method)

    for method in type.methods:
        if (not method.name.endswith('TypeInitialise') and
                not method.name.endswith('_Initialise')):
            try:
                py_class.append(py_method(type, method))
                py_class.append('')
            except UnsupportedParameterError:
                sys.stderr.write("Skipping routine with unsupported "
                    "parameter: %s\n" % method.name)

    for (name, get_method, set_method, docstring) in type_properties(type):
        py_class.append('    %s = property(%s, %s, None, """%s""")\n' %
            (lower_camel(name), get_method, set_method, docstring))

    return '\n'.join(py_class).rstrip()


def type_properties(type):
    """Returns a list of tuples representing properties of the type

    Each tuple has the property name, the get method name, the set method
    name, and the property docstring. There may be properties with only a
    get method or only a set method.
    """

    get_methods = dict(
        (method_name(type, m)[:-3], m)
        for m in type.methods
        if method_name(type, m).endswith("Get") and len(m.parameters) == 2)

    set_methods = dict(
        (method_name(type, m)[:-3], m)
        for m in type.methods
        if method_name(type, m).endswith("Set") and len(m.parameters) == 2)

    # Use comment from get method rather than set method if both exist
    all_properties = set_methods.copy()
    all_properties.update(get_methods)
    return [
        (p,
        p + 'Get' if p in get_methods else 'None',
        p + 'Set' if p in set_methods else 'None',
        property_docstring(all_properties[p]))
        for p in all_properties]


def property_docstring(property):
    """Return a docstring for a property, given the comment describing
    either the get or set routine."""

    comment_lines = property.comment_lines

    comment = '\n'.join(comment_lines)
    start_re = re.compile(
        r'^((sets\/changes)|(set)|(get)|(return))s?\s+', re.IGNORECASE)
    comment = start_re.sub('', comment, 1)
    try:
        return comment[0].upper() + comment[1:]
    except IndexError:
        return comment


def method_name(type, routine):
    "Return the name of a method of an object"""

    c_name = subroutine_c_names(routine)[0]
    if '_' in c_name:
        name = c_name.split('_')[-1]
    elif (c_name.startswith('CMISSFieldML') and
            not c_name.startswith('CMISSFieldMLIO')):
        # Special case for FieldML routines that start
        # with FieldML but take a CMISSFieldMLIOType, although
        # some start with CMISSFieldMLIO...
        name = c_name[len('CMISSFieldML'):]
    else:
        # Old code style
        name = c_name[len(type.name) - len('Type'):]
    if name == 'TypeFinalise':
        name = 'Finalise'
    return name


def py_method(type, routine):
    """Write subroutine as method of Python class"""

    if type.name == 'CMISSFieldMLIOType':
        type_name = 'CMISSFieldMLType'
    else:
        type_name = type.name

    name = method_name(type, routine)
    c_name = subroutine_c_names(routine)[0]
    create_start_name = type_name[:-len('Type')] + 'Create'

    if c_name.startswith('CMISSFieldMLOutputCreate'):
        # Have to account for this specifically as it doesn't
        # match other create start routines
        is_create_start = True
        self_parameter = routine.parameters[-1]
        all_parameters = routine.parameters[0:-1]
    elif (c_name.startswith(create_start_name) and
            routine.parameters[-1].var_type == Parameter.CUSTOM_TYPE and
            routine.parameters[-1].type_name  == type.name):
        is_create_start = True
        # Last parameter is self parameter
        self_parameter = routine.parameters[-1]
        all_parameters = routine.parameters[0:-1]
    else:
        is_create_start = False
        self_parameter = routine.parameters[0]
        all_parameters = routine.parameters[1:]

    (pre_code, py_args, swig_args) = process_parameters(all_parameters)

    if is_create_start:
        swig_args = swig_args + [self_parameter.name]
    else:
        swig_args = [self_parameter.name] + swig_args
    py_args = (['self'] + py_args)

    docstring = [remove_doxygen_commands(
            '\n        '.join(routine.comment_lines))]
    docstring.append('\n\n')
    docstring.append(' ' * 8)
    docstring.append("\n        ".join(
        parameters_docstring(all_parameters).splitlines()))
    docstring = ''.join(docstring).strip()

    method = ["    def %s(%s):" % (name, ', '.join(py_args))]
    method.append('        """%s\n        """\n' % docstring)
    method.append("        %s = self" % self_parameter.name)
    for line in pre_code:
        method.append("        %s" % line)
    method.append("        return _wrap_routine(_opencmiss_swig.%s, [%s])" %
        (c_name, ', '.join(swig_args)))

    return '\n'.join(method)


def routine_to_py(routine):
    c_name = subroutine_c_names(routine)[0]
    name = c_name[len(PREFIX):]

    docstring = [remove_doxygen_commands('\n    '.join(routine.comment_lines))]
    docstring.append('')
    docstring.append(' ' * 4 + '\n    '.join(
            parameters_docstring(routine.parameters).splitlines()))
    docstring = '\n'.join(docstring).strip()

    (pre_code, py_args, swig_args) = process_parameters(routine.parameters)

    py_routine = ["def %s(%s):" % (name, ', '.join(py_args))]
    py_routine.append('    """%s\n    """\n' % docstring)
    for line in pre_code:
        py_routine.append("    %s" % line)
    py_routine.append('    return _wrap_routine(_opencmiss_swig.%s, [%s])' %
                      (c_name, ', '.join(swig_args)))

    return '\n'.join(py_routine)


def check_parameter(parameter):
    """Checks whether a parameter is supported by the Python bindings

    Raises an UnsupportedParameterError if it isn't supported,
    otherwise does nothing.
    """

    if parameter.array_dims > 0 and parameter.pointer:
        raise UnsupportedParameterError("Python bindings don't support "
            "passing arrays by pointer")


def parameters_docstring(parameters):
    """Create docstring section for parameters and return values"""

    return_values = []
    docstring = []
    for param in parameters:
        if param.intent == 'OUT':
            return_values.append(param)
            if (param.array_dims == param.required_sizes == 1 and
                    param.var_type != Parameter.CHARACTER):
                docstring.append(':param %sSize: %s' %
                    (param.name, "Size of " + param.name + " to allocate."))
            elif (param.array_dims > 1 and
                    param.array_dims == param.required_sizes):
                docstring.append(':param %sSizes: %s' % (param.name,
                    "Tuple of dimensions of %s to allocate, with length %d." %
                    (param.name, param.array_dims)))
        else:
            docstring.append(':param %s: %s' %
                (param.name, replace_doxygen_commands(param)))
            docstring.append(':type %s: %s' %
                (param.name, param_type_comment(param)))
    return_comments = [return_comment(r) for r in return_values]
    if len(return_values) == 0:
        docstring.append(':rtype: None')
    elif len(return_values) == 1:
        docstring.append(':returns: %s. %s' % (return_values[0].name,
            return_comments[0]))
        docstring.append(':rtype: %s' % (param_type_comment(return_values[0])))
    else:
        docstring.append(':returns: (%s)' %
            (', '.join([c.rstrip('.')for c in return_comments])))
        docstring.append(':rtype: tuple. (%s)' % ', '.join(
            param_type_comment(r) for r in return_values))

    return '\n'.join([l.rstrip() for l in docstring])


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
    if param.var_type == Parameter.CHARACTER:
        if param.array_dims == 2:
            type = "Array of %ss" % type
    else:
        if param.array_dims == 1:
            if param.var_type == Parameter.CUSTOM_TYPE:
                type = "Array of %s objects" % type
            else:
                type = "Array of %ss" % type
        elif param.array_dims >= 1:
            if param.var_type == Parameter.CUSTOM_TYPE:
                type = "%dd array of %s objects" % (param.array_dims, type)
            else:
                type = "%dd array of %ss" % (param.array_dims, type)
    return type


def remove_doxygen_commands(comment):
    see_re = r'\.?\s*\\see\s*[^\s]*'
    match = re.search(see_re, comment)
    if match:
        comment = comment[0:match.start(0)] + comment[match.end(0):]
    return comment.strip()


def replace_doxygen_commands(param):
    """Replace doxygen see command with a reference to the Python enum class"""

    comment = param.comment

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
        doxygen_comment = remove_doxygen_commands(constant.comment)
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
            elif prefix == PREFIX + 'FieldML':
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


def process_parameters(parameters):
    """Processes list of parameters

    Adds any extra size parameters and returns parameters used by the
    python module function and parameters sent to the underlying swig
    module, with any extra parameters added.

    Because the Fortran API expects return arrays to be already allocated
    when passed in, the Python SWIG routines need the size of the array
    to allocate to be passed in, which is an additional input parameter.

    For strings, we don't know how long the string will be but these are
    always just labels at the moment, so a maximum output size is defined
    in the SWIG bindings.

    Returns a tuple, the first value is a list of strings of extra code
    called to set parameters.  The second is a list of parameters accepted
    by the python routine.  The third is a list of parameters sent through
    to the underlying SWIG routine.
    """

    pre_code = []
    python_parameters = []
    swig_parameters = []

    for param in parameters:
        check_parameter(param)
        if param.intent in ('IN', 'INOUT'):
            python_parameters.append(param.name)
            swig_parameters.append(param.name)
            if (param.array_dims == 1 and param.required_sizes == 0):
                pre_code.append("assert(len(%s) == %d)" %
                        (param.name, int(param.array_spec[0])))
        elif param.intent == 'OUT' and param.array_dims > 0:
            if (param.array_dims == 1 and
                    param.var_type == Parameter.CHARACTER):
                pass
            elif (param.array_dims == 1 and param.required_sizes == 1):
                python_parameters.append(param.name + 'Size')
                swig_parameters.append(param.name + 'Size')
            elif (param.array_dims > 1 and
                    param.required_sizes == param.array_dims):
                python_parameters.append(param.name + 'Sizes')
                swig_parameters.append(param.name + 'Sizes')
            elif (param.required_sizes == 0):
                pass
            else:
                sys.stderr.write("Warning: Output of array with dimension = "
                    "%d and required sizes = %d is not implemented.")

    return (pre_code, python_parameters, swig_parameters)


def lower_camel(s):
    try:
        return s[0].lower() + s[1:]
    except IndexError:
        return s
