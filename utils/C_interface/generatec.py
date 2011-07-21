#!/usr/bin/env python
"""
OpenCMISS C Interface Generation
--------------------------------

This python script generates the C interface for the OpenCMISS library from opencmiss.f90.

It finds integer constants from the source code and includes them in opencmiss.h

It generates C functions in opencmiss_c.f90 from subroutines in opencmiss.f90 and puts the function
declaration in opencmiss.h.

For interfaces where one routine takes parameters as scalars and another takes them as arrays,
only the routine that takes arrays is included. This is done by checking the routine names for
"Number0", "Number1", "Number01" etc. so relies on this naming convention to work correctly.

Limitations
-----------

- Doesn't support multi-dimensional arrays of CMISS types or Logicals, but this can be added if required.

- Doesn't account for the difference in storage order of multi-dimensional arrays between C and Fortran, except in CMISSC2FStrings.

"""

from __future__ import with_statement
import sys
import re
from operator import attrgetter

class LibrarySource(object):
    """
    Holds info on all the library source code
    """

    class SourceFile(object):
        """
        Info for an individual source file
        """

        def __init__(self,source_file,params_only=False):
            """
            Initialise SourceFile object

            Arguments:
            source_file -- Path to the source file
            """
            self.file_path = source_file
            self.public = []
            self.interfaces = {}
            self.subroutines = {}
            self.constants = {}
            self.types = {}
            self.parse_file(params_only)

        def parse_file(self,params_only=False):
            """
            Run through file once, getting everything we'll need
            """
            source_lines = _join_lines(open(self.file_path,'r').read()).splitlines(True)
            in_subroutine = False
            in_interface = False
            in_type = False
            re_subroutine_start = re.compile(r'^\s*(RECURSIVE\s+)?SUBROUTINE\s+([A-Z0-9_]+)\(',re.IGNORECASE)
            re_subroutine_end = re.compile(r'^\s*END\s+SUBROUTINE',re.IGNORECASE)
            re_interface_start = re.compile(r'^\s*INTERFACE\s+([A-Z0-9_]+)',re.IGNORECASE)
            re_interface_end = re.compile(r'^\s*END\s+INTERFACE',re.IGNORECASE)
            re_type_start = re.compile(r'^\s*TYPE\s+([A-Z0-9_]+)',re.IGNORECASE)
            re_type_end = re.compile(r'^\s*END\s+TYPE',re.IGNORECASE)
            re_public = re.compile(r'^\s*PUBLIC\s+([A-Z0-9_,\s]+)',re.IGNORECASE)
            re_constant = re.compile(r'^\s*INTEGER\([A-Z0-9\(\),_\s]+::\s*([A-Z0-9_]+)\s*=\s*([A-Z0-9_\-\.]+)',re.IGNORECASE)
            if params_only:
                for (lineno,line) in enumerate(source_lines):
                    #Integer parameter
                    match = re_constant.search(line)
                    if match:
                        name = match.group(1)
                        self.constants[name] = Constant(name,match.group(2),lineno,line)
            else:
                for (lineno,line) in enumerate(source_lines):
                    if line.strip() == '' or line.strip().startswith('!'):
                        continue
                    #If inside a subroutine, interface or type
                    if in_subroutine:
                        subroutine.lines.append(line)
                        if re_subroutine_end.search(line):
                            in_subroutine = False
                            continue
                    if in_interface:
                        interface.lines.append(line)
                        if re_interface_end.search(line):
                            in_interface = False
                            continue
                    if in_type:
                        type.lines.append(line)
                        if re_type_end.search(line):
                            in_type = False
                            continue
                    #Public declaration
                    match = re_public.search(line)
                    if match:
                        for symbol in match.group(1).split(','):
                            self.public.append(symbol.strip())
                        continue
                    #Integer parameter
                    match = re_constant.search(line)
                    if match:
                        name = match.group(1)
                        self.constants[name] = Constant(name,match.group(2),lineno,line)
                        continue
                    #Subroutine
                    match = re_subroutine_start.search(line)
                    if match:
                        name = match.group(2)
                        in_subroutine = True
                        subroutine = Subroutine(name,lineno,line)
                        self.subroutines[name] = subroutine
                        continue
                    #Interface
                    match = re_interface_start.search(line)
                    if match:
                        name = match.group(1)
                        in_interface = True
                        interface = Interface(name,self,lineno,line)
                        self.interfaces[name] = interface
                        continue
                    #Type
                    match = re_type_start.search(line)
                    if match:
                        name = match.group(1)
                        in_type = True
                        type = Type(name,lineno,line)
                        self.types[name] = type
                        continue


    def __init__(self,lib_source,source_files):
        """
        Load library information from source files

        Arguments:
        lib_source -- Path to library source file
        source_files -- List of other source files used by the library
        """
        self.lib_source = self.SourceFile(lib_source)
        self.sources = [self.SourceFile(source,params_only=True) for source in source_files]
        self.resolve_constants()
        self.public_types=[t for t in self.lib_source.types.values() if t.name in self.lib_source.public]
        self.public_types=sorted(self.public_types,key=attrgetter('name'))
        self.public_constants=[const for const in self.lib_source.constants.values() if const.name in self.lib_source.public]
        self.public_constants=sorted(self.public_constants,key=attrgetter('name'))
        self.public_subroutines=[routine for routine in self.lib_source.subroutines.values() if routine.name in self.lib_source.public]
        for interface in self.lib_source.interfaces.values():
            if interface.name in self.lib_source.public:
                self.public_subroutines+=[self.lib_source.subroutines[routine] for routine in interface.get_subroutines()]
        self.public_subroutines=sorted(self.public_subroutines,key=attrgetter('name'))

        #Remove CMISS...TypesCopy routines, as these are only used within the C bindings
        #Also remove CMISSGeneratedMeshSurfaceGet for now as it takes an allocatable array but will be removed soon anyways.
        self.public_subroutines = filter(lambda r: not (r.name.startswith('CMISSGeneratedMeshSurfaceGet') or r.name.endswith('TypesCopy')),self.public_subroutines)

    def resolve_constants(self):
        """
        Go through all public constants and work out their actual values
        """
        for pub in self.lib_source.public:
            if self.lib_source.constants.has_key(pub):
                self.get_constant_value(pub)

    def get_constant_value(self,constant):
        """
        Get the actual value for a constant from the source files

        Arguments:
        constant -- Name of the constant to get the value for
        """
        assignment = self.lib_source.constants[constant].assignment
        exhausted = False
        while (not self.lib_source.constants[constant].resolved) and (not exhausted):
            for (i,source) in enumerate(self.sources):
                if source.constants.has_key(assignment):
                    if source.constants[assignment].resolved:
                        self.lib_source.constants[constant].value=source.constants[assignment].value
                        self.lib_source.constants[constant].resolved=True
                        break
                    else:
                        assignment = source.constants[assignment].assignment
                        break
                if i == len(self.sources):
                    exhausted = True
        if not self.lib_source.constants[constant].resolved:
            sys.stderr.write("Warning: Couldn't resolve constant value: %s\n" % constant)

    def write_c_header(self,output):
        """
        Write opencmiss.h containing constants, typedefs and routine declarations

        Arguments:
        output -- File to write to
        """
        output.write('/*\n * opencmiss.h automatically generated from opencmiss.f90\n * Do not edit this file directly, instead edit opencmiss.f90 or generatec.py\n */\n\n' + \
            '#ifndef OPENCMISS_H\n' + \
            '#define OPENCMISS_H\n' + \
            '\n/*\n * Defines\n */\n\n' + \
            'const int CMISSTrue = 1;\n' + \
            'const int CMISSFalse = 0;\n' + \
            'const int CMISSNoError = 0;\n' + \
            'const int CMISSPointerIsNULL = -1;\n' + \
            'const int CMISSPointerNotNULL = -2;\n' + \
            'const int CMISSCouldNotAllocatePointer = -3;\n' + \
            'const int CMISSErrorConvertingPointer = -4;\n\n' + \
            '\n/*\n * Struct defs\n */\n\n')
        for t in self.public_types:
            output.write('struct %s_;\n' % t.name)
        output.write('\n/*\n * Type defs\n */\n\n')
        output.write('typedef int CMISSError;\n')
        for t in self.public_types:
            output.write('typedef struct %s_ *%s;\n' % (t.name,t.name))
        output.write('\n/*\n * Parameters\n */\n\n')
        for const in self.public_constants:
            output.write(const.to_c())
        output.write('\n/*\n * Routines\n */\n\n')
        for subroutine in self.public_subroutines:
            output.write(subroutine.to_c_declaration())
        output.write('#endif\n')

    def write_c_f90(self,output):
        """
        Write opencmiss_c.f90 containing Fortran routines

        Arguments:
        output -- File to write to
        """
        output.write('!\n! opencmiss_c.f90 automatically generated from opencmiss.f90\n! Do not edit this file directly, instead edit opencmiss.f90 or generatec.py\n!\n' + \
            'MODULE OPENCMISS_C\n\n' + \
            '  USE ISO_C_BINDING\n' + \
            '  USE ISO_VARYING_STRING\n' + \
            '  USE OPENCMISS\n' + \
            '  USE CMISS_FORTRAN_C\n\n' + \
            '  IMPLICIT NONE\n\n' + \
            '  PRIVATE\n\n' + \
            '  INTEGER(C_INT), PARAMETER :: CMISSTrue = 1\n' + \
            '  INTEGER(C_INT), PARAMETER :: CMISSFalse = 0\n' + \
            '  INTEGER(C_INT), PARAMETER :: CMISSNoError = 0\n' + \
            '  INTEGER(C_INT), PARAMETER :: CMISSPointerIsNULL = -1\n' + \
            '  INTEGER(C_INT), PARAMETER :: CMISSPointerNotNULL = -2\n' + \
            '  INTEGER(C_INT), PARAMETER :: CMISSCouldNotAllocatePointer = -3\n' + \
            '  INTEGER(C_INT), PARAMETER :: CMISSErrorConvertingPointer = -4\n\n')
        output.write('\n'.join(('  PUBLIC %s' % subroutine.c_f90_name for subroutine in self.public_subroutines)))
        output.write('\nCONTAINS\n\n')
        for subroutine in self.public_subroutines:
            output.write(subroutine.to_c_f90())
        output.write('END MODULE OPENCMISS_C')


class Constant(object):
    """
    Information on a public constant
    """
    def __init__(self,name,assignment,lineno,line):
        """
        Initialise Constant

        Arguments:
        name -- Variable name
        assignment -- Value or another variable assigned to this variable
        lineno -- Line number in the source file where this variable is defined
        line -- Contents of line defining this Constant
        """
        self.name = name
        self.lineno = lineno
        self.line = line
        self.assignment = assignment
        try:
            self.value = int(self.assignment)
            self.resolved = True
        except ValueError:
            try:
                self.value = float(self.assignment)
                self.resolved = True
            except ValueError:
                self.value = None
                self.resolved = False
    def to_c(self):
        """
        Return the C definition of this constant
        """
        if self.resolved:
            return 'static int %s = %d;\n' % (self.name,self.value)
        else:
            return ''


class Interface(object):
    """
    Information on an interface
    """
    def __init__(self,name,source,lineno,line):
        self.name = name
        self.source = source
        self.lineno = lineno
        self.line = line
        self.lines = [line]

    def get_subroutines(self):
        """
        Find the subroutines for an interface

        Choose the one with the highest number if there are options. This
        corresponds to the routine that takes array parameters

        Returns a list of subroutines
        """
        subroutines = []
        obj_routines = []
        number_routines = []
        other_dim_routines = []
        other_routines = []
        routine_re=re.compile(r'MODULE PROCEDURE ([A-Z0-9_]+)',re.IGNORECASE)
        varying_string_re=re.compile(r'VSC*(Obj|Number|)[0-9]*$',re.IGNORECASE)
        for line in self.lines:
            match = routine_re.search(line)
            if match:
                routine_name = match.group(1)
                obj_match = re.search(r'ObjC?([0-9]*)C?$',routine_name)
                number_match = re.search(r'NumberC?([0-9]*)C?$',routine_name)
                other_dim_match = re.search(r'([0-9]+)C?$',routine_name)
                if varying_string_re.search(routine_name):
                    #Don't include routines using varying_string parameters
                    pass
                elif obj_match:
                    obj_routines.append((routine_name,obj_match.group(1)))
                elif number_match:
                    number_routines.append((routine_name,number_match.group(1)))
                elif other_dim_match:
                    number_routines.append((routine_name,other_dim_match.group(1)))
                else:
                    other_routines.append(routine_name)
        #find routines with highest number (parameters passed as arrays rather than scalars)
        if len(obj_routines) > 0:
            max_routine = obj_routines[0]
            for routine in obj_routines[1:]:
                try:
                    if int(routine[1]) > int(max_routine[1]):
                        max_routine = routine
                except ValueError:
                    subroutines.append(routine[0])
                    self.source.subroutines[routine[0]].interface = self
            subroutines.append(max_routine[0])
            self.source.subroutines[max_routine[0]].interface = self
        if len(number_routines) > 0:
            max_routine = number_routines[0]
            for routine in number_routines[1:]:
                try:
                    if int(routine[1]) > int(max_routine[1]):
                        max_routine = routine
                except ValueError:
                    subroutines.append(routine[0])
                    self.source.subroutines[routine[0]].interface = self
            subroutines.append(max_routine[0])
            self.source.subroutines[max_routine[0]].interface = self
        if len(other_dim_routines) > 0:
            max_routine = other_dim_routines[0]
            for routine in number_routines[1:]:
                if int(routine[1]) > int(max_routine[1]):
                    max_routine = routine
            subroutines.append(max_routine[0])
            self.source.subroutines[max_routine[0]].interface = self
        if len(other_routines) > 0:
            for routine in other_routines:
                subroutines.append(routine)
                self.source.subroutines[routine].interface = self
        return subroutines


class Subroutine(object):
    """
    Store information for a subroutine
    """
    def __init__(self,name,lineno,line):
        self.name = name
        self.lineno = lineno
        self.line = line
        self.lines = [line]
        self.parameters = None
        self.interface = None
        self._set_c_name()

    def get_parameters(self):
        """
        Get details of the subroutine parameters

        Sets the Subroutines parameters property as a list of all parameters
        """
        self.parameters = []
        match = re.search(r'^\s*(RECURSIVE\s+)?SUBROUTINE\s+([A-Z0-9_]+)\(([A-Z0-9_,\s]*)\)',self.line,re.IGNORECASE)
        parameters = [p.strip() for p in match.group(3).split(',')]
        try:
            parameters.remove('Err')
        except ValueError:
            sys.stderr.write("Warning: Routine doesn't take Err parameter: %s\n" % self.name)
        for param in parameters:
            param_pattern = r"""
            ^\s*([A-Z_]+\s*(\([A-Z_=\*0-9]+\))?) # parameter type at start of line, followed by possible extra specification in brackets
            \s*([A-Z0-9\s_\(\):,\s]+)?\s*        # extra specifications such as intent
            ::
            [A-Z_,\s\(\):]*                      # Allow for other parameters to be included on the same line
            [,\s:]                               # Make sure we matched the full parameter name
            %s                                   # Parameter name
            (\([0-9,:]+\))?                      # Array dimensions if present
            [,\s$]                               # Whitespace, comma or end of line to make sure we've matched the full parameter name
            """ % param
            param_re = re.compile(param_pattern,re.IGNORECASE|re.VERBOSE)
            for line in self.lines:
                match = param_re.search(line)
                if match:
                    param_type = match.group(1)
                    type_pt2 = match.group(2)
                    extra_stuff = match.group(3)
                    if match.group(4) is not None:
                        array = match.group(4).replace('(','').replace(')','')
                    else:
                        array = ''
                    self.parameters.append(Parameter(param,self,param_type,type_pt2,extra_stuff,array))
                    break
            if not match:
                raise RuntimeError, "Couldn't find parameter %s for subroutine %s" % (param,self.name)
        # work out parameters passed to c_f90 routine
        self.c_parameters = []
        for param in self.parameters:
            self.c_parameters.extend(param.to_c())

    def _set_c_name(self):
        """
        Get the name of the routine as used from C and as used in opencmiss_c.f90

        Sets the subroutines c_name and c_f90_name parameters
        """
        #Note that some interfaces have a 'C' variant with C in the name already:
        if re.search(r'Obj[0-9]*C?$',self.name):
            #For routines with a Region and Interface variant, the Region one is default
            #and Region is removed from the name
            self.c_name = re.sub(r'(Region)?C?(Obj)[0-9]*C?$','',self.name)
            self.c_f90_name = re.sub(r'(Region)?C?Obj[0-9]*C?$','C',self.name)
        elif re.search(r'Number[0-9]*C?$',self.name):
            self.c_name = re.sub(r'C?Number[0-9]*C?$',r'Num',self.name)
            self.c_f90_name = re.sub(r'C?Number[0-9]*C?$',r'CNum',self.name)
        else:
            self.c_name = re.sub(r'C?[0-9]*$','',self.name)
            self.c_f90_name = re.sub(r'C?[0-9]*$',r'C',self.name)
        #Make sure names are different, might have matched a C at the end of the subroutine name
        if self.c_f90_name == self.name:
            self.c_f90_name = self.name+'C'

    def to_c_declaration(self):
        """
        Returns the function declaration in C
        """
        if self.parameters is None:
            self.get_parameters()
        output = 'CMISSError %s(' % self.c_name
        output += ',\n    '.join(self.c_parameters)
        output = output+');\n\n'
        return output

    def local_c_f90_vars(self):
        """
        Returns a list of extra local variables required for this subroutine, for use in converting
        to and from C variables
        """
        local_variables = [param.local_variables for param in self.parameters]
        return list(_chain_iterable(local_variables))

    def pre_call(self):
        """
        Returns a list of lines to add before calling the Fortran routine for converting
        to and from C
        """
        pre_lines = [param.pre_call for param in self.parameters]
        return list(_chain_iterable(pre_lines))

    def post_call(self):
        """
        Returns a list of lines to add after calling the Fortran routine for converting
        to and from C
        """
        #reverse to keep everything in order, as some if statements need to line up
        #from the pre_call lines
        post_lines = [param.post_call for param in reversed(self.parameters)]
        return list(_chain_iterable(post_lines))

    def to_c_f90(self):
        """
        Returns the C function implemented in Fortran for opencmiss_c.f90
        """
        if self.parameters is None:
            self.get_parameters()
        c_f90_params = [param.c_f90_names() for param in self.parameters]
        c_f90_declarations = [param.c_f90_declaration() for param in self.parameters]
        if self.interface is not None:
            function_call = self.interface.name
        else:
            function_call = self.name

        output =  '  FUNCTION %s(%s) &\n    & BIND(C, NAME="%s")\n\n' %(self.c_f90_name,','.join(c_f90_params),self.c_name)
        output += '    !Argument variables\n'
        output += '    %s\n' % '\n    '.join(c_f90_declarations)
        output += '    !Function return variable\n'
        output += '    INTEGER(C_INT) :: %s\n' % self.c_f90_name
        output += '    !Local variables\n'

        content = []
        content.extend(self.local_c_f90_vars())
        content.append('')
        content.append('%s = CMISSNoError' % self.c_f90_name)
        content.extend(self.pre_call())
        content.append('CALL %s(%s)' % (function_call,','.join([p.call_name() for p in self.parameters]+[self.c_f90_name])))
        content.extend(self.post_call())
        output += _indent_lines(content,2,4)

        output += '\n    RETURN\n\n'
        output += '  END FUNCTION %s\n\n' % self.c_f90_name
        output += '  !\n'
        output += '  !'+'='*129+'\n'
        output += '  !\n\n'
        output = '\n'.join([_fix_length(line) for line in output.split('\n')])
        return output


class Parameter(object):
    """
    Information on a subroutine parameter
    """
    #Parameter types enum:
    (INTEGER, \
    FLOAT, \
    DOUBLE, \
    CHARACTER, \
    LOGICAL, \
    CUSTOM_TYPE) \
        = range(6)

    #Corresponding variable types for C
    CTYPES = ('int','float','double','char','int',None)

    #Variable types as used in opencmiss_c.f90
    F90TYPES = ('INTEGER(C_INT)','REAL(C_FLOAT)','REAL(C_DOUBLE)','CHARACTER(LEN=1,KIND=C_CHAR)','INTEGER(C_INT)','TYPE(C_PTR)')

    def __init__(self,name,routine,param_type,type_pt2,extra_stuff,array=''):
        """
        Initialise a parameter

        Arguments:
        name -- Parameter name
        routine -- Pointer back to the subroutine this parameter belongs to
        param_type -- String from the parameter declaration
        type_pt2 -- Any extra parameter type specification, eg "(DP)" for a real
        extra_stuff -- Any extra parameter properties listed after the type, including intent
        array -- The array dimensions included after the parameter name if they exist, otherwise an empty string
        """
        self.name = name
        self.routine = routine
        self.pointer = False
        intent = None
        if extra_stuff is not None:
            match = re.search(r'INTENT\(([A-Z]+)\)?',extra_stuff,re.IGNORECASE)
            if match is not None:
                intent = match.group(1)
            if extra_stuff.find('DIMENSION') > -1:
                sys.stderr.write("Warning: Ignoring DIMENSION specification on parameter %s of routine %s\n" % (self.name,routine.name))
                sys.stderr.write("         Using DIMENSION goes against the OpenCMISS style guidelines.\n")
            if extra_stuff.find('POINTER') > -1:
                self.pointer = True
        #Get parameter intent
        if intent is None:
            #cintent is the intent used in opencmiss_c.f90, which may be different to the intent in opencmiss.f90
            self.intent = 'INOUT'
            self.cintent = 'INOUT'
            sys.stderr.write("Warning: No intent for parameter %s of routine %s\n" % (self.name,routine.name))
        else:
            self.intent = intent
            self.cintent = intent
        #Get array dimensions and work out how many dimension sizes are variable
        if array != '':
            self.array_spec = [a.strip() for a in array.split(',')]
            self.array_dims = len(self.array_spec)
            self.required_sizes = self.array_spec.count(':')
            if self.array_dims > 0 and not self.pointer:
                #Need to pass C pointer by value
                self.cintent = 'IN'
        else:
            self.array_spec = []
            self.array_dims = 0
            self.required_sizes = 0
        #Work out the type of parameter
        param_type = param_type.upper()
        if param_type.startswith('INTEGER'):
            self.var_type = Parameter.INTEGER
        elif param_type.startswith('REAL'):
            self.precision = type_pt2.replace('(','').replace(')','')
            if self.precision == 'DP':
                self.var_type = Parameter.DOUBLE
            else:
                self.var_type = Parameter.FLOAT
        elif param_type.startswith('CHARACTER'):
            self.var_type = Parameter.CHARACTER
            #Add extra dimension, 1D array of strings in Fortran is a 2D array of chars in C
            self.array_spec.append(':')
            self.array_dims += 1
            self.required_sizes += 1
            #Need to pass C pointer by value
            self.cintent = 'IN'
        elif param_type.startswith('LOGICAL'):
            self.var_type = Parameter.LOGICAL
        elif param_type.startswith('TYPE'):
            self.var_type = Parameter.CUSTOM_TYPE
            self.type_name = type_pt2.replace('(','').replace(')','')
            if self.array_dims == 0 and self.intent == 'INOUT':
                #Should actually be in, as we need to pass the pointer by value
                self.cintent = 'IN'
        else:
            sys.stderr.write("Error: Unknown type %s for routine %s\n" % (param_type,routine.name))
            self.var_type = None
            self.type_name = param_type
        self._set_size_list()
        self._set_conversion_lines()

    def _set_size_list(self):
        """
        Get the list of dimension sizes for an array, as constants or variable names

        Sets the size_list and required_size_list properties
        required_size_list does not include any dimensions that are constant
        """
        self.size_list = []
        self.required_size_list = []
        i=0
        for dim in self.array_spec:
            if dim == ':':
                if self.required_sizes == 1:
                    self.size_list.append('%sSize' % (self.name))
                elif self.var_type == Parameter.CHARACTER:
                    try:
                        self.size_list.append(['%sNumStrings' % self.name, '%sStringLength' % self.name][i])
                    except IndexError:
                        raise ValueError, ">2D arrays of strings not supported"
                else:
                    self.size_list.append('%sSize%d' % (self.name,i+1))
                i += 1
                self.required_size_list.append(self.size_list[-1])
            else:
                self.size_list.append(dim)

    def _set_conversion_lines(self):
        """
        Get any pointer or string conversions/checks as well as extra local variables required

        Sets pre_call and post_call properties, which are lines to add before and after calling
        the routine in opencmiss.f90
        """
        self.local_variables = []
        self.pre_call = []
        self.post_call = []

        #CMISS Types
        if self.var_type == Parameter.CUSTOM_TYPE:
            if self.array_dims > 0:
                self.local_variables.append('TYPE(%s), POINTER :: %s(%s)' % (self.type_name,self.name,','.join([':']*self.array_dims)))
                self.local_variables.append('INTEGER(C_INT) :: Err')
            else:
                self.local_variables.append('TYPE(%s), POINTER :: %s' % (self.type_name,self.name))
            # If we're in a CMISS...TypeInitialise routine, then objects get allocated
            # in Fortran and we need to convert pointers to C pointers before returning them.
            # For all other routines the pointer to the buffer object doesn't change so we
            # ignore the intent and always check for association before calling the Fortran routine.
            if self.routine.name.endswith('TypeInitialise'):
                self.local_variables.append('INTEGER(C_INT) :: Err')
                self.pre_call.extend(('IF(C_ASSOCIATED(%sPtr)) THEN' % self.name,
                    '%s = CMISSPointerNotNULL' % self.routine.c_f90_name,
                    'ELSE',
                    'NULLIFY(%s)' % self.name,
                    'ALLOCATE(%s, STAT = Err)' % self.name,
                    'IF(Err /= 0) THEN',
                    '%s = CMISSCouldNotAllocatePointer' % self.routine.c_f90_name,
                    'ELSE'))

                self.post_call.extend(('%sPtr=C_LOC(%s)' % (self.name,self.name),
                    'ENDIF',
                    'ENDIF'))
            elif self.routine.name.endswith('TypeFinalise'):
                self.pre_call.extend(('IF(C_ASSOCIATED(%sPtr)) THEN' % self.name,
                    'CALL C_F_POINTER(%sPtr,%s)' % (self.name,self.name),
                    'IF(ASSOCIATED(%s)) THEN' % self.name))

                self.post_call.extend(('DEALLOCATE(%s)' % self.name,
                    '%sPtr = C_NULL_PTR' % self.name,
                    'ELSE',
                    '%s = CMISSErrorConvertingPointer' % self.routine.c_f90_name,
                    'ENDIF',
                    'ELSE',
                    '%s = CMISSPointerIsNULL' % self.routine.c_f90_name,
                    'ENDIF'))
            else:
                if self.array_dims > 0:
                    self.pre_call.extend(('IF(C_ASSOCIATED(%sPtr)) THEN' % self.name,
                        'ALLOCATE(%s(%sSize),STAT=Err)' % (self.name,self.name),
                        'IF(Err == 0) THEN'))
                    if self.intent == 'IN':
                        # Passing an array of CMISS Types to Fortran
                        self.pre_call.append('CALL %ssCopy(%s,%sSize,%sPtr,%s)' % (self.type_name,self.name,self.name,self.name,self.routine.c_f90_name))
                    else:
                        # Getting an array of CMISS Types from Fortran and setting an array of C pointers
                        self.local_variables.append('INTEGER(C_INT) :: %sIndex' % self.name)
                        self.local_variables.append('TYPE(C_PTR), POINTER :: %sCPtrs(:)' % self.name)
                        self.post_call.extend(('CALL C_F_POINTER(%sPtr,%sCPtrs,[%s])' % (self.name,self.name,','.join(self.size_list)),
                            'DO %sIndex=1,%sSize' % (self.name,self.name),
                            '%sCPtrs(%sIndex) = C_LOC(%s(%sIndex))' % ((self.name,) * 4),
                            'ENDDO'))
                    self.post_call.extend(('ELSE',
                        '%s = CMISSCouldNotAllocatePointer' % self.routine.c_f90_name,
                        'ENDIF',
                        'ELSE',
                        '%s = CMISSPointerIsNULL' % self.routine.c_f90_name,
                        'ENDIF'))
                else:
                    self.pre_call.extend(('IF(C_ASSOCIATED(%s)) THEN' % self.c_f90_name(),
                        'CALL C_F_POINTER(%sPtr,%s)' % (self.name, self.name),
                        'IF(ASSOCIATED(%s)) THEN' % self.name))
                    self.post_call.extend(('ELSE',
                        '%s = CMISSErrorConvertingPointer' % self.routine.c_f90_name,
                        'ENDIF',
                        'ELSE',
                        '%s = CMISSPointerIsNULL' % self.routine.c_f90_name,
                        'ENDIF'))
        #Character arrays
        elif self.var_type == Parameter.CHARACTER:
            if self.array_dims > 1:
                #Fortran array has one less dimension
                char_sizes = '(%s)' % ','.join(self.size_list[:-1])
            else:
                char_sizes = ''
            self.local_variables.append('CHARACTER(LEN=%s-1) :: Fortran%s%s' % (self.size_list[-1],self.name,char_sizes))
            self.local_variables.append('%s, POINTER :: %sCChars(%s)' % (Parameter.F90TYPES[self.var_type],self.name,','.join(':'*self.array_dims)))
            if self.intent == 'IN':
                #reverse to account for difference in storage order
                self.pre_call.append('CALL C_F_POINTER(%s,%sCChars,[%s])' % (self.name,self.name,','.join(reversed(self.size_list))))
                if self.array_dims > 1:
                    char_sizes = '(%s)' % ','.join(self.size_list[:-1])
                    self.pre_call.append('CALL CMISSC2FStrings(%sCChars,Fortran%s)' % (self.name,self.name))
                else:
                    char_sizes = ''
                    self.pre_call.append('CALL CMISSC2FString(%sCChars,Fortran%s)' % (self.name,self.name))
            else:
                if self.array_dims > 1:
                    raise ValueError, "output of strings >1D not implemented"
                self.post_call.append('CALL C_F_POINTER(%s,%sCChars,[%s])' % (self.name,self.name,self.size_list[0]))
                self.post_call.append('CALL CMISSF2CString(Fortran%s,%sCChars)' % (self.name,self.name))
        #Arrays of floats, integers or logicals
        elif self.array_dims > 0:
            if self.var_type == Parameter.CHARACTER:
                self.local_variables.append('CHARACTER, POINTER :: %s(%s)' % (self.name,','.join([':']*self.array_dims)))
            else:
                self.local_variables.append('%s, POINTER :: %s(%s)' % (Parameter.F90TYPES[self.var_type],self.name,','.join([':']*self.array_dims)))
            if self.pointer == True and self.cintent == 'OUT':
                #we are setting the value of a pointer
                if self.var_type == Parameter.LOGICAL:
                    self.pre_call.append('NULLIFY(%sLogical)' % self.name)
                else:
                    self.pre_call.append('NULLIFY(%s)' % self.name)
                self.post_call.extend(('%sPtr = C_LOC(%s(1))' % (self.name,self.name),
                    '%sSize = SIZE(%s,1)' % (self.name,self.name),
                    'IF(.NOT.C_ASSOCIATED(%sPtr)) THEN' % self.name,
                    '%s = CMISSErrorConvertingPointer' % self.routine.c_f90_name,
                    'ENDIF'))
            else:
                #pointer is pointing to allocated memory that is being set
                self.pre_call.extend(('IF(C_ASSOCIATED(%s)) THEN' % self.c_f90_name(),
                    'CALL C_F_POINTER(%s,%s,[%s])' % (self.c_f90_name(),self.name,','.join(self.size_list)),
                    'IF(ASSOCIATED(%s)) THEN' % self.name))
                self.post_call.extend(('ELSE',
                    '%s = CMISSErrorConvertingPointer' % self.routine.c_f90_name,
                    'ENDIF',
                    'ELSE',
                    '%s = CMISSPointerIsNULL' % self.routine.c_f90_name,
                    'ENDIF'))
        # Convert from logicals to C integers
        if self.var_type == Parameter.LOGICAL and self.intent != 'IN':
            if self.array_dims > 0:
                #todo if ever required: support more than one dimension
                self.local_variables.append('LOGICAL, POINTER :: %sLogical(%s)' % (self.name,','.join([':']*self.array_dims)))
                self.local_variables.append('INTEGER(C_INT) :: %sLogicalIndex' % (self.name))
                post_call = ['DO %sLogicalIndex=1,SIZE(%sLogical,1)' % (self.name,self.name),
                    'IF(%sLogical(%sLogicalIndex)) THEN' % (self.name,self.name),
                    '%s(%sLogicalIndex) = CMISSTrue' % (self.name,self.name),
                    'ELSE',
                    '%s(%sLogicalIndex) = CMISSFalse' % (self.name,self.name),
                    'ENDIF',
                    'ENDDO']
            else:
                self.local_variables.append('LOGICAL :: %sLogical' % (self.name))
                post_call = ['IF(%sLogical) THEN' % self.name,
                    '%s = CMISSTrue' % self.name,
                    'ELSE',
                    '%s = CMISSFalse' % self.name,
                    'ENDIF']
            if self.pointer == True and self.intent == 'OUT':
                post_call = ['ALLOCATE(%s(SIZE(%sLogical,1)))' % (self.name,self.name)]+post_call
            self.post_call = post_call+self.post_call
        # Convert from C integers to logicals
        elif self.var_type == Parameter.LOGICAL and self.intent == 'IN':
            if self.array_dims > 0:
                #todo if ever required: support more than one dimension
                self.local_variables.append('LOGICAL, POINTER :: %sLogical(%s)' % (self.name,','.join([':']*self.array_dims)))
                self.local_variables.append('INTEGER(C_INT) :: %sLogicalIndex' % (self.name))
                self.pre_call.extend(['ALLOCATE(%sLogical(%s))' % (self.name,','.join(self.size_list)),
                    'DO %sLogicalIndex=1,%sSize' % (self.name,self.name),
                    'IF(%s(%sLogicalIndex) == CMISSTrue) THEN' % (self.name,self.name),
                    '%sLogical(%sLogicalIndex) = .TRUE.' % (self.name,self.name),
                    'ELSE',
                    '%sLogical(%sLogicalIndex) = .FALSE.' % (self.name,self.name),
                    'ENDIF',
                    'ENDDO'])
            else:
                self.local_variables.append('LOGICAL :: %sLogical' % (self.name))
                self.pre_call.extend(['IF(%s == CMISSTrue) THEN' % self.name,
                    '%sLogical = .TRUE.' % self.name,
                    'ELSE',
                    '%sLogical = .FALSE.' % self.name,
                    'ENDIF'])

    def c_f90_name(self):
        """
        Return the name of the parameter as used by the routine in opencmiss_c.f90
        """
        c_f90_name = self.name
        if self.var_type == Parameter.CUSTOM_TYPE or \
                (self.array_dims > 0 and self.var_type != Parameter.CHARACTER):
            c_f90_name += 'Ptr'
        return c_f90_name

    def c_f90_names(self):
        """
        Return param name + name of size param if it exists, separated by a comma for use in the
        function declaration in Fortran
        """
        return ','.join([s for s in self.required_size_list]+[self.c_f90_name()])

    def c_f90_declaration(self):
        """
        Return the parameter declaration for use in the subroutine in opencmiss_c.f90
        """
        c_f90_name=self.c_f90_name()
        output = ''

        #pass by value?
        if self.cintent == 'IN':
            value = 'VALUE, '
        else:
            value = ''

        #possible size parameter
        if self.pointer == True and self.cintent == 'OUT':
            size_type = 'INTEGER(C_INT), INTENT(OUT)'
        else:
            size_type = 'INTEGER(C_INT), VALUE, INTENT(IN)'
        output += ''.join([size_type+' :: '+size_name+'\n    ' for size_name in self.required_size_list])

        if self.array_dims > 0:
            output += 'TYPE(C_PTR), %sINTENT(%s) :: %s' % (value,self.cintent,c_f90_name)
        else:
            output += '%s, %sINTENT(%s) :: %s' % (Parameter.F90TYPES[self.var_type],value,self.cintent,c_f90_name)
        return output

    def call_name(self):
        """
        Return the parameter name to be used in the opencmiss.f90 subroutine call

        Used to pass the name of a converted variable
        """
        output = self.name
        if self.var_type == Parameter.LOGICAL:
            output += 'Logical'
        elif self.var_type == Parameter.CHARACTER:
            output = 'Fortran'+output
        return output

    def to_c(self):
        """
        Calculate C parameter declaration for opencmiss.h

        For arrays, return two parameters with the first being the size
        """
        param = self.name
        #pointer argument?
        if self.array_dims > 0 or self.var_type == Parameter.CHARACTER \
                or self.cintent == 'OUT':
            param = '*'+param
        if self.cintent == 'OUT' and self.pointer == True:
            #add another * as we need a pointer to a pointer, to modify the pointer value
            param = '*'+param
        #parameter type
        if self.var_type != Parameter.CUSTOM_TYPE:
            param = Parameter.CTYPES[self.var_type]+' '+param
        else:
            param = self.type_name+' '+param
        #const?
        if self.intent == 'IN':
            param = 'const '+param
        #size?
        if self.pointer == True and self.cintent == 'OUT':
            #Size is an output
            size_type = 'int *'
        else:
            size_type = 'const int '
        return tuple([size_type+size_name for size_name in self.required_size_list]+[param])


class Type(object):
    """
    Information on a Fortran type
    """
    def __init__(self,name,lineno,line):
        """
        Initialise type

        Arguments:
        name -- Type name
        lineno -- Line number in source where this is defined
        line -- Contents of first line where this type is defined
        """
        self.name = name
        self.lineno = lineno
        self.line = line
        self.lines = [line]


def _join_lines(source):
    """
    Remove Fortran line continuations
    """
    return re.sub(r'[\t ]*&[\t ]*\n[\t ]*&[\t ]*',' ',source)

def _fix_length(line,max_length=132):
    """
    Add Fortran line continuations to break up long lines

    Tries to put the line continuation after a comma
    """
    max_length=132
    line = line.rstrip()
    #account for comments
    commentsplit=line.split('!')
    if len(commentsplit) < 2:
        content = line
        comment = ''
    else:
        content = commentsplit[0]
        comment = '!'.join(commentsplit[1:])
    if content.strip() == '':
        return line
    remaining_content = content
    indent = _get_indent(content)
    content = ''
    while len(remaining_content) > max_length:
        break_pos = remaining_content.rfind(',',0,130)+1
        if break_pos < 0:
            sys.stderr.write("Error: Couldn't truncate line: %s\n" % line)
            exit(1)
        content += remaining_content[0:break_pos]+' &\n'
        remaining_content = indent+'  & '+remaining_content[break_pos:]
    content = content + remaining_content
    if comment:
        content = content+'!'+comment
    return content

def _get_indent(line):
    """
    Return the indentation in front of a line
    """
    indent = line.replace(line.lstrip(),'')
    return indent

def _chain_iterable(iterables):
    """
    Implement itertools.chain.from_iterable to be compatible with Python 2.5
    """
    for it in iterables:
        for element in it:
            yield element

def _indent_lines(lines,indent_size=2,initial_indent=4):
    """
    Indent function content to show nesting of if, else and endif statements
    """
    output = ''
    indent = 0
    for line in lines:
        if line.startswith('ELSE') \
            or line.startswith('ENDIF') \
            or line.startswith('ENDDO'):
            indent -= 1
        output += ' '*(initial_indent+indent*indent_size)+line+'\n'
        if line.startswith('IF') \
            or line.startswith('ELSE') \
            or line.startswith('DO'):
            indent += 1
    return output


if __name__ == '__main__':
    import os
    if len(sys.argv) == 4:
        (cm_path,opencmiss_h_path,opencmiss_c_f90_path) = sys.argv[1:]
    else:
        sys.stderr.write('Usage: %s cm_path opencmiss_h_path opencmiss_c_f90_path\n' % sys.argv[0])
        exit(1)

    cm_source_path=cm_path+os.sep+'src'
    sources = [cm_source_path+os.sep+file_name \
            for file_name in os.listdir(cm_source_path) \
            if file_name.endswith('.f90') and file_name != 'opencmiss.f90']

    library = LibrarySource(cm_source_path+os.sep+'opencmiss.f90',sources)

    with open(opencmiss_h_path,'w') as opencmissh:
        library.write_c_header(opencmissh)
    with open(opencmiss_c_f90_path,'w') as opencmisscf90:
        library.write_c_f90(opencmisscf90)
