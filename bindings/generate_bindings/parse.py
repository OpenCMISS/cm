from __future__ import with_statement
import sys, os
import re
from operator import attrgetter

class LibrarySource(object):
    """Holds info on all the library source code"""

    class SourceFile(object):
        """Info for an individual source file"""

        class SectionFinder(object):
            """Match a section within a source file"""

            def __init__(self,source_file):
                self.match = None
                self.lineno = 0
                self.lines = []
                self.source_file = source_file

            def check_for_end(self,line):
                if self.end_re.search(line):
                    self.finish()
                    self.lines = []
                    return True
                return False

            def check_for_start(self,lineno,line):
                match = self.start_re.search(line)
                if match:
                    self.match = match
                    self.lineno = lineno
                    self.lines.append(line)
                    return True
                return False

        class LineFinder(object):
            """Match a line within a source file"""

            def __init__(self,source_file):
                self.source_file = source_file

            def check_match(self,line,lineno):
                match = self.line_re.search(line)
                if match:
                    self.add(match,lineno)

        class SubroutineFinder(SectionFinder):
            start_re = re.compile(r'^\s*(RECURSIVE\s+)?SUBROUTINE\s+([A-Z0-9_]+)\(',re.IGNORECASE)
            end_re = re.compile(r'^\s*END\s*SUBROUTINE',re.IGNORECASE)

            def finish(self):
                name = self.match.group(2)
                self.source_file.subroutines[name] = Subroutine(name,self.lineno,self.lines,self.source_file)

        class InterfaceFinder(SectionFinder):
            start_re = re.compile(r'^\s*INTERFACE\s+([A-Z0-9_]+)',re.IGNORECASE)
            end_re = re.compile(r'^\s*END\s*INTERFACE',re.IGNORECASE)

            def finish(self):
                name = self.match.group(1)
                self.source_file.interfaces[name] = Interface(name,self.lineno,self.lines,self.source_file)

        class TypeFinder(SectionFinder):
            start_re = re.compile(r'^\s*TYPE\s+([A-Z0-9_]+)',re.IGNORECASE)
            end_re = re.compile(r'^\s*END\s*TYPE',re.IGNORECASE)

            def finish(self):
                name = self.match.group(1)
                self.source_file.types[name] = Type(name,self.lineno,self.lines,self.source_file)

        class PublicFinder(LineFinder):
            line_re = re.compile(r'^\s*PUBLIC\s+([A-Z0-9_,\s]+)',re.IGNORECASE)

            def add(self,match,lineno):
                for symbol in match.group(1).split(','):
                    self.source_file.public.append(symbol.strip())

        class ConstantFinder(LineFinder):
            line_re = re.compile(r'^\s*INTEGER\([A-Z0-9\(\),_\s]+::\s*([A-Z0-9_]+)\s*=\s*([A-Z0-9_\-\.]+)[^!]*(!<.*$)?',re.IGNORECASE)

            def add(self,match,lineno):
                name = match.group(1)
                assignment = match.group(2)
                if match.group(3) is None:
                    doxy = ''
                else:
                    doxy = match.group(3)[2:].strip()
                self.source_file.constants[name] = Constant(name,lineno,assignment,doxy)

        class DoxygenGroupingFinder(LineFinder):
            #match at least one whitespace character before the ! to make sure
            #we don't get stuff from the file header
            line_re = re.compile(r'^\s+!\s*>\s*(\\(addtogroup|brief|see)|@[\{\}])(.*$)',re.IGNORECASE)

            def add(self,match,lineno):
                line = match.group(1)
                if match.group(3) is not None:
                    line += match.group(3)
                self.source_file.doxygen_groupings.append(DoxygenGrouping(lineno,line))

        def __init__(self,source_file,params_only=False):
            """Initialise SourceFile object

            Arguments:
            source_file -- Path to the source file
            """

            self.file_path = source_file
            self.public = []
            self.doxygen_groupings = []
            self.interfaces = {}
            self.subroutines = {}
            self.constants = {}
            self.types = {}
            self.parse_file(params_only)

        def parse_file(self,params_only=False):
            """Run through file once, getting everything we'll need"""

            source_lines = _join_lines(open(self.file_path,'r').read()).splitlines(True)
            if not params_only:
                #only keep the source_lines if we need them
                self.source_lines = source_lines

            #Set the things we want to find
            line_finders = []
            section_finders = []
            line_finders.append(self.ConstantFinder(self))
            if not params_only:
                line_finders.extend((
                    self.PublicFinder(self),
                    self.DoxygenGroupingFinder(self)))
                section_finders.extend((
                    self.SubroutineFinder(self),
                    self.InterfaceFinder(self),
                    self.TypeFinder(self)))

            #Find them
            current_section = None
            for (lineno,line) in enumerate(source_lines):
                if current_section is not None:
                    current_section.lines.append(line)
                    if current_section.check_for_end(line):
                        current_section = None
                else:
                    for line_finder in line_finders:
                        line_finder.check_match(line,lineno)

                    for section in section_finders:
                        if section.check_for_start(lineno,line):
                            current_section = section
                            break

    def __init__(self,cm_path):
        """Load library information from source files

        Arguments:
        cm_path -- Path to OpenCMISS cm directory
        """

        self.lib_source = self.SourceFile(os.sep.join((cm_path,'src','opencmiss.f90')))
        cm_source_path = cm_path+os.sep+'src'
        source_files = [cm_source_path+os.sep+file_name \
                for file_name in os.listdir(cm_source_path) \
                if file_name.endswith('.f90') and file_name != 'opencmiss.f90']
        self.sources = [self.SourceFile(source,params_only=True) for source in source_files]

        self.resolve_constants()

        #Get all public types, constants and routines to include
        #Store all objects to be output in a dictionary with line number as key
        self.public_objects = {}
        for t in self.lib_source.types.values():
            if t.name in self.lib_source.public:
                self.public_objects[t.lineno] = t

        for const in self.lib_source.constants.values():
            if const.name in self.lib_source.public:
                self.public_objects[const.lineno] = const

        self.public_subroutines=[routine for routine in self.lib_source.subroutines.values() \
            if routine.name in self.lib_source.public]

        for interface in self.lib_source.interfaces.values():
            if interface.name in self.lib_source.public:
                self.public_subroutines += [self.lib_source.subroutines[routine] \
                        for routine in interface.get_subroutines()]

        self.public_subroutines=sorted(self.public_subroutines,key=attrgetter('name'))
        #Remove CMISS...TypesCopy routines, as these are only used within the C bindings
        #Also remove CMISSGeneratedMeshSurfaceGet for now as it takes an allocatable array but will be removed soon anyways.
        self.public_subroutines = filter(lambda r: not (r.name.startswith('CMISSGeneratedMeshSurfaceGet') or r.name.endswith('TypesCopy')),self.public_subroutines)

        for routine in self.public_subroutines:
            self.public_objects[routine.lineno] = routine

        for doxygen_grouping in self.lib_source.doxygen_groupings:
            self.public_objects[doxygen_grouping.lineno] = doxygen_grouping

    def resolve_constants(self):
        """Go through all public constants and work out their actual values"""

        for pub in self.lib_source.public:
            if self.lib_source.constants.has_key(pub):
                self.get_constant_value(pub)

    def get_constant_value(self,constant):
        """Get the actual value for a constant from the source files

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
                if i == (len(self.sources) - 1):
                    exhausted = True
        if not self.lib_source.constants[constant].resolved:
            sys.stderr.write("Warning: Couldn't resolve constant value: %s\n" % constant)

    def write_c_header(self,output):
        """Write opencmiss.h containing constants, typedefs and routine declarations

        Arguments:
        output -- File to write to
        """

        output.write('/*\n * opencmiss.h automatically generated from opencmiss.f90\n * Do not edit this file directly, instead edit opencmiss.f90 or generatec.py\n */\n\n' + \
            '#ifndef OPENCMISS_H\n' + \
            '#define OPENCMISS_H\n' + \
            '\n/*\n * Defines\n */\n\n' + \
            'const int CMISSNoError = 0;\n' + \
            'const int CMISSPointerIsNULL = -1;\n' + \
            'const int CMISSPointerNotNULL = -2;\n' + \
            'const int CMISSCouldNotAllocatePointer = -3;\n' + \
            'const int CMISSErrorConvertingPointer = -4;\n\n' + \
            'typedef %s CMISSBool;\n' % _logical_type() + \
            'const CMISSBool CMISSTrue = 1;\n' + \
            'const CMISSBool CMISSFalse = 0;\n\n' + \
            'typedef int CMISSError;\n\n')
        for lineno in sorted(self.public_objects.keys()):
            output.write(self.public_objects[lineno].to_c_header())
        output.write('\n#endif\n')

    def write_c_f90(self,output):
        """Write opencmiss_c.f90 containing Fortran routines

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
    """Information on a public constant"""

    def __init__(self,name,lineno,assignment,doxygen_comment):
        """Initialise Constant

        Arguments:
        name -- Variable name
        assignment -- Value or another variable assigned to this variable
        doxygen_comment -- Contents of the doxygen comment describing the constant
        """

        self.name = name
        self.lineno = lineno
        self.assignment = assignment
        self.doxygen_comment = doxygen_comment
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

    def to_c_header(self):
        """Return the C definition of this constant"""

        if self.resolved:
            if self.doxygen_comment != '':
                return 'const int %s = %d; /*<%s */\n' % (self.name,self.value,self.doxygen_comment)
            else:
                return 'const int %s = %d;\n' % (self.name,self.value)
        else:
            return ''


class Interface(object):
    """Information on an interface"""

    def __init__(self,name,lineno,lines,source_file):
        """Initialise an interface

        Arguments:
        name -- Interface name
        lineno -- Line number where the interface starts
        lines -- Contents of interface as a list of lines
        source_file -- Source file containing the interface
        """

        self.name = name
        self.lineno = lineno
        self.lines = lines
        self.source = source_file

    def get_subroutines(self):
        """Find the subroutines for an interface

        Choose the one with the highest number if there are options. This
        corresponds to the routine that takes array parameters

        Returns a list of subroutines
        """

        all_subroutines = []
        routine_re=re.compile(r'MODULE PROCEDURE ([A-Z0-9_]+)',re.IGNORECASE)
        varying_string_re=re.compile(r'VSC*(Obj|Number|)[0-9]*$',re.IGNORECASE)

        for line in self.lines:
            match = routine_re.search(line)
            if match:
                routine_name = match.group(1)
                if varying_string_re.search(routine_name):
                    #Don't include routines using varying_string parameters
                    pass
                else:
                    all_subroutines.append(routine_name)

        subroutines = self._get_array_routines(all_subroutines)

        for routine in subroutines:
            self.source.subroutines[routine].interface = self

        return subroutines

    def _get_array_routines(self,routine_list):
        """Return a list of the routines that take array parameters if there
        is an option between passing an array or a scalar

        Arguments:
        routine_list -- List of subroutine names
        """

        routine_groups = {}
        routines = []

        #Group routines depending on their name, minus any number indicating
        #whether they take a scalar or array
        for routine in routine_list:
            routine_group = re.sub('\d','0',routine)
            if routine_groups.has_key(routine_group):
                routine_groups[routine_group].append(routine)
            else:
                routine_groups[routine_group] = [routine]

        for group in routine_groups.keys():
            max_number = -1
            for routine in routine_groups[group]:
                try:
                    number = int(filter(str.isdigit,routine))
                    if number > max_number:
                        array_routine = routine
                except ValueError:
                    #only one routine in group
                    array_routine = routine
            routines.append(array_routine)

        return routines


class Subroutine(object):
    """Store information for a subroutine"""

    def __init__(self,name,lineno,lines,source_file):
        self.name = name
        self.lineno = lineno
        self.lines = lines
        self.source_file = source_file
        self.parameters = None
        self.interface = None
        self._set_c_name()
        self._get_comments()

    def get_parameters(self):
        """Get details of the subroutine parameters

        Sets the Subroutines parameters property as a list of all parameters, excluding the Err parameter
        """

        def filter_match(string):
            if string is None: return ''
            else: return string.strip()

        self.parameters = []
        match = re.search(r'^\s*(RECURSIVE\s+)?SUBROUTINE\s+([A-Z0-9_]+)\(([A-Z0-9_,\s]*)\)',self.lines[0],re.IGNORECASE)
        parameters = [p.strip() for p in match.group(3).split(',')]
        try:
            parameters.remove('Err')
        except ValueError:
            sys.stderr.write("Warning: Routine doesn't take Err parameter: %s\n" % self.name)

        for param in parameters:
            param_pattern = r"""
            ^\s*([A-Z_]+\s*(\(([A-Z_=\*0-9]+)\))?) # parameter type at start of line, followed by possible type parameters in brackets
            \s*([A-Z0-9\s_\(\):,\s]+)?\s*          # extra specifications such as intent
            ::
            [A-Z_,\s\(\):]*                        # Allow for other parameters to be included on the same line
            [,\s:]                                 # Make sure we matched the full parameter name
            %s                                     # Parameter name
            (\(([0-9,:]+)\))?                      # Array dimensions if present
            [,\s$]                                 # Whitespace, comma or end of line to make sure we've matched the full parameter name
            [^!]*(!<(.*)$)?                        # Doxygen comment
            """ % param
            param_re = re.compile(param_pattern,re.IGNORECASE|re.VERBOSE)

            for line in self.lines:
                match = param_re.search(line)
                if match:
                    param_type = match.group(1)
                    (type_params,extra_stuff,array,doxygen) = (filter_match(match.group(i)) for i in (3,4,6,8))
                    self.parameters.append(Parameter(param,self,param_type,type_params,extra_stuff,array,doxygen))
                    break
            if not match:
                raise RuntimeError, "Couldn't find parameter %s for subroutine %s" % (param,self.name)

    def _set_c_name(self):
        """Get the name of the routine as used from C and as used in opencmiss_c.f90

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

    def _get_comments(self):
        """Sets the comment_lines property, a list of comments above the subroutine definition"""

        self.comment_lines = []
        line_num = self.lineno - 1
        while self.source_file.source_lines[line_num].strip().startswith('!>'):
            self.comment_lines.append(self.source_file.source_lines[line_num].strip()[2:].strip())
            line_num -= 1
        self.comment_lines.reverse()

    def to_c_header(self):
        """Returns the function declaration in C"""

        (swig_start,swig_end) = self.swig_lines()
        output = '\n'+swig_start

        if self.parameters is None:
            self.get_parameters()
        output += '/*>'
        output += '\n *>'.join(self.comment_lines)
        output += ' */\n'
        output += 'CMISSError %s(' % self.c_name

        c_parameters = _chain_iterable([p.to_c() for p in self.parameters])
        comments = _chain_iterable([p.doxygen_comments() for p in self.parameters])
        output += ',\n    '.join(['%s /*<%s */' % (p,c) for (p,c) in zip(c_parameters,comments)])
        output += ');\n'

        output += swig_end

        return output

    def swig_lines(self):
        """Return lines used before and after subroutine for SWIG interfaces
        """
        if self.name.endswith('TypeInitialise'):
            type = self.name[0:-len('Initialise')]
            name = type[0:-len('Type')]
            start_lines = '#ifdef SWIG\n  %%apply CMISSDummyInitialiseType *CMISSDummy{%s *%s};\n#endif\n' % (type,name)
            end_lines = '#ifdef SWIG\n  %%clear %s *%s;\n#endif\n' % (type,name)
            return (start_lines, end_lines)
        return ('','')

    def local_c_f90_vars(self):
        """Returns a list of extra local variables required for this subroutine, for use in converting
        to and from C variables
        """

        local_variables = [param.local_variables for param in self.parameters]
        return list(_chain_iterable(local_variables))

    def pre_call(self):
        """Returns a list of lines to add before calling the Fortran routine
        for converting to and from C
        """

        pre_lines = [param.pre_call for param in self.parameters]
        return list(_chain_iterable(pre_lines))

    def post_call(self):
        """Returns a list of lines to add after calling the Fortran routine
        for converting to and from C
        """

        #reverse to keep everything in order, as some if statements need to line up
        #from the pre_call lines
        post_lines = [param.post_call for param in reversed(self.parameters)]
        return list(_chain_iterable(post_lines))

    def to_c_f90(self):
        """Returns the C function implemented in Fortran for opencmiss_c.f90"""

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
    """Information on a subroutine parameter"""

    #Parameter types enum:
    (INTEGER, \
    FLOAT, \
    DOUBLE, \
    CHARACTER, \
    LOGICAL, \
    CUSTOM_TYPE) \
        = range(6)

    #Corresponding variable types for C
    CTYPES = ('int','float','double','char','CMISSBool',None)

    #Variable types as used in opencmiss_c.f90
    F90TYPES = ('INTEGER(C_INT)','REAL(C_FLOAT)','REAL(C_DOUBLE)','CHARACTER(LEN=1,KIND=C_CHAR)','LOGICAL','TYPE(C_PTR)')

    def __init__(self,name,routine,param_type,type_params,extra_stuff,array,doxygen):
        """Initialise a parameter

        Arguments:
        name -- Parameter name
        routine -- Pointer back to the subroutine this parameter belongs to
        param_type -- String from the parameter declaration
        type_params -- Any parameters for parameter type, eg "DP" for a real
        extra_stuff -- Any extra parameter properties listed after the type, including intent
        array -- The array dimensions included after the parameter name if they exist, otherwise an empty string
        doxygen -- The doxygen comment after the parameteter
        """

        self.name = name
        self.routine = routine
        self.pointer = False
        self.doxygen = doxygen
        intent = None

        if extra_stuff != '':
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
            self.precision = type_params
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
            self.type_name = type_params
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
        """Get the list of dimension sizes for an array, as constants or variable names

        Sets the size_list, required_size_list and size_doxygen properties
        required_size_list does not include any dimensions that are constant
        size_doxygen has the same length as required_size_list
        """

        self.size_list = []
        self.required_size_list = []
        self.size_doxygen = []
        i=0
        for dim in self.array_spec:
            if dim == ':':
                if self.required_sizes == 1:
                    self.size_list.append('%sSize' % (self.name))
                    if self.var_type == Parameter.CHARACTER:
                        self.size_doxygen.append('Length of %s string' % self.name)
                    else:
                        self.size_doxygen.append('Length of %s' % self.name)
                elif self.var_type == Parameter.CHARACTER:
                    try:
                        self.size_list.append(['%sNumStrings' % self.name, '%sStringLength' % self.name][i])
                        self.size_doxygen.append(['Number of strings in %s' % self.name, 'Length of strings in %s' % self.name][i])
                    except IndexError:
                        raise ValueError, ">2D arrays of strings not supported"
                else:
                    self.size_list.append('%sSize%d' % (self.name,i+1))
                    self.size_doxygen.append('Size of dimension %d of %s' % (i+1,self.name))
                i += 1
                self.required_size_list.append(self.size_list[-1])
            else:
                self.size_list.append(dim)

    def _set_conversion_lines(self):
        """Get any pointer or string conversions/checks as well as extra local variables required

        Sets pre_call and post_call properties, which are lines to add before and after calling
        the routine in opencmiss.f90
        """

        self.local_variables = []
        self.pre_call = []
        self.post_call = []

        #CMISS Types
        if self.var_type == Parameter.CUSTOM_TYPE:
            if self.array_dims > 0:
                self.local_variables.append('TYPE(%s), TARGET :: %s(%s)' % (self.type_name,self.name,','.join(self.size_list)))
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
                    self.pre_call.append('IF(C_ASSOCIATED(%sPtr)) THEN' % self.name)
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


    def c_f90_name(self):
        """Return the name of the parameter as used by the routine in opencmiss_c.f90"""

        c_f90_name = self.name
        if self.var_type == Parameter.CUSTOM_TYPE or \
                (self.array_dims > 0 and self.var_type != Parameter.CHARACTER):
            c_f90_name += 'Ptr'
        return c_f90_name

    def c_f90_names(self):
        """Return param name + name of size param if it exists,
        separated by a comma for use in the function declaration in Fortran
        """

        return ','.join([s for s in self.required_size_list]+[self.c_f90_name()])

    def c_f90_declaration(self):
        """Return the parameter declaration for use in the subroutine in opencmiss_c.f90"""

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
        """Return the parameter name to be used in the opencmiss.f90 subroutine call

        Used to pass the name of a converted variable
        """

        output = self.name
        if self.var_type == Parameter.CHARACTER:
            output = 'Fortran'+output
        return output

    def to_c(self):
        """Calculate C parameter declaration for opencmiss.h

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

    def doxygen_comments(self):
        """Return a list of doxygen comments corresponding to the list of
        parameters returned by to_c
        """

        return self.size_doxygen+[self.doxygen]


class Type(object):
    """Information on a Fortran type"""

    def __init__(self,name,lineno,lines,source_file):
        """Initialise type

        Arguments:
        name -- Type name
        lineno -- Line number in source where this is defined
        lines -- Contents of lines where this type is defined
        """

        self.name = name
        self.lineno = lineno
        self.lines = lines
        self.source_file = source_file
        self._get_comments()

    def _get_comments(self):
        """Sets the comment_lines property, a list of comments above the type definition"""

        self.comment_lines = []
        line_num = self.lineno - 1
        while self.source_file.source_lines[line_num].strip().startswith('!>'):
            self.comment_lines.append(self.source_file.source_lines[line_num].strip()[2:].strip())
            line_num -= 1
        self.comment_lines.reverse()

    def to_c_header(self):
        """Return the struct and typedef definition for use in opencmiss.h"""

        output = 'struct %s_;\n' % self.name
        output += '/*>'
        output += '\n *>'.join(self.comment_lines)
        output += ' */\n'
        output += 'typedef struct %s_ *%s;\n\n' % (self.name,self.name)
        return output


class DoxygenGrouping(object):
    """Store a line used for grouping in Doxygen"""

    def __init__(self,lineno,line):
        self.lineno = lineno
        self.line = line.strip()

    def to_c_header(self):
        """Return the doxygen comment for use in opencmiss.h"""
        return '/*>'+self.line+' */\n'


def _join_lines(source):
    """Remove Fortran line continuations"""

    return re.sub(r'[\t ]*&[\t ]*\n[\t ]*&[\t ]*',' ',source)


def _fix_length(line,max_length=132):
    """Add Fortran line continuations to break up long lines

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
    """Return the indentation in front of a line"""

    indent = line.replace(line.lstrip(),'')
    return indent


def _chain_iterable(iterables):
    """Implement itertools.chain.from_iterable so we can use it in Python 2.5"""

    for it in iterables:
        for element in it:
            yield element


def _indent_lines(lines,indent_size=2,initial_indent=4):
    """Indent function content to show nesting of if, else and endif statements"""

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

def _logical_type():
    """Return the C type to match Fortran logical type depending on the compiler used"""
    #Both ifortran and gfortran use 4 bytes
    #uint32_t is optional so might not be defined
    return "unsigned int"


