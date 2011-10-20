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
            line_re = re.compile(r'^\s*PUBLIC\s*:*\s*([A-Z0-9_,\s]+)',re.IGNORECASE)

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
            routine.get_parameters()

        #todo: work out which routines belong to classes
        self.unbound_routines = self.public_subroutines

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
        is an option between passing an array or a scalar. All other routines
        are also returned.

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
        self.owner_type = None
        self._get_comments()

    def get_parameters(self):
        """Get details of the subroutine parameters

        Sets the Subroutines parameters property as a list of all parameters, excluding the Err parameter
        """

        def filter_match(string):
            if string is None: return ''
            else: return string.strip()

        self.parameters = []
        match = re.search(r'^\s*(RECURSIVE\s+)?SUBROUTINE\s+([A-Z0-9_]+)\(([A-Z0-9_,\*\s]*)\)',self.lines[0],re.IGNORECASE)
        parameters = [p.strip() for p in match.group(3).split(',')]
        try:
            parameters.remove('Err')
        except ValueError:
            try:
                parameters.remove('err')
            except ValueError:
                sys.stderr.write("Warning: Routine doesn't take Err parameter: %s\n" % self.name)

        for param in parameters:
            param_pattern = r"""
            ^\s*([A-Z_]+\s*(\(([A-Z_=,\*0-9]+)\))?)# parameter type at start of line, followed by possible type parameters in brackets
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

    def _get_comments(self):
        """Sets the comment_lines property, a list of comments above the subroutine definition"""

        self.comment_lines = []
        line_num = self.lineno - 1
        while self.source_file.source_lines[line_num].strip().startswith('!>'):
            self.comment_lines.append(self.source_file.source_lines[line_num].strip()[2:].strip())
            line_num -= 1
        self.comment_lines.reverse()


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
        self.type_name = None
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


class DoxygenGrouping(object):
    """Store a line used for grouping in Doxygen"""

    def __init__(self,lineno,line):
        self.lineno = lineno
        self.line = line.strip()


def _join_lines(source):
    """Remove Fortran line continuations"""

    return re.sub(r'[\t ]*&[\t ]*\n[\t ]*&[\t ]*',' ',source)


