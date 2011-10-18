import os
from parse import *
from c import subroutine_c_name

def generate_docstring(routine):
    docstring = '\n'.join(routine.comment_lines)
    docstring += '\n\n'
    routine.get_parameters()
    for param in routine.parameters:
        docstring += ':param %s: %s\n' % (param.name, param.doxygen)
    return docstring.strip()


def generate(cm_path,args):
    """
    Create a list of CMISS types for the Python bindings to use
    and generate docstrings for routines using the doxygen comments
    and argument lists
    """
    types_file = open(os.sep.join((cm_path,'bindings','python','opencmiss','_types.py')),'w')
    docstrings_file = open(os.sep.join((cm_path,'bindings','python','opencmiss','_docstrings.py')),'w')

    library = LibrarySource(cm_path)

    types_file.write('types = [\n')
    docstrings_file.write('docstrings = {\n')

    for pub in library.public_objects.values():
        if isinstance(pub,Type):
            if pub.name.startswith('CMISS') and pub.name.endswith('Type'):
                types_file.write("  '"+pub.name[5:-4]+"',\n")

    for routine in library.public_subroutines:
        docstrings_file.write('  "%s":"""%s""",\n' % (subroutine_c_name(routine)[0], generate_docstring(routine)))

    types_file.write(']\n')
    docstrings_file.write('}\n')

    types_file.close()
    docstrings_file.close()

