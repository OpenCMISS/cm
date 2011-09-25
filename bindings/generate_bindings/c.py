import os
from parse import LibrarySource

def generate(cm_path,args):
    opencmiss_h_path,opencmiss_c_f90_path = args

    library = LibrarySource(cm_path)

    with open(opencmiss_h_path,'w') as opencmissh:
        library.write_c_header(opencmissh)
    with open(opencmiss_c_f90_path,'w') as opencmisscf90:
        library.write_c_f90(opencmisscf90)

