#! /bin/tcsh -f
#
# This shell file will call the given makefile with the list of files specified
# on the command line, once per file, with the extension of the file changed 
# from <new_extension>.
#
# Usage: 
#   call_makefile.sh makefile new_extension list_of_files
# Created: 
#   Glen Harris, 8 August 1996
# Updates: 
#   Glen Harris, 16 August 1996  Removed duplicates from call
#
set makefile_name=$1
shift
set new_extension=$1
shift
set temp_file=call_makefile.tmp
if ${1}x != x then
	rm -f ${temp_file}
	foreach file ($*)
		echo ${file:r}.${new_extension} >> ${temp_file}
	end
#       Kill the 'Target is up to date' warning output by make
	make -f ${makefile_name} `sort -u ${temp_file}` |& grep -v "is up to date"
	rm -f ${temp_file}
endif

