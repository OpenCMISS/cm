#!/bin/sh
#
# Shell file for generating fig files (for the first time) from
# gnu files.
#
# Usage:
#   genfigs *.gnu
# Created:
#   Chris Bradley 9 March 1996
# Updates:
#
for filename
do
	${OPENCMISS_ROOT}/cm/doc/latex/genfigs1.sh $filename
done
