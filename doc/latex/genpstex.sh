#!/bin/sh
#
# Shell file for generating pstex files (for the first time)
#
# Usage:
#   genpstex [figs]/plots *.fig 
# Created:
#   Chris Bradley
# Updates:
#   Chris Bradley 10/3/96 Added figs/plots option
#
type=$1
shift
for filename
do
	${OPENCMISS_ROOT}/cm/doc/latex/genpstex1.sh $type $filename
done
