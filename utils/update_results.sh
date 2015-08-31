#!/bin/bash

#update the expected results in the 42 testing points directories.
#works from the current directory and copies all exnode or exelem
#files into the expected_results directory

#check we're in 42 point testing dir
FINDDIR=$(pwd | grep 42TestingPoints)
if [ -z $FINDDIR ]; then
    echo "Must be run in 42TestingPoints directory" 1>&2
    exit 1
fi

#find files that aren't in an "output" directory
find . \
         -type d -name '*.svn*' -prune -o \
         -type d -name 'expected_results' -prune -o \
         -type d -name 'output' -prune -o \
         -type f -name "*.exnode" -execdir cp {} expected_results/ \; -o \
         -type f -name "*.exelem" -execdir cp {} expected_results/ \;
         
#now find files in an "output" directory
find . \
         -type d -name '*.svn*' -prune -o \
         -type d -name 'expected_results' -prune -o \
         -type f -wholename "*/output/*.exnode" -execdir cp {} ../expected_results/ \; -o \
         -type f -wholename "*/output/*.exelem" -execdir cp {} ../expected_results/ \;
         
