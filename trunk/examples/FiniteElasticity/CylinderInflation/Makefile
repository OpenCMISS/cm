# -*- makefile -*-
#
# For use with GNU make.
#
# $Id: Makefile 27 2007-07-24 16:52:51Z cpb $
#
#----------------------------------------------------------------------------------------------------------------------------------
# Makefile for compiling openCMISS examples
#
# Original by Chris Bradley adapted from the CMISS Makefile by Karl Tomlinson 
# Changes:
#
#----------------------------------------------------------------------------------------------------------------------------------
#
# LICENSE
#
# Version: MPL 1.1/GPL 2.0/LGPL 2.1
#
# The contents of this file are subject to the Mozilla Public License
# Version 1.1 (the "License"); you may not use this file except in
# compliance with the License. You may obtain a copy of the License at
# http://www.mozilla.org/MPL/
#
# Software distributed under the License is distributed on an "AS IS"
# basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
# License for the specific language governing rights and limitations
# under the License.
#
# The Original Code is openCMISS
#
# The Initial Developer of the Original Code is University of Auckland,
# Auckland, New Zealand and University of Oxford, Oxford, United
# Kingdom. Portions created by the University of Auckland and University
# of Oxford are Copyright (C) 2007 by the University of Auckland and
# the University of Oxford. All Rights Reserved.
#
# Contributor(s):
#
# Alternatively, the contents of this file may be used under the terms of
# either the GNU General Public License Version 2 or later (the "GPL"), or
# the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
# in which case the provisions of the GPL or the LGPL are applicable instead
# of those above. If you wish to allow use of your version of this file only
# under the terms of either the GPL or the LGPL, and not to allow others to
# use your version of this file under the terms of the MPL, indicate your
# decision by deleting the provisions above and replace them with the notice
# and other provisions required by the GPL or the LGPL. If you do not delete
# the provisions above, a recipient may use your version of this file under
# the terms of any one of the MPL, the GPL or the LGPL.
#	
#----------------------------------------------------------------------------------------------------------------------------------

MAKEFLAGS = --no-builtin-rules --warn-undefined-variables

#----------------------------------------------------------------------------------------------------------------------------------

ifndef OPENCMISS_ROOT
  OPENCMISS_ROOT := ../..
  GLOBAL_ROOT := $(CURDIR)/../..
else
  GLOBAL_ROOT := ${OPENCMISS_ROOT}/cm
endif

#----------------------------------------------------------------------------------------------------------------------------------


EXAMPLE_PATH = FiniteElasticity/

EXAMPLE_NAME = CylinderInflation

include $(GLOBAL_ROOT)/ExampleMakefile


#----------------------------------------------------------------------------------------------------------------------------------
# Chris is going to hate this, but I need a rule for solid_mechanics_IO_cylinder_geometry.f90 so redefine stuff here
OBJECTS += $(OBJECT_DIR)/solid_mechanics_IO_cylinder_geometry.o
OBJECTS += $(OBJECT_DIR)/CylinderInflationAnalytic.o

$(OBJECT_DIR)/solid_mechanics_IO_cylinder_geometry.o    :       $(EXAMPLE_SOURCE_DIR)/solid_mechanics_IO_cylinder_geometry.f90 \
        $(OBJECT_DIR)/constants.o \
        $(OBJECT_DIR)/computational_environment.o \
        $(OBJECT_DIR)/kinds.o \
        $(OBJECT_DIR)/solid_mechanics_IO_routines.o \
        $(OBJECT_DIR)/strings.o \
        $(OBJECT_DIR)/types.o \

$(OBJECT_DIR)/CylinderInflationAnalytic.o    :       $(EXAMPLE_SOURCE_DIR)/CylinderInflationAnalytic.f90 \
        $(OBJECT_DIR)/constants.o \
        $(OBJECT_DIR)/kinds.o \
        $(OBJECT_DIR)/types.o \

$(OBJECT_DIR)/$(BASE_EXAMPLE_NAME).o    :       $(EXAMPLE_SOURCE_DIR)/$(BASE_EXAMPLE_NAME).f90 \
        $(OBJECT_DIR)/analytic_analysis_routines.o \
        $(OBJECT_DIR)/base_routines.o \
        $(OBJECT_DIR)/basis_routines.o \
        $(OBJECT_DIR)/boundary_condition_routines.o \
        $(OBJECT_DIR)/cmiss.o \
        $(OBJECT_DIR)/cmiss_mpi.o \
        $(OBJECT_DIR)/computational_environment.o \
        $(OBJECT_DIR)/constants.o \
        $(OBJECT_DIR)/control_loop_routines.o \
        $(OBJECT_DIR)/coordinate_routines.o \
        $(OBJECT_DIR)/CylinderInflationAnalytic.o \
        $(OBJECT_DIR)/distributed_matrix_vector.o \
        $(OBJECT_DIR)/domain_mappings.o \
        $(OBJECT_DIR)/equations_routines.o \
        $(OBJECT_DIR)/equations_set_constants.o \
        $(OBJECT_DIR)/equations_set_routines.o \
        $(OBJECT_DIR)/field_routines.o \
        $(OBJECT_DIR)/field_IO_routines.o \
        $(OBJECT_DIR)/input_output.o \
        $(OBJECT_DIR)/iso_varying_string.o \
        $(OBJECT_DIR)/kinds.o \
        $(OBJECT_DIR)/lists.o \
        $(OBJECT_DIR)/mesh_routines.o \
        $(OBJECT_DIR)/node_routines.o \
        $(OBJECT_DIR)/problem_constants.o \
        $(OBJECT_DIR)/problem_routines.o \
        $(OBJECT_DIR)/region_routines.o \
        $(OBJECT_DIR)/solid_mechanics_IO_routines.o \
        $(OBJECT_DIR)/solid_mechanics_IO_cylinder_geometry.o \
        $(OBJECT_DIR)/solver_routines.o \
        $(OBJECT_DIR)/test_framework_routines.o \
        $(OBJECT_DIR)/timer_f.o \
        $(OBJECT_DIR)/types.o 

