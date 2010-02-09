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
# The Original Code is OpenCMISS
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
  OPENCMISS_ROOT := ../../../..
  GLOBAL_ROOT := $(CURDIR)/../../../..
else
  GLOBAL_ROOT := ${OPENCMISS_ROOT}/cm
endif

#----------------------------------------------------------------------------------------------------------------------------------

EXAMPLE_PATH=FluidMechanics/Darcy/

EXAMPLE_NAME = Static

include $(GLOBAL_ROOT)/NewExampleMakefile

# customisation begins here ---
MAIN_OBJECT_DIR = $(GLOBAL_ROOT)/object/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(MPI)/$(COMPILER)
MAIN_SOURCE_DIR = $(GLOBAL_ROOT)/src

# let's assume fluid_mechanics_IO_routines has been built in main_object, copy them over
# redefine preliminaries
preliminaries: $(OBJECT_DIR) $(EXE_DIR)
	cp $(MAIN_OBJECT_DIR)/fluid_mechanics_IO_routines.o $(OBJECT_DIR)
	cp $(MAIN_OBJECT_DIR)/fluid_mechanics_io_routines.mod $(OBJECT_DIR)

info:
	@echo "OBJECTS:"
	@echo $(OBJECTS)
	@echo "OBJECT_DIR:"
	@echo $(OBJECT_DIR)
	@echo "MAIN_OBJECT_DIR:"
	@echo $(MAIN_OBJECT_DIR)
	@echo "EXAMPLE_SOURCE_DIR:"
	@echo $(EXAMPLE_SOURCE_DIR)
	@echo "MAIN_SOURCE_DIR:"
	@echo $(MAIN_SOURCE_DIR)
