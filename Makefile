# -*- makefile -*-
#
# For use with GNU make.
#
#
#----------------------------------------------------------------------------------------------------------------------------------
# Makefile for compiling OpenCMISS library
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
# Auckland, New Zealand, the University of Oxford, Oxford, United
# Kingdom and King's College, London, United Kingdom. Portions created
# by the University of Auckland, the University of Oxford and King's
# College, London are Copyright (C) 2007-2010 by the University of
# Auckland, the University of Oxford and King's College, London.
# All Rights Reserved.
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
  OPENCMISS_ROOT = $(CURDIR)/../
endif

GLOBAL_CM_ROOT = $(CURDIR)

include $(GLOBAL_CM_ROOT)/utils/MakefileCommon.inc

#----------------------------------------------------------------------------------------------------------------------------------

SOURCE_DIR = $(GLOBAL_CM_ROOT)/src
OBJECT_DIR := $(MAIN_OBJECT_DIR)
BASE_LIB_NAME = OpenCMISS
MODULE_DIR := $(OBJECT_DIR)
MOD_INC_NAME := opencmiss.mod
MOD_INCLUDE := $(INC_DIR)/$(MOD_INC_NAME)
MOD_SOURCE_INC := $(OBJECT_DIR)/$(MOD_INC_NAME)
HEADER_INC_NAME := opencmiss.h
HEADER_INCLUDE := $(INC_DIR)/$(HEADER_INC_NAME)
C_F90_SOURCE := $(SOURCE_DIR)/opencmiss_c.f90
BINDINGS_DIR = $(GLOBAL_CM_ROOT)/bindings
BINDINGS_GENERATE_SCRIPT := $(BINDINGS_DIR)/generate_bindings
LIB_NAME := lib$(BASE_LIB_NAME)$(EXE_ABI_SUFFIX)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX).a
LIBRARY := $(LIB_DIR)/$(LIB_NAME)

C_INCLUDE_DIRS := $(SOURCE_DIR)
F_INCLUDE_DIRS := $(MODULE_DIR)

CPPFLAGS += $(EXTERNAL_INCLUDE_PATH)
FPPFLAGS += $(EXTERNAL_INCLUDE_PATH)
ELFLAGS += $(EXTERNAL_LIB_PATH)
DLFLAGS += $(EXTERNAL_LIB_PATH)
DLFLAGS += $(EXTERNAL_LIBRARIES)

CPPFLAGS += $(addprefix -I, $(C_INCLUDE_DIRS) )
FPPFLAGS += $(addprefix -I, $(F_INCLUDE_DIRS) )

.SUFFIXES:	.f90	.c

main: preliminaries \
	$(LIBRARY) \
	$(MOD_INCLUDE) \
	$(MOD_FIELDML_TARGET) \
	$(HEADER_INCLUDE)

PREPROCESSED_OBJECTS = 

$(OBJECT_DIR)/%.o : $(SOURCE_DIR)/%.f90 $(OBJECT_DIR)/.directory
	( cd $(OBJECT_DIR) && $(FC) -o $@ $(FFLAGS) $(FPPFLAGS) -c $< )

$(OBJECT_DIR)/%.o : $(SOURCE_DIR)/%.c $(OBJECT_DIR)/.directory
	( cd $(OBJECT_DIR) && $(CC) -o $@ $(CFLAGS) $(CPPFLAGS) -c $< )

$(PREPROCESSED_OBJECTS) : $(OBJECT_DIR)/%.o : $(SOURCE_DIR)/%.f90 $(OBJECT_DIR)/.directory
	( m4 --prefix-builtins $< > $(subst .o,-expanded.f90,$@) && cd $(OBJECT_DIR) && $(FC) -o $@ $(FFLAGS) $(FPPFLAGS) -c $(subst .o,-expanded.f90,$@) )

# Target to create directories (but as the changing mTime of directories confuses make, we create a hidden file in it and reference it instead of the directory)
%/.directory:
	( mkdir -p $(@D) && touch $@ )

ifeq ($(USEFIELDML),true)
    FIELDML_OBJECT =  \
      $(OBJECT_DIR)/fieldml_util_routines.o \
      $(OBJECT_DIR)/fieldml_input_routines.o \
      $(OBJECT_DIR)/fieldml_output_routines.o \
      $(OBJECT_DIR)/fieldml_types.o
else
    FIELDML_OBJECT = #
endif

#ifeq ($(COMPILER),intel) # TODO: temporarily disable intel build for opencmiss.f90 and etc.
#    FIELDML_OBJECT = #
#    MOD_INCLUDE := #
#    MOD_SOURCE_INC := #
#    MOD_FIELDML_TARGET := #
#    WRAPPER_OBJECTS = #   
#else
    WRAPPER_OBJECTS =  \
    $(OBJECT_DIR)/opencmiss.o \
    $(OBJECT_DIR)/opencmiss_c.o
#endif

OBJECTS = $(OBJECT_DIR)/advection_diffusion_equation_routines.o \
	$(OBJECT_DIR)/analytic_analysis_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/bioelectric_routines.o \
	$(OBJECT_DIR)/biodomain_equation_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/blas.o \
	$(OBJECT_DIR)/Burgers_equation_routines.o \
	$(OBJECT_DIR)/classical_field_routines.o \
	$(OBJECT_DIR)/cmiss.o \
	$(OBJECT_DIR)/cmiss_c.o \
	$(OBJECT_DIR)/cmiss_cellml.o \
	$(OBJECT_DIR)/cmiss_fortran_c.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/cmiss_parmetis.o \
	$(OBJECT_DIR)/cmiss_petsc.o \
	$(OBJECT_DIR)/cmiss_petsc_types.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/Darcy_equations_routines.o \
	$(OBJECT_DIR)/Darcy_pressure_equations_routines.o \
	$(OBJECT_DIR)/finite_elasticity_Darcy_routines.o \
	$(OBJECT_DIR)/finite_elasticity_fluid_pressure_routines.o \
	$(OBJECT_DIR)/bioelectric_finite_elasticity_routines.o \
	$(OBJECT_DIR)/data_point_routines.o \
	$(OBJECT_DIR)/data_projection_routines.o \
	$(OBJECT_DIR)/diffusion_advection_diffusion_routines.o \
	$(OBJECT_DIR)/diffusion_diffusion_routines.o \
	$(OBJECT_DIR)/diffusion_equation_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/distributed_matrix_vector_IO.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/elasticity_routines.o \
	$(OBJECT_DIR)/electromechanics_routines.o \
	$(OBJECT_DIR)/electrophysiology_cell_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/equations_set_routines.o \
	$(OBJECT_DIR)/external_dae_solver_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/field_IO_routines.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/FieldExport.o \
	$(OBJECT_DIR)/fitting_routines.o \
	$(OBJECT_DIR)/generated_mesh_routines.o \
	$(OBJECT_DIR)/Hamilton_Jacobi_equations_routines.o \
	$(OBJECT_DIR)/Helmholtz_equations_routines.o \
	$(OBJECT_DIR)/history_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_routines.o \
	$(OBJECT_DIR)/interface_conditions_routines.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/interface_equations_routines.o \
	$(OBJECT_DIR)/interface_mapping_routines.o \
	$(OBJECT_DIR)/interface_matrices_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/Laplace_equations_routines.o \
	$(OBJECT_DIR)/linear_elasticity_routines.o \
	$(OBJECT_DIR)/linkedlist_routines.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/monodomain_equations_routines.o \
	$(OBJECT_DIR)/multi_compartment_transport_routines.o \
	$(OBJECT_DIR)/multi_physics_routines.o \
	$(OBJECT_DIR)/Navier_Stokes_equations_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(WRAPPER_OBJECTS) \
	$(OBJECT_DIR)/Poiseuille_equations_routines.o \
	$(OBJECT_DIR)/Poisson_equations_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/problem_routines.o \
	$(OBJECT_DIR)/reaction_diffusion_equation_routines.o \
	$(OBJECT_DIR)/reaction_diffusion_IO_routines.o \
	$(OBJECT_DIR)/region_routines.o \
	$(OBJECT_DIR)/Stokes_equations_routines.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/solver_mapping_routines.o \
	$(OBJECT_DIR)/solver_matrices_routines.o \
	$(OBJECT_DIR)/sorting.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/test_framework_routines.o \
	$(OBJECT_DIR)/timer_c.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/trees.o \
	$(OBJECT_DIR)/types.o \
	$(OBJECT_DIR)/util_array.o \
	$(FIELDML_OBJECT) \
	$(PREPROCESSED_OBJECTS)

ifeq ($(OPERATING_SYSTEM),linux)# Linux
  MACHINE_OBJECTS = $(OBJECT_DIR)/machine_constants_linux.o
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    MACHINE_OBJECTS = $(OBJECT_DIR)/machine_constants_aix.o
  else# windows
    MACHINE_OBJECTS = $(OBJECT_DIR)/machine_constants_win32.o
  endif
endif

OBJECTS += $(MACHINE_OBJECTS)

preliminaries: $(OBJECT_DIR)/.directory \
	$(INC_DIR)/.directory \
	$(LIB_DIR)/.directory

$(LIBRARY) : $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

$(MOD_INCLUDE) : $(MOD_SOURCE_INC) $(INC_DIR)/.directory
	cp $(MOD_SOURCE_INC) $@

MOD_FIELDML: $(FIELDML_OBJECT) $(INC_DIR)/.directory
	cp $(OBJECT_DIR)/fieldml_input_routines.mod $(INC_DIR)/fieldml_input_routines.mod
	cp $(OBJECT_DIR)/fieldml_output_routines.mod $(INC_DIR)/fieldml_output_routines.mod
	cp $(OBJECT_DIR)/fieldml_util_routines.mod $(INC_DIR)/fieldml_util_routines.mod
	cp $(OBJECT_DIR)/fieldml_types.mod $(INC_DIR)/fieldml_types.mod

$(HEADER_INCLUDE) $(C_F90_SOURCE): $(SOURCE_DIR)/opencmiss.f90  $(BINDINGS_GENERATE_SCRIPT)/parse.py $(BINDINGS_GENERATE_SCRIPT)/c.py
	python $(BINDINGS_GENERATE_SCRIPT) $(GLOBAL_CM_ROOT) C $(HEADER_INCLUDE) $(C_F90_SOURCE)

# Place the list of dependencies for the objects here.
#
# ----------------------------------------------------------------------------

ifeq ($(OPERATING_SYSTEM),aix)

   #Need to disable argument list checking for MPI calls which may have multiple types for the same parameters
   $(OBJECT_DIR)/computational_environment.o : DBGCF_FLGS = -qfullpath -C -qflttrap=inv:en
   $(OBJECT_DIR)/distributed_matrix_vector.o : DBGCF_FLGS = -qfullpath -C -qflttrap=inv:en
   $(OBJECT_DIR)/field_IO_routines.o : DBGCF_FLGS = -qfullpath -C -qflttrap=inv:en
   $(OBJECT_DIR)/analytic_analysis_routines.o : DBGCF_FLGS = -qfullpath -C -qflttrap=inv:en
   $(OBJECT_DIR)/data_projection_routines.o : DBGCF_FLGS = -qfullpath -C -qflttrap=inv:en

   #Need to disable argument list checking for c interface modules to allow for the c->fortran char->integer string conversion
   $(OBJECT_DIR)/timer_c.o : DBGCF_FLGS = -qfullpath -C -qflttrap=inv:en

   #Need to to auto-promote single precision constants to doubles
   $(OBJECT_DIR)/finite_elasticity_routines.o : DBGCF_FLGS = -qdpc

endif

$(OBJECT_DIR)/advection_diffusion_equation_routines.o	:	$(SOURCE_DIR)/advection_diffusion_equation_routines.f90 \
	$(OBJECT_DIR)/analytic_analysis_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/analytic_analysis_routines.o  : $(SOURCE_DIR)/analytic_analysis_routines.f90 \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/base_routines.o	:	$(SOURCE_DIR)/base_routines.f90 \
	$(MACHINE_OBJECTS) \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/basis_routines.o	:	$(SOURCE_DIR)/basis_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/blas.o	:	$(SOURCE_DIR)/blas.f90 \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/bioelectric_routines.o	:	$(SOURCE_DIR)/bioelectric_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/biodomain_equation_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/monodomain_equations_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/biodomain_equation_routines.o	:	$(SOURCE_DIR)/biodomain_equation_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_IO_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/boundary_condition_routines.o  : $(SOURCE_DIR)/boundary_condition_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/linkedlist_routines.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/classical_field_routines.o	:	$(SOURCE_DIR)/classical_field_routines.f90 \
	$(OBJECT_DIR)/Hamilton_Jacobi_equations_routines.o \
	$(OBJECT_DIR)/Helmholtz_equations_routines.o \
	$(OBJECT_DIR)/Laplace_equations_routines.o \
	$(OBJECT_DIR)/Poisson_equations_routines.o \
	$(OBJECT_DIR)/advection_diffusion_equation_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/diffusion_equation_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/reaction_diffusion_equation_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/cmiss.o	:	$(SOURCE_DIR)/cmiss.f90 \
	$(MACHINE_OBJECTS) \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/generated_mesh_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_routines.o \
	$(OBJECT_DIR)/region_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/cmiss_c.o	:	$(SOURCE_DIR)/cmiss_c.c

$(OBJECT_DIR)/cmiss_cellml.o	:	$(SOURCE_DIR)/cmiss_cellml.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/cmiss_fortran_c.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/cmiss_cellml_dummy.o	:	$(SOURCE_DIR)/cmiss_cellml_dummy.f90

$(OBJECT_DIR)/cmiss_fortran_c.o	:	$(SOURCE_DIR)/cmiss_fortran_c.f90

$(OBJECT_DIR)/cmiss_mpi.o	:	$(SOURCE_DIR)/cmiss_mpi.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o

$(OBJECT_DIR)/cmiss_parmetis.o	:	$(SOURCE_DIR)/cmiss_parmetis.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/cmiss_petsc.o	:	$(SOURCE_DIR)/cmiss_petsc.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/cmiss_petsc_types.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/cmiss_petsc_types.o	:	$(SOURCE_DIR)/cmiss_petsc_types.f90 \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/computational_environment.o	:	$(SOURCE_DIR)/computational_environment.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/cmiss_petsc.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/constants.o	:	$(SOURCE_DIR)/constants.f90 \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/control_loop_routines.o	:	$(SOURCE_DIR)/control_loop_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/coordinate_routines.o	:	$(SOURCE_DIR)/coordinate_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/Darcy_equations_routines.o	:	$(SOURCE_DIR)/Darcy_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/Darcy_pressure_equations_routines.o	:	$(SOURCE_DIR)/Darcy_pressure_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/fieldml_input_routines.o: $(SOURCE_DIR)/fieldml_input_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/cmiss.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fieldml_types.o \
	$(OBJECT_DIR)/fieldml_util_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/region_routines.o \
	$(OBJECT_DIR)/strings.o

$(OBJECT_DIR)/fieldml_output_routines.o: $(SOURCE_DIR)/fieldml_output_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fieldml_types.o \
	$(OBJECT_DIR)/fieldml_util_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/region_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/fieldml_util_routines.o: $(SOURCE_DIR)/fieldml_util_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fieldml_types.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/region_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/fieldml_types.o: $(SOURCE_DIR)/fieldml_types.f90 \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/finite_elasticity_Darcy_routines.o	:	$(SOURCE_DIR)/finite_elasticity_Darcy_routines.f90 \
	$(OBJECT_DIR)/Darcy_equations_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/finite_elasticity_fluid_pressure_routines.o	:	$(SOURCE_DIR)/finite_elasticity_fluid_pressure_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/bioelectric_finite_elasticity_routines.o	:	$(SOURCE_DIR)/bioelectric_finite_elasticity_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/biodomain_equation_routines.o \
	$(OBJECT_DIR)/bioelectric_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_IO_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/data_point_routines.o	:	$(SOURCE_DIR)/data_point_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/data_projection_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/trees.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/data_projection_routines.o	:	$(SOURCE_DIR)/data_projection_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/sorting.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/trees.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/diffusion_advection_diffusion_routines.o	:	$(SOURCE_DIR)/diffusion_advection_diffusion_routines.f90 \
	$(OBJECT_DIR)/advection_diffusion_equation_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/diffusion_equation_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/diffusion_diffusion_routines.o	:	$(SOURCE_DIR)/diffusion_diffusion_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/diffusion_equation_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/diffusion_equation_routines.o	:	$(SOURCE_DIR)/diffusion_equation_routines.f90 \
	$(OBJECT_DIR)/analytic_analysis_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/Burgers_equation_routines.o	:	$(SOURCE_DIR)/Burgers_equation_routines.f90 \
	$(OBJECT_DIR)/analytic_analysis_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/distributed_matrix_vector.o	:	$(SOURCE_DIR)/distributed_matrix_vector.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/cmiss_petsc.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/linkedlist_routines.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/distributed_matrix_vector_IO.o	:	$(SOURCE_DIR)/distributed_matrix_vector_IO.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/domain_mappings.o	:	$(SOURCE_DIR)/domain_mappings.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/elasticity_routines.o	:	$(SOURCE_DIR)/elasticity_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/linear_elasticity_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/electromechanics_routines.o	:	$(SOURCE_DIR)/electromechanics_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/electrophysiology_cell_routines.o	:	$(SOURCE_DIR)/electrophysiology_cell_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/equations_routines.o	:	$(SOURCE_DIR)/equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/equations_mapping_routines.o	:	$(SOURCE_DIR)/equations_mapping_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/equations_matrices_routines.o	:	$(SOURCE_DIR)/equations_matrices_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/linkedlist_routines.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/equations_set_constants.o	:	$(SOURCE_DIR)/equations_set_constants.f90 \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/equations_set_routines.o	:	$(SOURCE_DIR)/equations_set_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/bioelectric_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/classical_field_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/elasticity_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fitting_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/monodomain_equations_routines.o \
	$(OBJECT_DIR)/multi_physics_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/external_dae_solver_routines.o	:	$(SOURCE_DIR)/external_dae_solver_routines.c \
	$(SOURCE_DIR)/external_dae_solver_routines.h

$(OBJECT_DIR)/field_routines.o	:	$(SOURCE_DIR)/field_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/field_IO_routines.o	:	$(SOURCE_DIR)/field_IO_routines.f90 \
	$(MACHINE_OBJECTS) \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o \
	$(SOURCE_DIR)/FieldExportConstants.h

$(OBJECT_DIR)/finite_elasticity_routines.o	:	$(SOURCE_DIR)/finite_elasticity_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/cmiss_petsc.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/generated_mesh_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/fluid_mechanics_routines.o	:	$(SOURCE_DIR)/fluid_mechanics_routines.f90 \
	$(OBJECT_DIR)/Burgers_equation_routines.o \
	$(OBJECT_DIR)/Darcy_equations_routines.o \
	$(OBJECT_DIR)/Darcy_pressure_equations_routines.o \
	$(OBJECT_DIR)/Navier_Stokes_equations_routines.o \
	$(OBJECT_DIR)/Poiseuille_equations_routines.o \
	$(OBJECT_DIR)/Stokes_equations_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/fluid_mechanics_IO_routines.o	:	$(SOURCE_DIR)/fluid_mechanics_IO_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/FieldExport.o	:	$(SOURCE_DIR)/FieldExport.c \
	$(SOURCE_DIR)/FieldExportConstants.h

$(OBJECT_DIR)/input_output.o	:	$(SOURCE_DIR)/input_output.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o

$(OBJECT_DIR)/fitting_routines.o	:	$(SOURCE_DIR)/fitting_routines.f90 \
	$(OBJECT_DIR)/Darcy_equations_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/generated_mesh_routines.o	:	$(SOURCE_DIR)/generated_mesh_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/Helmholtz_equations_routines.o	:	$(SOURCE_DIR)/Helmholtz_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/history_routines.o	:	$(SOURCE_DIR)/history_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/interface_routines.o	:	$(SOURCE_DIR)/interface_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/generated_mesh_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_conditions_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/interface_conditions_constants.o	:	$(SOURCE_DIR)/interface_conditions_constants.f90 \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/interface_conditions_routines.o	:	$(SOURCE_DIR)/interface_conditions_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/interface_equations_routines.o \
	$(OBJECT_DIR)/interface_mapping_routines.o \
	$(OBJECT_DIR)/interface_matrices_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/interface_equations_routines.o	:	$(SOURCE_DIR)/interface_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/interface_mapping_routines.o \
	$(OBJECT_DIR)/interface_matrices_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/interface_mapping_routines.o	:	$(SOURCE_DIR)/interface_mapping_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/interface_matrices_routines.o	:	$(SOURCE_DIR)/interface_matrices_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/iso_varying_string.o	:	$(SOURCE_DIR)/iso_varying_string.f90

$(OBJECT_DIR)/kinds.o	:	$(SOURCE_DIR)/kinds.f90

$(OBJECT_DIR)/Hamilton_Jacobi_equations_routines.o	:	$(SOURCE_DIR)/Hamilton_Jacobi_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/Laplace_equations_routines.o	:	$(SOURCE_DIR)/Laplace_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/linear_elasticity_routines.o	:	$(SOURCE_DIR)/linear_elasticity_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/lists.o	:	$(SOURCE_DIR)/lists.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o


$(OBJECT_DIR)/linkedlist_routines.o	:	$(SOURCE_DIR)/linkedlist_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/machine_constants_aix.o	:	$(SOURCE_DIR)/machine_constants_aix.f90 \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/machine_constants_linux.o	:	$(SOURCE_DIR)/machine_constants_linux.f90 \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/machine_constants_win32.o	:	$(SOURCE_DIR)/machine_constants_win32.f90 \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/maths.o	:	$(SOURCE_DIR)/maths.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/matrix_vector.o	:	$(SOURCE_DIR)/matrix_vector.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/linkedlist_routines.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/mesh_routines.o	:	$(SOURCE_DIR)/mesh_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/cmiss_mpi.o \
	$(OBJECT_DIR)/cmiss_parmetis.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/trees.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/monodomain_equations_routines.o    :       $(SOURCE_DIR)/monodomain_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/electrophysiology_cell_routines.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fitting_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/multi_compartment_transport_routines.o	:	$(SOURCE_DIR)/multi_compartment_transport_routines.f90 \
	$(OBJECT_DIR)/advection_diffusion_equation_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/diffusion_equation_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/maths.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/multi_physics_routines.o	:	$(SOURCE_DIR)/multi_physics_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/bioelectric_finite_elasticity_routines.o \
	$(OBJECT_DIR)/diffusion_advection_diffusion_routines.o \
	$(OBJECT_DIR)/diffusion_diffusion_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/finite_elasticity_Darcy_routines.o \
	$(OBJECT_DIR)/finite_elasticity_fluid_pressure_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/multi_compartment_transport_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/Navier_Stokes_equations_routines.o	:	$(SOURCE_DIR)/Navier_Stokes_equations_routines.f90 \
	$(OBJECT_DIR)/analytic_analysis_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/node_routines.o	:	$(SOURCE_DIR)/node_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/trees.o \
	$(OBJECT_DIR)/types.o

$(MOD_SOURCE_INC) : $(OBJECT_DIR)/opencmiss.o

$(OBJECT_DIR)/opencmiss.o	:	$(SOURCE_DIR)/opencmiss.f90 \
	$(FIELDML_OBJECT) \
	$(OBJECT_DIR)/Hamilton_Jacobi_equations_routines.o \
	$(OBJECT_DIR)/analytic_analysis_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/cmiss.o \
	$(OBJECT_DIR)/cmiss_cellml.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/data_point_routines.o \
	$(OBJECT_DIR)/data_projection_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/equations_set_routines.o \
	$(OBJECT_DIR)/field_IO_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/generated_mesh_routines.o \
	$(OBJECT_DIR)/history_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/interface_conditions_routines.o \
	$(OBJECT_DIR)/interface_equations_routines.o \
	$(OBJECT_DIR)/interface_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/problem_routines.o \
	$(OBJECT_DIR)/region_routines.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/opencmiss_c.o	:	$(C_F90_SOURCE) \
	$(OBJECT_DIR)/cmiss_fortran_c.o \
	$(OBJECT_DIR)/opencmiss.o

$(OBJECT_DIR)/Poiseuille_equations_routines.o	:	$(SOURCE_DIR)/Poiseuille_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o


$(OBJECT_DIR)/Poisson_equations_routines.o	:	$(SOURCE_DIR)/Poisson_equations_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/problem_constants.o	:	$(SOURCE_DIR)/problem_constants.f90 \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/problem_routines.o	:	$(SOURCE_DIR)/problem_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/bioelectric_routines.o \
	$(OBJECT_DIR)/classical_field_routines.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/elasticity_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/equations_set_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/finite_elasticity_routines.o \
	$(OBJECT_DIR)/fitting_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_conditions_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/multi_physics_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_matrices_routines.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o



$(OBJECT_DIR)/reaction_diffusion_equation_routines.o	:	$(SOURCE_DIR)/reaction_diffusion_equation_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/reaction_diffusion_IO_routines.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/reaction_diffusion_IO_routines.o	:	$(SOURCE_DIR)/reaction_diffusion_IO_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/region_routines.o	:	$(SOURCE_DIR)/region_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/cmiss_cellml.o \
	$(OBJECT_DIR)/coordinate_routines.o \
	$(OBJECT_DIR)/data_point_routines.o \
	$(OBJECT_DIR)/equations_set_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/generated_mesh_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_routines.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/mesh_routines.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/solver_routines.o	:	$(SOURCE_DIR)/solver_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/cmiss_cellml.o \
	$(OBJECT_DIR)/cmiss_petsc.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/external_dae_solver_routines.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_mapping_routines.o \
	$(OBJECT_DIR)/solver_matrices_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/solver_mapping_routines.o	:	$(SOURCE_DIR)/solver_mapping_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/computational_environment.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/lists.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/solver_matrices_routines.o	:	$(SOURCE_DIR)/solver_matrices_routines.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/interface_conditions_constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/sorting.o	:	$(SOURCE_DIR)/sorting.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/Stokes_equations_routines.o	:	$(SOURCE_DIR)/Stokes_equations_routines.f90 \
	$(OBJECT_DIR)/analytic_analysis_routines.o \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/basis_routines.o \
	$(OBJECT_DIR)/boundary_condition_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/control_loop_routines.o \
	$(OBJECT_DIR)/distributed_matrix_vector.o \
	$(OBJECT_DIR)/domain_mappings.o \
	$(OBJECT_DIR)/equations_mapping_routines.o \
	$(OBJECT_DIR)/equations_matrices_routines.o \
	$(OBJECT_DIR)/equations_routines.o \
	$(OBJECT_DIR)/equations_set_constants.o \
	$(OBJECT_DIR)/field_routines.o \
	$(OBJECT_DIR)/fluid_mechanics_IO_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/node_routines.o \
	$(OBJECT_DIR)/problem_constants.o \
	$(OBJECT_DIR)/solver_routines.o \
	$(OBJECT_DIR)/strings.o \
	$(OBJECT_DIR)/timer_f.o \
	$(OBJECT_DIR)/types.o

$(OBJECT_DIR)/strings.o	:	$(SOURCE_DIR)/strings.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o

$(OBJECT_DIR)/test_framework_routines.o	:	$(SOURCE_DIR)/test_framework_routines.f90 \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/matrix_vector.o \
	$(OBJECT_DIR)/strings.o

$(OBJECT_DIR)/timer_c.o		:	$(SOURCE_DIR)/timer_c.c

$(OBJECT_DIR)/timer_f.o	:	$(SOURCE_DIR)/timer_f.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o

$(OBJECT_DIR)/trees.o	:	$(SOURCE_DIR)/trees.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/input_output.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/strings.o

$(OBJECT_DIR)/types.o	:	$(SOURCE_DIR)/types.f90 \
	$(OBJECT_DIR)/cmiss_petsc_types.o \
	$(OBJECT_DIR)/constants.o \
	$(OBJECT_DIR)/iso_varying_string.o \
	$(OBJECT_DIR)/kinds.o \
	$(OBJECT_DIR)/linkedlist_routines.o \
	$(OBJECT_DIR)/trees.o

$(OBJECT_DIR)/util_array.o   :       $(SOURCE_DIR)/util_array.f90 \
	$(OBJECT_DIR)/base_routines.o \
	$(OBJECT_DIR)/types.o

# ----------------------------------------------------------------------------
#
# SWIG bindings to other languages

GENERATED_INTERFACE = $(BINDINGS_DIR)/python/opencmiss.i
PYTHON_INTERFACE = $(BINDINGS_DIR)/python/opencmiss_py.i
PYTHON_MODULE = $(BINDINGS_DIR)/python/opencmiss/CMISS.py
PYTHON_MODULE_SO = $(BINDINGS_DIR)/python/opencmiss/_opencmiss_swig.so
PYTHON_WRAPPER = $(BINDINGS_DIR)/python/opencmiss/opencmiss_wrap.c
PYTHON_INCLUDES = $(shell python-config --includes)

python: $(PYTHON_MODULE) $(PYTHON_MODULE_SO)

$(GENERATED_INTERFACE): $(BINDINGS_GENERATE_SCRIPT)/parse.py $(BINDINGS_GENERATE_SCRIPT)/swig.py $(SOURCE_DIR)/opencmiss.f90
	python $(BINDINGS_GENERATE_SCRIPT) $(GLOBAL_CM_ROOT) SWIG $@

$(PYTHON_MODULE): $(BINDINGS_GENERATE_SCRIPT)/parse.py $(BINDINGS_GENERATE_SCRIPT)/python.py $(SOURCE_DIR)/opencmiss.f90
	python $(BINDINGS_GENERATE_SCRIPT) $(GLOBAL_CM_ROOT) Python

$(PYTHON_WRAPPER) : $(PYTHON_INTERFACE) $(GENERATED_INTERFACE) $(HEADER_INCLUDE)
# Remove opencmiss_swig.py after running SWIG as we generate our own Python wrapper code
	( cd $(BINDINGS_DIR)/python/opencmiss && swig -python -o $@ -module opencmiss_swig -outdir . -I$(INC_DIR) $(PYTHON_INTERFACE) && rm opencmiss_swig.py )

$(PYTHON_MODULE_SO) : $(LIBRARY) $(PYTHON_WRAPPER) $(OBJECTS)
	( cd $(BINDINGS_DIR)/python && $(CC) -c $(PYTHON_WRAPPER) $(CFLAGS) $(CPPFLAGS) -I$(INC_DIR) $(PYTHON_INCLUDES) -o opencmiss/opencmiss_wrap.o )
	( cd $(BINDINGS_DIR)/python && $(CC) opencmiss/opencmiss_wrap.o $(OBJECTS) $(DLFLAGS) -o $(PYTHON_MODULE_SO) )

# ----------------------------------------------------------------------------
#
# clean and clobber for removing objects and executable.

clean:
	@echo "Cleaning house ..."
	rm -rf $(OBJECT_DIR) $(LIBRARY) $(MOD_INCLUDE) $(HEADER_INCLUDE) $(C_F90_SOURCE)

allclean:
	@echo "Cleaning house ..."
	rm -rf object/* lib/*

clobber: clean
	rm -f $(LIBRARY)

externallibs:
	$(MAKE) --no-print-directory -f $(EXTERNAL_CM_ROOT)/packages/Makefile DEBUG=$(DEBUG) ABI=$(ABI)

debug opt debug64 opt64:
	$(MAKE) --no-print-directory DEBUG=$(DEBUG) ABI=$(ABI)

debug debug64: DEBUG=true
opt opt64: DEBUG=false
ifneq (,$(filter $(MACHNAME),ia64 x86_64))# ia64 or x86_64
   debug opt: ABI=64
else
   debug opt: ABI=32
endif
debug64 opt64: ABI=64

all: debug opt
all64: debug64 opt64

ifndef SIZE
  SIZE=small
endif

ifndef DIR
  DIR=.
endif

test: main
	@echo "================================================================================================"
	@echo "By default, this will test small set of examples. (All tests in the nightly build)"
	@echo "To test large set of examples (All tests in the nightly build), please set SIZE=large"
	@echo "To test examples inside a given directory, please set DIR=<DIR>, e.g. DIR=ClassicalField/Laplace"
	@echo "You need to install nosetests in order to run the tests."
	@echo "Please go to http://readthedocs.org/docs/nose/ for details."
	@echo "For detailed logfiles, go to <OPENCMISS_ROOT>/build/logs directory."
	@echo "================================================================================================"
	COMPILER=$(COMPILER) SIZE=${SIZE} DIR=$(DIR) ABI=${ABI} nosetests ${OPENCMISSEXAMPLES_ROOT}/noseMain.py:test_example

#-----------------------------------------------------------------------------

help:
	@echo "           Compile a library version of OpenCMISS"
	@echo "           ======================================"
	@echo
	@echo "Examples of usage:   "
	@echo
	@echo "	gmake"
	@echo "	gmake OPT= ABI=32"
	@echo "	gmake PROF="
	@echo "	gmake debug64"
	@echo
	@echo "Options: (The former is the default unless specified.)"
	@echo
	@echo "	(DEBUG=|OPT=)"
	@echo "	MPI=(mpich2|intel|openmpi|mvapich2|cray)"
	@echo "	PROF=(true|)"
	@echo "	MPIPROF=(true|)"
	@echo "	ABI=(32|64)"
	@echo "	COMPILER=(intel|gnu|ibm|cray)"
	@echo "	USECELLML=(false|true)"
	@echo "	USEFIELDML=(false|true)"
	@echo
	@echo "Available targets:                            "
	@echo
	@echo "	clean"
	@echo "		Remove generated files associated with a single"
	@echo "		version."
	@echo
	@echo "	clobber"
	@echo "		Remove all files associated with a single version."
	@echo
	@echo "	help"
	@echo "		Display this message."
	@echo
	@echo "	debug opt debug64 opt64"
	@echo "		Compile the specified version with automatic setting"
	@echo "		of DEBUG, ABI, and MP."
	@echo
	@echo "	all"
	@echo "		Compile all versions."
	@echo
	@echo "	all64"
	@echo "		Compile all 64-bit versions."
	@echo
	@echo "	python"
	@echo "		Build the Python bindings. To install them, change directory to bindings/python and run 'python setup.py install'."
	@echo
	@echo "	externallibs"
	@echo "		Compile the external libraries."
	@echo
	@echo "	test"
	@echo "		Build, execute and check for test."
	@echo

