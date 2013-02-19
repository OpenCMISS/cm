# -*- makefile -*-
#
# For use with GNU make.
#
#
#----------------------------------------------------------------------------------------------------------------------------------
# Define variables used by both Makefile and ExampleMakefile
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

ifndef OPENCMISS_ROOT
  OPENCMISS_ROOT := ../.
  OC_GLOBAL_ROOT := $(OPENCMISS_ROOT)
else
  OC_GLOBAL_ROOT := ${OPENCMISS_ROOT}
endif
ifndef OC_CM_GLOBAL_ROOT
  OC_CM_GLOBAL_ROOT := $(OCE_GLOBAL_ROOT)/cm
endif
ifndef OC_CELLML_GLOBAL_ROOT
  OC_CELLML_GLOBAL_ROOT := $(OCE_GLOBAL_ROOT)/cellml
endif

ifndef OPENCMISSEXTRAS_ROOT
  OPENCMISSEXTRAS_ROOT := ../../OpenCMISSExtras
  OCE_GLOBAL_ROOT := $(OPENCMISSEXTRAS_ROOT)
else
  OCE_GLOBAL_ROOT := ${OPENCMISSEXTRAS_ROOT}
endif

OCE_MAKEINC_ROOT := $(OCE_GLOBAL_ROOT)/utils/make

include $(OCE_MAKEINC_ROOT)/Makefile.inc
include $(OCE_MAKEINC_ROOT)/Makefile_Compilers.inc
include $(OCE_MAKEINC_ROOT)/Makefile_MPI.inc
include $(OCE_MAKEINC_ROOT)/Makefile_BlasLapack.inc
include $(OCE_MAKEINC_ROOT)/Makefile_BlacsScalapack.inc
include $(OCE_MAKEINC_ROOT)/Makefile_HDF5.inc
include $(OCE_MAKEINC_ROOT)/Makefile_NetCDF.inc
include $(OCE_MAKEINC_ROOT)/Makefile_PAPI.inc
include $(OCE_MAKEINC_ROOT)/Makefile_CUDA.inc
include $(OCE_MAKEINC_ROOT)/Makefile_OSOperations.inc

ifndef OCE_CM_INSTALL_DIR
  OCE_CM_INSTALL_DIR := $(CM_INSTALL_DIR)
endif
ifndef OCE_CELLML_INSTALL_DIR
  OCE_CELLML_INSTALL_DIR := $(CELLML_INSTALL_DIR)
endif
ifndef OCE_FIELDML_INSTALL_DIR
  OCE_FIELDML_INSTALL_DIR := $(FIELDML_INSTALL_DIR)
endif
ifndef OCE_COMMON_INSTALL_DIR
  OCE_COMMON_INSTALL_DIR := $(COMMON_INSTALL_DIR)
endif
ifndef OCE_MPI_INSTALL_DIR
  OCE_MPI_INSTALL_DIR := $(MPI_INSTALL_DIR)
endif
ifndef OC_CELLML_ROOT
  OC_CELLML_ROOT := ${OC_GLOBAL_ROOT}/cellml
endif
OC_CELLML_INSTALL_DIR := $(OC_CELLML_ROOT)/$(ENVIRONMENT)$(BUILD_TYPE)$(NOMPI_TOOLCHAIN)

#----------------------------------------------------------------------------------------------------------------------------------

MAIN_OBJECT_DIR := $(OC_CM_GLOBAL_ROOT)/object/$(ENVIRONMENT)$(BUILD_TYPE)$(TOOLCHAIN)
INC_DIR := $(OC_CM_GLOBAL_ROOT)/include/$(ENVIRONMENT)$(BUILD_TYPE)$(TOOLCHAIN)
LIB_DIR := $(OC_CM_GLOBAL_ROOT)/lib/$(ENVIRONMENT)$(BUILD_TYPE)$(TOOLCHAIN)

#----------------------------------------------------------------------------------------------------------------------------------
# compiling commands

MPIFC = $(OCE_MPIF90)
MPICC = $(OCE_MPICC)
FC = $(MPIFC)
CC = $(MPICC)
AR = ar
EXE_LINK = $(FC)
DSO_LINK = ld

DBGCF_FLGS = -g #OPT=false flags for C and fortran
# Option lists
# (suboption lists become more specific so that later ones overrule previous)
CFLAGS = $(strip $(CFL_FLGS) $(CFE_FLGS) $(CF_FLGS) $(C_FLGS))
FFLAGS = $(strip $(CFL_FLGS) $(CFE_FLGS) $(CF_FLGS) $(F_FLGS))
ELFLAGS = $(strip $(CFL_FLGS) $(L_FLGS) $(CFE_FLGS))
DLFLAGS = $(strip $(CFL_FLGS) $(L_FLGS) $(D_FLGS))
ifneq ($(DEBUG),false)
  CFLAGS += $(strip $(DBGCF_FLGS) $(DBGC_FLGS))
  FFLAGS += $(strip $(DBGCF_FLGS) $(DBGF_FLGS))
  DLFLAGS += $(strip $(DBGCF_FLGS) $(DBGC_FLGS))
  CPPFLAGS += -DDEBUG
else
  ifneq ($(MPIPROF),false)
    CFLAGS += $(strip $(DBGCF_FLGS) $(DBGC_FLGS))
    FFLAGS += $(strip $(DBGCF_FLGS) $(DBGF_FLGS))
    CPPFLAGS += -DDEBUG
  else
    CFLAGS += $(strip $(OPTCFE_FLGS) $(OPTCF_FLGS) $(OPTC_FLGS))
    FFLAGS += $(strip $(OPTCFE_FLGS) $(OPTCF_FLGS) $(OPTF_FLGS))
    ELFLAGS += $(OPTCFE_FLGS)
  endif
endif
ifneq ($(MP),false)
  CFLAGS += $(MP_FLGS)
  FFLAGS += $(MP_FLGS)
endif

ARFLAGS = -crsv

# suboption lists
CFL_FLGS =#    flags for C fortran and linking
L_FLGS =#      flags for linking only
CFE_FLGS =#    flags for C fortran and linking executables only
CF_FLGS = -c#  flags for C and fortran only
C_FLGS =#      flags for C only
F_FLGS =#      flags for fortran only
D_FLGS = -shared#     for linking dynamic shared objects only
DBGC_FLGS =#   OPT=false flags for C only
DBGF_FLGS =#   OPT=false flags for fortran only
OPTCFE_FLGS =# OPT=true flags for C and fortran and linking executables
OPTCF_FLGS = -O#OPT=true flags for C and fortran only
OPTC_FLGS =#   OPT=true flags for C only
OPTF_FLGS =#   OPT=true flags for fortran only

# The list of objects may be too long for the operating system limit on
# argument length so the list of objects is stored in a file.  This linker
# arguments for this file depend on the linker.  If the linker cannot
# interpret such a file then try to use the shell and hope the list isn't too
# long.
olist_args = `cat $1`

CC = $(OCE_MPICC)
FC = $(OCE_MPIF90)

#----------------------------------------------------------------------------------------------------------------------------------
ifeq ($(OPERATING_SYSTEM),linux)
  OPTCF_FLGS =# Use separate flags for fortran and c
  olist_args = $1
  # Set the flags for the various different CC compilers
  ifeq ($(COMPILER),gnu)
    # Position independent code is required for the Python bindings
    C_FLGS += -pipe -m$(ABI) -fPIC
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        C_FLGS += -march=nocona
        F_FLGS += -march=nocona
      endif
    endif
    DBGC_FLGS += -O0 -fbounds-check
    OPTC_FLGS = -O3 -funroll-all-loops
    ifeq ($(PROF),false)
      ifneq ($(filter $(INSTRUCTION),i686 x86_64),)# i686 or x86_64
        OPTC_FLGS += -momit-leaf-frame-pointer
        OPTF_FLGS += -momit-leaf-frame-pointer
      endif
    else
      C_FLGS += -g -pg -fprofile-arcs -ftest-coverage
      ELFLGS += -g -pg -fprofile-arcs -ftest-coverage
      F_FLGS += -g -pg -fprofile-arcs -ftest-coverage
      ELFLAGS += -pg -fprofile-arcs -ftest-coverage
    endif
    F_FLGS += -pipe -m$(ABI) -fno-second-underscore -Wall -x f95-cpp-input -fPIC
    # Restrict line length to 132
    F_FLGS += -ffree-line-length-132
    # for now change max identifier length. Should restrict to 63 (F2003) in future
    F_FLGS += -fmax-identifier-length=63
    DBGF_FLGS += -O0 -ffpe-trap=invalid,zero -fbacktrace
    ifdef COMPILER_VERSION
      ifeq ($(COMPILER_VERSION),4.4)
        DBGF_FLGS += -fbounds-check
      else
        # For any more recent versions
        DBGF_FLGS += -fcheck=all
      endif
    endif
    OPTF_FLGS = -O3 -Wuninitialized -funroll-all-loops
    ELFLAGS += -m$(ABI)
  else
    ifeq ($(COMPILER),intel)
      # Turn on all warnings
      C_FLGS += -Wall -m$(ABI) -fPIC
      ifeq ($(MACHNAME),x86_64)
        ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
          ifneq ($(shell grep Duo /proc/cpuinfo 2>/dev/null),)
            ifneq ($(shell grep "Core(TM)2" /proc/cpuinfo 2>/dev/null),)
              C_FLGS += -xSSSE3# Core2 Duo
              F_FLGS += -xSSSE3# Core2 Duo
            else
              C_FLGS += -msse3# Core Duo
              F_FLGS += -msse3# Core Duo
            endif
          else
            C_FLGS += -xSSE3# for sse3 (90nm Pentium 4 series)
            F_FLGS += -xSSE3# for sse3 (90nm Pentium 4 series)
          endif
        else
          C_FLGS += -msse2# Pentium4 compatible (?sse2)
          F_FLGS += -msse2# Pentium4 compatible (?sse2)
        endif
      endif
      ifeq ($(filter-out i%86,$(MACHNAME)),)
        ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
          C_FLGS += -xSSE2# for Pentium 4
          F_FLGS += -xSSE2# for Pentium 4
        endif
      endif
      DBGC_FLGS += -O0
      OPTC_FLGS = -O3 -ansi_alias
      ifneq ($(PROF),false)
        C_FLGS += -g -pg
        F_FLGS += -g -pg
        ELFLAGS += -pg
      endif
      ifeq ($(MPIPROF),true)
        ifeq ($(MPI),mpich2)
          C_FLGS += -Wl,--export-dynamic
          DBGC_FLGS += -Wl,--export-dynamic
          F_FLAS += -Wl,--export-dynamic
          DBGF_FLGS += -Wl,--export-dynamic
        else
          F_FLGS += -tcollect
          C_FLGS += -tcollect
        endif
      endif
      # turn on preprocessing,
      # turn on warnings,
      # warn about non-standard Fortran 95
      F_FLGS += -cpp -warn all -m$(ABI) -fPIC
      DBGF_FLGS += -O0 -check all -traceback -debug all
      OPTF_FLGS = -O3
#     MP_FLGS = -openmp
      ELFLAGS += -nofor_main -m$(ABI) -traceback
    else
      ifeq ($(COMPILER),ibm)
        CFLAGS += -qinfo=gen:ini:por:pro:trd:tru:use
        C_FLGS += -q$(ABI) -qarch=auto -qhalt=e
        # -qinitauto for C is bytewise: 7F gives large integers.
        DBGC_FLGS += -qfullpath -C -qflttrap=inv:en -qinitauto=7F
        OPTC_FLGS += -O3
        # for trailing _ on fortran symbols
        CPPFLAGS += -Dunix
        F_FLGS += -q$(ABI) -qarch=auto -qhalt=e -qextname -qsuffix=cpp=f90
        ELFLAGS += -q$(ABI)
        ifeq ($(ABI),64)
          F_FLGS += -qwarn64
        endif
        ifeq ($(DEBUG),false)
          MP_FLGS = -qsmp=omp
        else
          MP_FLGS = -qsmp=omp:noopt
        endif
        # -qinitauto for Fortran 7FF7FFFF is a large integer or NaNQ real*4 or NaNS real*8
        DBGF_FLGS += -qfullpath -C -qflttrap=inv:en -qextchk -qinitauto=7FF7FFFF
        OPTF_FLGS += -O3
      else
        ifeq ($(COMPILER),cray)
          DBGC_FLGS += -O0 -g
          OPTC_FLGS = -O3
          ifeq ($(PROF),true)
            C_FLGS += -g -pg
            F_FLGS += -g -pg
            ELFLAGS += -pg
          endif
          DBGF_FLGS += -O0 -g
          OPTF_FLGS = -O3
	else
	  $(error COMPILER is not valid
	endif
      endif
    endif
  endif

  # Avoid versioning problems with libgcc_s by linking statically.

  # libgcc2.c from gcc 3.4.4 says:
  # In addition to the permissions in the GNU General Public License, the
  # Free Software Foundation gives you unlimited permission to link the
  # compiled version of this file into combinations with other programs,
  # and to distribute those combinations without any restriction coming
  # from the use of this file.

  # (With dynamic version, should copy libgcc_s.so.N if copying libstdc++.so.N)
  ELFLAGS += -static-libgcc

  # Use the BSD timers
  CPPFLAGS += -DBSD_TIMERS
endif
ifeq ($(OPERATING_SYSTEM),win32)
  F_FLGS += -fno-second-underscore
  OPTCF_FLGS = -O2
  ELFLAGS += -Wl,-static
  # Use the ANSI C timers
  CPPFLAGS += -DANSI_C_TIMERS
  olist_args = $1
endif
ifeq ($(OPERATING_SYSTEM),aix)
  F_FLGS += -qsuffix=cpp=f90 -qnoextname
  CFLAGS += -qinfo=gen:ini:por:pro:trd:tru:use
  ELFLAGS += -q$(ABI)
  CFE_FLGS += -q$(ABI) -qarch=auto -qhalt=e
  L_FLGS += -b$(ABI)
  D_FLGS = -G -bexpall -bnoentry
  ifeq ($(ABI),32)
    # Without -bmaxdata, the only one 256M virtual segment is available for
    # data.
    # In AIX 5.3, 0xAFFFFFFF is the largest value we can use here and still
    # use global shared libraries. (see aixprggd/genprogc/lrg_prg_support.htm)
    # However, 0xAFFFFFFF/dsa causes the system to crash on loading of perl
    # modules (File::Find and Time::HiRes atleast).  0x80000000 seems to work.
    # dsa allows segments to be allocated dynamically for shmat/mmap or data
    # as required.
    ELFLAGS += -bmaxdata:0x80000000/dsa
  else
    CF_FLGS += -qwarn64
    # It seems that somewhere between AIX 5.1 and 5.3 the kernel loader
    # started modifying a process's soft data resource limit to match to match
    # its maxdata value (if non-zero).  As 32-bit applications need a non-zero
    # maxdata value to access more than 256M of data many applications
    # (including perl) will cause the data limit to be lowered to a 32-bit
    # addressable value.  As cmiss is likely to be a child of such 32-bit
    # processes, to access more than 32-bit addressable memory, it either
    # needs to raise its data limit or use its own maxdata value.
    # max heap size is 0x06FFFFFFFFFFFFF8
    # Arbitrary.  0x0000100000000000 should provide ~16TB.
    ELFLAGS += -bmaxdata:0x0000100000000000
  endif
  ifeq ($(DEBUG),false)
    MP_FLGS = -qsmp=omp
  else
    MP_FLGS = -qsmp=omp:noopt
  endif
  # Should -qflttrap=nans be used as well or instead of -qflttrap=inv:en?
  DBGCF_FLGS += -qfullpath -C -qflttrap=inv:en -qextchk
  # -qinitauto for Fortran: 7FF7FFFF is a large integer or NaNQ real*4 or NaNS real*8
  # -qinitauto for C is bytewise: 7F gives large integers.
  DBGF_FLGS += -qinitauto=7FF7FFFF
  DBGC_FLGS += -qinitauto=7F
  OPTCF_FLGS = -O3
  OPTC_FLGS += -qnoignerrno
  olist_args = -f $1
  # Use the BSD timers
  CPPFLAGS += -DBSD_TIMERS
endif

# Check that call function works (for searchdirs, olist_args, etc.)
ifeq ($(call olist_args,test),)
  $(error call function not available.  Use GNU make version 3.78 or newer)
endif

#SLEPc
SLEPC_INCLUDE_PATH =#
SLEPC_LIBRARIES =#
SLEPC_LIB_PATH =#
ifeq (${OPENCMISS_USE_SLEPC},true)
  SLEPC_LIBRARIES += -lslepc
  SLEPC_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#TAO
TAO_INCLUDE_PATH =#
TAO_LIBRARIES =#
TAO_LIB_PATH =#
ifeq (${OPENCMISS_USE_TAO},true)
  TAO_LIBRARIES += -ltao
  TAO_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#PETSc
PETSC_INCLUDE_PATH =#
PETSC_LIBRARIES =#
PETSC_LIB_PATH =#
ifeq (${OPENCMISS_USE_PETSC},true)
  ifeq ($(OPERATING_SYSTEM),linux)# Linux
    PETSC_INCLUDE_PATH += $(addprefix -I, $(OCE_CM_INSTALL_DIR)/ )
    PETSC_INCLUDE_PATH += $(addprefix -I, $(OCE_CM_INSTALL_DIR)/include )
    PETSC_INCLUDE_PATH += $(addprefix -I, $(OCE_CM_INSTALL_DIR)/conf )
    PETSC_TSLIB = $(call searchdirs, $(OCE_CM_INSTALL_DIR)/lib, libpetscts.a )
    ifeq (,$(PETSC_TSLIB))
      PETSC_LIBRARIES += -lpetsc -lMrm -lXmu -lXm -lXt -lX11 -lstdc++ -ldl #Prob newer version of PETSc which has rolled all the libraries into one.
    else
      PETSC_LIBRARIES += -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lpetsccontrib -lMrm -lXmu -lXm -lXt -lX11 -lstdc++ -ldl
    endif
    PETSC_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
    ifeq ($(ABI),64)
      PETSC_LIB_PATH += $(addprefix -L, /usr/X11R6/lib64 )
      PETSC_LIB_PATH += $(addprefix -L, /usr/lib64 )
      PETSC_LIB_PATH += $(addprefix -L, /usr/lib64/openmotif )
    else
      PETSC_LIB_PATH += $(addprefix -L, /usr/X11R6/lib )
      PETSC_LIB_PATH += $(addprefix -L, /usr/lib )
      PETSC_LIB_PATH += $(addprefix -L, /usr/lib/openmotif )
    endif
  else
    ifeq ($(OPERATING_SYSTEM),aix)# AIX
      PETSC_INCLUDE_PATH += $(addprefix -I, $(OCE_CM_INSTALL_DIR)/ )
      PETSC_INCLUDE_PATH += $(addprefix -I, $(OCE_CM_INSTALL_DIR)/include )
      PETSC_INCLUDE_PATH += $(addprefix -I, $(OCE_CM_INSTALL_DIR)/conf )
      PETSC_TSLIB = $(call searchdirs, $(OCE_CM_INSTALL_DIR)/lib, libpetscts.a )
      ifeq (,$(PETSC_TSLIB))
        PETSC_LIBRARIES += -lpetsc #Prob newer version of PETSc which has rolled all the libraries into one.
      else
        PETSC_LIBRARIES += -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lpetsccontrib -lMrm -lXmu -lXm -lXt -lX11 -lstdc++
      endif
      PETSC_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
      PETSC_LIB_PATH += $(addprefix -L, /usr/X11R6/lib )
    else# windows
      PETSC_LIB_PATH += $(addprefix -L, /home/users/local/lib )
      PETSC_INCLUDE_PATH += $(addprefix -I, /home/users/local )
    endif
  endif
endif

#UMFPACK
UMFPACK_INCLUDE_PATH =#
UMFPACK_LIBRARIES =#
UMFPACK_LIB_PATH =#
ifeq (${OPENCMISS_USE_UMFPACK},true)
  UMFPACK_LIBRARIES += -lumfpack -lcholmod -lamd -lcamd -lccolamd -lcolamd -lsuitesparseconfig
  UMFPACK_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#SUNDIALS
SUNDIALS_INCLUDE_PATH =#
SUNDIALS_LIBRARIES =#
SUNDIALS_LIB_PATH =#
ifeq (${OPENCMISS_USE_SUNDIALS},true)
  SUNDIALS_LIBRARIES += -lsundials_cvode -lsundials_nvecparallel -lsundials_nvecserial
  SUNDIALS_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#HYPRE
HYPRE_INCLUDE_PATH =#
HYPRE_LIBRARIES =#
HYPRE_LIB_PATH =#
ifeq (${OPENCMISS_USE_HYPRE},true)
  HYPRE_LIBRARIES += -lHYPRE
  HYPRE_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#PaStiX and Scotch
PASTIX_INCLUDE_PATH =#
PASTIX_LIBRARIES =#
PASTIX_LIB_PATH =#
ifeq (${OPENCMISS_USE_PASTIX},true)
  PASTIX_LIBRARIES += -lpastix -lptscotch -lptscotcherr -lptscotcherrexit -lscotch -lscotcherr -lscotcherrexit -lrt
  PASTIX_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#SuperLU_DIST and SuperLU
SUPERLU_INCLUDE_PATH =#
SUPERLU_LIBRARIES =#
SUPERLU_LIB_PATH =#
ifeq (${OPENCMISS_USE_SUPERLU},true)
  SUPERLU_LIBRARIES += -lsuperlu -lsuperlu_dist
  SUPERLU_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#PLAPACK
PLAPACK_INCLUDE_PATH =#
PLAPACK_LIBRARIES =#
PLAPACK_LIB_PATH =#
ifeq (${OPENCMISS_USE_PLAPACK},true)
  PLAPACK_LIBRARIES += -lPLAPACK
  PLAPACK_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#MUMPS, ScaLAPACK and BLACS
MUMPS_INCLUDE_PATH =#
MUMPS_LIBRARIES =#
MUMPS_LIB_PATH =#
ifeq (${OPENCMISS_USE_MUMPS},true)
  MUMPS_LIBRARIES += -ldmumps -lmumps_common $(call convertlistlink,$(OCE_SCALAPACKLIBS)) -lblacs $(call convertlistlink,$(OCE_BLACSLIBS)) -lpord
  MUMPS_LIB_PATH += $(addprefix -L, $(OCE_BLACSLIBPATH) )
  MUMPS_LIB_PATH += $(addprefix -L, $(OCE_SCALAPACKLIBPATH) )
  MUMPS_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#ParMETIS
PARMETIS_INCLUDE_PATH = #
PARMETIS_LIBRARIES = #
PARMETIS_LIB_PATH =#
ifeq (${OPENCMISS_USE_PARMETIS},true)
  PARMETIS_LIBRARIES += -lparmetis -lmetis
  PARMETIS_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#CUBPACK
CUBPACK_INCLUDE_PATH =#
CUBPACK_LIBRARIES = #
CUBPAK_LIB_PATH =#
ifeq (${OPENCMISS_USE_CUBPACK},true)
  CUBPACK_LIBRARIES = -lcubpack
  CUBPACK_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#TAU
TAU_INCLUDE_PATH =#
TAU_LIBRARIES =#
TAU_LIB_PATH =#
ifeq (${OPENCMISS_USE_TAU},true)
  CPPFLAGS += -DTAUPROF -DOPENCMISS_USE_TAU
  FPPFLAGS += -DTAUPROF -DOPENCMISS_USE_TAU
  TAU_LIBRARIES += -lTAU_tf
  ifeq (${OPENCMISS_USE_PAPI},true)
    TAU_LIBRARIES += -lTauMpi-memory-phase-papi-mpi -ltau-memory-phase-papi-mpi
  else
    TAU_LIBRARIES += -lTauMpi-memory-phase-mpi -ltau-memory-phase-mpi
  endif
  TAU_LIB_PATH += $(addprefix -L, $(OCE_CM_INSTALL_DIR)/lib )
endif

#PAPI
PAPI_INCLUDE_PATH =#
PAPI_LIBRARIES =#
PAPI_LIB_PATH =#
ifeq (${OPENCMISS_USE_PAPI},true)
  PAPI_LIBRARIES += $(call convertlistlink,$(OCE_PAPILIBS))
  PAPI_INCLUDE_PATH += $(addprefix -I, $(OCE_PAPIINCLUDEPATH) )
  PAPI_LIB_PATH += $(addprefix -L, $(OCE_PAPILIBPATH) )
endif

#CELLML
CELLML_INCLUDE_PATH =#
CELLML_LIBRARIES =#
CELLML_LIB_PATH =#
ifeq (${OPENCMISS_USE_CELLML},true)
  CPPFLAGS += -DUSECELLML -DOPENCMISS_USE_CELLML
  FPPFLAGS += -DUSECELLML -DOPENCMISS_USE_CELLML
  CELLML_LIBRARIES += -loccellml -lcellml -lccgs -lcevas -lcuses -lmalaes -lvacss -lannotools -lxml2 -lz
  CELLML_INCLUDE_PATH += $(addprefix -I, $(OC_CELLML_INSTALL_DIR)/include/ )
  CELLML_INCLUDE_PATH += $(addprefix -I, $(OCE_CELLML_INSTALL_DIR)/include )
  CELLML_LIB_PATH += $(addprefix -L, $(OC_CELLML_INSTALL_DIR)/lib )
  CELLML_LIB_PATH += $(addprefix -L, $(OCE_CELLML_INSTALL_DIR)/lib )
  CELLML_LIB_PATH += $(addprefix -L, $(OCE_COMMON_INSTALL_DIR)/lib )
endif

#FIELDML
FIELDML_INCLUDE_PATH =#
MOD_FIELDML_TARGET = #
FIELDML_LIBRARIES =#
FIELDML_LIB_PATH =#
ifeq (${OPENCMISS_USE_FIELDML},true)
  CPPFLAGS += -DUSEFIELDML -DOPENCMISS_USE_FIELDML
  FPPFLAGS += -DUSEFIELDML -DOPENCMISS_USE_FIELDML
  MOD_FIELDML_TARGET = MOD_FIELDML
  FIELDML_LIBRARIES += -lfieldml_api -lfieldml_io_api -lxml2 -lz
  FIELDML_INCLUDE_PATH += $(addprefix -I, $(OCE_FIELDML_INSTALL_DIR)/include )
  FIELDML_LIB_PATH += $(addprefix -L, $(OCE_FIELDML_INSTALL_DIR)/lib )
  FIELDML_LIB_PATH += $(addprefix -L, $(OCE_COMMON_INSTALL_DIR)/lib )
endif

#NetCDF
NETCDF_INCLUDE_PATH =#
NETCDF_LIBRARIES =#
NETCDF_LIB_PATH =#
ifeq (${OPENCMISS_USE_NETCDF},true)
  NETCDF_LIBRARIES += $(call convertlistlink,$(OCE_NETCDFLIBS))
  NETCDF_LIB_PATH += $(addprefix -L, $(OCE_NETCDFLIBPATH) )
  NETCDF_INCLUDE_PATH += $(addprefix -I, $(OCE_NETCDFINCLUDEPATH) )
endif

#HDF5
HDF5_INCLUDE_PATH =#
HDF5_LIBRARIES =#
HDF5_LIB_PATH =#
ifeq (${OPENCMISS_USE_HDF5},true)
  HDF5_LIBRARIES += $(call convertlistlink,$(OCE_HDF5LIBS))
  HDF5_LIB_PATH += $(addprefix -L, $(OCE_HDF5LIBPATH) )
  HDF5_INCLUDE_PATH += $(addprefix -I, $(OCE_HDF5INCLUDEPATH) )
endif

#CUDA
CUDA_INCLUDE_PATH =#
CUDA_LIBRARIES =#
CUDA_LIB_PATH =#
ifeq (${OPENCMISS_USE_CUDA},true)
  CUDA_LIBRARIES += $(call convertlistlink,$(OCE_CUDALIBS))
  CUDA_INCLUDE_PATH += $(addprefix -L, $(OCE_CUDALIBPATH) )
  CUDA_INCLUDE_PATH += $(addprefix -I, $(OCE_CUDAINCLUDEPATH) )
endif

#MPI
MPI_LIB_PATH =#
MPI_LIBRARIES =#
MPI_INCLUDE_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  MPI_LIBRARIES += $(call convertlistlink,$(OCE_MPILIBS))
  MPI_LIB_PATH += $(addprefix -L, $(OCE_MPILIBPATH) )
  MPI_INCLUDE_PATH += $(addprefix -I, $(OCE_MPIINCLUDEPATH) )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    MPI_LIBRARIES += $(call convertlistlink,$(OCE_MPILIBS))
    MPI_LIB_PATH += $(addprefix -L, $(OCE_MPILIBPATH) )
    MPI_INCLUDE_PATH += $(addprefix -I, $(OCE_MPIINCLUDEPATH) )
  else# windows
    MPI_LIBRARIES += $(OCE_MPILIBRARIES)
    MPI_LIB_PATH += $(addprefix -L, $(OCE_MPILIBPATH) )
    MPI_INCLUDE_PATH += $(addprefix -I, $(OCE_MPIINSTALLPATH) )
  endif
endif

#BLAS/lapack
BLASLAPACK_LIB_PATH = $(addprefix -L, $(OCE_LAPACKLIBPATH) )
BLASLAPACK_LIB_PATH += $(call convertlistlibpaths,$(OCE_BLASLIBPATHF))
BLASLAPACK_LIBRARIES = $(call convertlistlink,$(OCE_LAPACKLIBS) $(OCE_BLASLIBSF))

EXTERNAL_LIBRARIES += $(strip $(SLEPC_LIBRARIES) $(TAO_LIBRARIES) $(PETSC_LIBRARIES) $(UMFPACK_LIBRARIES) $(SUNDIALS_LIBRARIES) $(HYPRE_LIBRARIES) $(PASTIX_LIBRARIES) $(SUPERLU_LIBRARIES) $(PLAPACK_LIBRARIES) $(MUMPS_LIBRARIES) $(PARMETIS_LIBRARIES) $(CUBPACK_LIBRARIES) $(TAU_LIBRARIES) $(PAPI_LIBRARIES) $(CELLML_LIBRARIES) $(FIELDML_LIBRARIES) $(NETCDF_LIBRARIES) $(HDF5_LIBRARIES) $(CUDA_LIBRARIES) $(MPI_LIBRARIES) $(BLASLAPACK_LIBRARIES))
#Note put MPI path in before NetCDF and HDF5 incase we have system NetCDF and HDF5 but OCE MPI.
EXTERNAL_LIB_PATH += $(strip $(CELLML_LIB_PATH) $(FIELDML_LIB_PATH) $(SLEPC_LIB_PATH) $(TAO_LIB_PATH) $(PETSC_LIB_PATH) $(UMFPACK_LIB_PATH) $(SUNDIALS_LIB_PATH) $(HYPRE_LIB_PATH) $(PASTIX_LIB_PATH) $(SUPERLU_LIB_PATH) $(PLAPACK_LIB_PATH) $(MUMPS_LIB_PATH) $(PARMETIS_LIB_PATH) $(CUBPACK_LIB_PATH) $(TAU_LIB_PATH) $(PAPI_LIB_PATH) $(MPI_LIB_PATH) $(NETCDF_LIB_PATH) $(HDF5_LIB_PATH) $(CUDA_LIB_PATH) $(BLASLAPACK_LIB_PATH))
EXTERNAL_INCLUDE_PATH += $(strip $(SLEPC_INCLUDE_PATH) $(TAO_INCLUDE_PATH) $(PETSC_INCLUDE_PATH) $(UMFPACK_INCLUDE_PATH) $(SUNDIALS_INCLUDE_PATH) $(HYPRE_INCLUDE_PATH) $(PASTIX_INCLUDE_PATH) $(SUPERLU_INCLUDE_PATH) $(PLAPACK_INCLUDE_PATH) $(MUMPS_INCLUDE_PATH) $(PARMETIS_INCLUDE_PATH) $(CUBPACK_INCLUDE_PATH) $(TAU_INCLUDE_PATH)  $(PAPI_INCLUDE_PATH) $(CELLML_INCLUDE_PATH) $(FIELDML_INCLUDE_PATH) $(NETCDF_INCLUDE_PATH) $(MPI_INCLUDE_PATH) $(HDF5_INCLUDE_PATH) $(CUDA_INCLUDE_PATH) )

# vim:filetype=make
