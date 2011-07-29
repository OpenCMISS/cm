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

ifndef OPENCMISSEXTRAS_ROOT
  OPENCMISSEXTRAS_ROOT := ../../opencmissextras
endif

ifndef EXTERNAL_CM_ROOT
  EXTERNAL_CM_ROOT := ${OPENCMISSEXTRAS_ROOT}/cm/external
endif
ifndef EXTERNAL_CELLML_ROOT
  EXTERNAL_CELLML_ROOT := ${OPENCMISSEXTRAS_ROOT}/cellml
endif
ifndef EXTERNAL_FIELDML_ROOT
  EXTERNAL_FIELDML_ROOT := ${OPENCMISSEXTRAS_ROOT}/fieldml
endif
ifndef EXTERNAL_COMMON_ROOT
  EXTERNAL_COMMON_ROOT := ${OPENCMISSEXTRAS_ROOT}/common
endif

ifndef GLOBAL_CELLML_ROOT
  GLOBAL_CELLML_ROOT := ${GLOBAL_CM_ROOT}/../cellml
endif
ifndef GLOBAL_FIELDML_ROOT
  GLOBAL_FIELDML_ROOT := ${OPENCMISSEXTRAS_ROOT}/fieldml
endif
ifndef GLOBAL_UTILS_ROOT
  GLOBAL_UTILS_ROOT := ${OPENCMISSEXTRAS_ROOT}/utils
endif

include $(GLOBAL_UTILS_ROOT)/Makefile.inc

#----------------------------------------------------------------------------------------------------------------------------------

ifndef MPI
  MPI := mpich2
endif

ifndef USECELLML
  USECELLML := false
endif

ifndef USEFIELDML
  USEFIELDML := true
endif

ifeq ($(MPI),mpich2)
  MPI := mpich2
else
  ifeq ($(MPI),intel)
    ifeq ($(OPERATING_SYSTEM),linux)
      ifdef I_MPI_ROOT
        MPI := intel
      else
        $(error Intel MPI libraries not setup)
      endif
    else
      $(error can only use intel mpi with Linux)
    endif
  else
    ifeq ($(MPI),openmpi)
      MPI := openmpi
    else
      ifeq ($(MPI),mvapich2)
        MPI := mvapich2
      else
        ifeq ($(MPI),cray)
          MPI := cray
        else
          ifeq ($(MPI),poe)
            MPI := poe
          else
            $(error unknown MPI type - $(MPI))
          endif
        endif
      endif
    endif
  endif
endif

ifeq ($(MPIPROF),true)
  ifeq ($(MPI),intel)
    ifndef VT_ROOT
      $(error intel trace collector not setup)
    endif
    ifndef VT_ADD_LIBS
      $(error intel trace collector not setup)
    endif
  endif
endif

ifeq ($(OPERATING_SYSTEM),linux)# Linux
  ifdef COMPILER_VERSION
    ifndef EXTERNAL_CM_DIR
      EXTERNAL_CM_DIR := $(EXTERNAL_CM_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(MPI)/$(COMPILER)_$(COMPILER_VERSION)
    endif
    EXTERNAL_CELLML_DIR := $(EXTERNAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)_$(COMPILER_VERSION)
    EXTERNAL_FIELDML_DIR := $(EXTERNAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)_$(COMPILER_VERSION)
    EXTERNAL_COMMON_DIR := $(EXTERNAL_COMMON_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)_$(COMPILER_VERSION)
    MAIN_OBJECT_DIR := $(GLOBAL_CM_ROOT)/object/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(MPI)/$(COMPILER)_$(COMPILER_VERSION)
    INC_DIR := $(GLOBAL_CM_ROOT)/include/$(BIN_ARCH_DIR)/$(MPI)/$(COMPILER)_$(COMPILER_VERSION)
    LIB_DIR := $(GLOBAL_CM_ROOT)/lib/$(BIN_ARCH_DIR)/$(MPI)/$(COMPILER)_$(COMPILER_VERSION)
  else
    ifndef EXTERNAL_CM_DIR
      EXTERNAL_CM_DIR := $(EXTERNAL_CM_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(MPI)/$(COMPILER)
    endif
    EXTERNAL_CELLML_DIR := $(EXTERNAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)
    EXTERNAL_FIELDML_DIR := $(EXTERNAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)
    EXTERNAL_COMMON_DIR := $(EXTERNAL_COMMON_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)
    MAIN_OBJECT_DIR := $(GLOBAL_CM_ROOT)/object/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(MPI)/$(COMPILER)
    INC_DIR := $(GLOBAL_CM_ROOT)/include/$(BIN_ARCH_DIR)/$(MPI)/$(COMPILER)
    LIB_DIR := $(GLOBAL_CM_ROOT)/lib/$(BIN_ARCH_DIR)/$(MPI)/$(COMPILER)
  endif
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    ifndef EXTERNAL_CM_DIR
      EXTERNAL_CM_DIR := $(EXTERNAL_CM_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(MPI)/$(COMPILER)
    endif
    EXTERNAL_CELLML_DIR := $(EXTERNAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)
    EXTERNAL_FIELDML_DIR := $(EXTERNAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)
    EXTERNAL_COMMON_DIR := $(EXTERNAL_COMMON_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
  else# windows
    ifndef EXTERNAL_CM_DIR
      EXTERNAL_CM_DIR := $(EXTERNAL_CM_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
    endif
    EXTERNAL_CELLML_DIR := $(EXTERNAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
    EXTERNAL_FIELDML_DIR := $(EXTERNAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
    EXTERNAL_COMMON_DIR := $(EXTERNAL_COMMON_ROOT)/$(LIB_ARCH_DIR)$(DEBUG_SUFFIX)$(PROF_SUFFIX)
  endif
  MAIN_OBJECT_DIR := $(GLOBAL_CM_ROOT)/object/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(MPI)/$(COMPILER)
  INC_DIR := $(GLOBAL_CM_ROOT)/include/$(BIN_ARCH_DIR)/$(MPI)/$(COMPILER)
  LIB_DIR := $(GLOBAL_CM_ROOT)/lib/$(BIN_ARCH_DIR)/$(MPI)/$(COMPILER)
endif

#----------------------------------------------------------------------------------------------------------------------------------
# compiling commands

ifeq ($(MPI),mpich2)
  MPIFC = mpif90
  MPICC = mpicc
else
  MPIFC = mpiifort
  MPICC = mpiicc
endif
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
ifneq ($(TAUPROF),false)
  CPPFLAGS += -DTAUPROF #Necessary?
  FPPFLAGS += -DTAUPROF
endif
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

#----------------------------------------------------------------------------------------------------------------------------------
ifeq ($(OPERATING_SYSTEM),linux)
  OPTCF_FLGS =# Use separate flags for fortran and c
  olist_args = $1

  CC = gcc
  FC = gfortran

  ifeq ($(COMPILER),intel)

    #Use Intel compilers if available (icc -V sends output to STDERR and exits with error).
    ifneq (,$(shell icc -V 2>&1 | grep -i intel))
      CC = icc
    endif
    ifneq (,$(shell ifort -V 2>&1 | grep -i intel))
      FC = ifort
    endif

  endif

  # Set the flags for the various different CC compilers
  ifeq ($(CC),gcc)# gcc
    C_FLGS += -pipe -m$(ABI)
    # Position independent code is actually only required for objects
    # in shared libraries but debug version may be built as shared libraries.
    DBGC_FLGS += -fPIC
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        C_FLGS += -march=nocona
      endif
    endif
    DBGC_FLGS += -O0 -fbounds-check
    OPTC_FLGS = -O3 -funroll-all-loops
    ifeq ($(PROF),false)
      ifneq ($(filter $(INSTRUCTION),i686 x86_64),)# i686 or x86_64
        OPTC_FLGS += -momit-leaf-frame-pointer
      endif
    else
      C_FLGS += -g -pg -fprofile-arcs -ftest-coverage
      ELFLGS += -g -pg -fprofile-arcs -ftest-coverage
    endif
  endif
  ifeq ($(CC),icc)
    # Turn on all warnings
    C_FLGS += -Wall -m$(ABI)
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        ifneq ($(shell grep Duo /proc/cpuinfo 2>/dev/null),)
          ifneq ($(shell grep "Core(TM)2" /proc/cpuinfo 2>/dev/null),)
            C_FLGS += -xT# Core2 Duo
          else
            C_FLGS += -x0# Core Duo
          endif
        else
          C_FLGS += -xP# for sse3 (90nm Pentium 4 series)
        endif
      else
        C_FLGS += -xW# Pentium4 compatible (?sse2)
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
        C_FLGS += -xN# for Pentium 4
      endif
    endif
    DBGC_FLGS += -O0
    OPTC_FLGS = -O3 -ansi_alias
    ifneq ($(PROF),false)
      C_FLGS += -g -pg
      ELFLAGS += -pg
    endif
    ifeq ($(MPIPROF),true)
      ifeq ($(MPI),mpich2)
        C_FLGS += -Wl,--export-dyanmic
        DBGC_FLGS += -Wl,--export-dyanmic
      else
        C_FLGS += -tcollect
      endif
    endif
  endif
  ifeq ($(filter-out xlc%,$(CC)),)# xlc* C compiler
    CFLAGS += -qinfo=gen:ini:por:pro:trd:tru:use
    C_FLGS += -q$(ABI) -qarch=auto -qhalt=e
    # -qinitauto for C is bytewise: 7F gives large integers.
    DBGC_FLGS += -qfullpath -C -qflttrap=inv:en -qinitauto=7F
    OPTC_FLGS += -O3
    # for trailing _ on fortran symbols
    CPPFLAGS += -Dunix
  endif
  ifeq ($(CC),cc)# cray cc
    DBGC_FLGS += -O0 -g
    OPTC_FLGS = -O3
    ifeq ($(PROF),true)
      C_FLGS += -g -pg
    endif
  endif

  # Set the flags for the various different Fortran compilers
  ifeq ($(FC),gfortran)
    #FC = /home/users/local/packages/gfortran/irun/bin/gfortran
    # -fstatck-check
    F_FLGS += -pipe -m$(ABI) -fno-second-underscore -Wall -x f95-cpp-input
    # Restrict line length to 132
    F_FLGS += -ffree-line-length-132
    # for now change max identifier length. Should restrict to 63 (F2003) in future
    F_FLGS += -fmax-identifier-length=63
    # Position independent code is actually only required for objects
    # in shared libraries but debug version may be built as shared libraries.
    DBGF_FLGS += -fPIC
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        F_FLGS += -march=nocona
      endif
    endif
    DBGF_FLGS += -O0 -ffpe-trap=invalid,zero -fbacktrace
    ifdef COMPILER_VERSION
      ifeq ($(COMPILER_VERSION),4.5)
        DBGF_FLGS += -fcheck=all
      endif
      ifeq ($(COMPILER_VERSION),4.4)
        DBGF_FLGS += -fbounds-check
      endif
    endif
    OPTF_FLGS = -O3 -Wuninitialized -funroll-all-loops
    ifeq ($(PROF),false)
      ifneq ($(filter $(INSTRUCTION),i686 x86_64),)# i686 or x86_64
        OPTF_FLGS += -momit-leaf-frame-pointer
      endif
    else
      F_FLGS += -g -pg -fprofile-arcs -ftest-coverage
      ELFLAGS += -pg -fprofile-arcs -ftest-coverage
    endif
    ELFLAGS += -m$(ABI)
  endif
  ifeq ($(FC),g95)
    F_FLAGS += -fno-second-underscore -Wall -m$(ABI) -std=f2003
    DBGF_FLGS += -fPIC
    DBGF_FLGS += -O0 -fbounds-check
    OPTF_FLGS = -O3 -Wuninitialized -funroll-all-loops
    ELFLAGS += -m$(ABI)
    #$(error g95 not implemented)
  endif
  ifeq ($(FC),ifort)
    # turn on preprocessing,
    # turn on warnings,
    # warn about non-standard Fortran 95
    F_FLGS += -cpp -warn all -m$(ABI)
    ifeq ($(MACHNAME),x86_64)
      ifneq ($(shell grep Intel /proc/cpuinfo 2>/dev/null),)
        ifneq ($(shell grep Duo /proc/cpuinfo 2>/dev/null),)
          ifneq ($(shell grep "Core(TM)2" /proc/cpuinfo 2>/dev/null),)
            F_FLGS += -xT# Core2 Duo
          else
            F_FLGS += -x0# Core Duo
          endif
        else
          F_FLGS += -xP# for sse3 (90nm Pentium 4 series)
        endif
      else
        F_FLGS += -xW# Pentium4 compatible (?sse2)
      endif
    endif
    ifeq ($(filter-out i%86,$(MACHNAME)),)
      ifneq ($(shell grep sse2 /proc/cpuinfo 2>/dev/null),)
        F_FLGS += -xN# for Pentium 4
      endif
    endif
    DBGF_FLGS += -O0 -check all -traceback -debug all
    OPTF_FLGS = -O3
    ifneq ($(PROF),false)
      F_FLGS += -g -pg
      ELFLAGS += -pg
    endif
    ifeq ($(MPIPROF),true)
      ifeq ($(MPI),mpich2)
        F_FLAS += -Wl,--export-dyanmic
        DBGF_FLGS += -Wl,--export-dyanmic
      else
        F_FLGS += -tcollect
      endif
    endif
#    MP_FLGS = -openmp
    ELFLAGS += -nofor_main -m$(ABI) -traceback
  endif
  ifeq ($(filter-out xlf%,$(FC)),)# xlf* fortran compiler
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
  endif
  ifeq ($(FC),ftn)
    DBGF_FLGS += -O0 -g
    OPTF_FLGS = -O3
    ifeq ($(PROF),true)
      F_FLGS += -g -pg
      ELFLAGS += -pg
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
  FC = gfortran
  F_FLGS += -fno-second-underscore
  OPTCF_FLGS = -O2
  ELFLAGS += -Wl,-static
  # Use the ANSI C timers
  CPPFLAGS += -DANSI_C_TIMERS
  olist_args = $1
endif
ifeq ($(OPERATING_SYSTEM),aix)
  ifeq ($(MP),false)
    FC = mpxlf95
    CC = xlc
  else
    FC = mpxlf95_r
    CC = xlc_r
  endif
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

# This returns an empty string if not found
searchdirs = $(firstword $(wildcard $(addsuffix /$(strip $2),$1)))
# This still returns the name of the desired file if not found and so is useful for error checking and reporting.
searchdirsforce = $(firstword $(wildcard $(addsuffix /$(strip $2),$1)) $2)

# Check that call function works (for searchdirs, olist_args, etc.)
ifeq ($(call olist_args,test),)
  $(error call function not available.  Use GNU make version 3.78 or newer)
endif

#TAO
TAO_INCLUDE_PATH =#
TAO_LIBRARIES =#
TAO_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  TAO_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
  TAO_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  TAO_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/finclude )
  TAO_LIBRARIES = -ltaofortran -ltaopetsc -ltao
  TAO_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    TAO_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
    TAO_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
    TAO_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/finclude )
    TAO_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
  else# windows
    TAO_INCLUDE_PATH = $(addprefix -I, /home/users/local/ )
    TAO_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#PETSc
PETSC_INCLUDE_PATH =#
PETSC_LIBRARIES =#
PETSC_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  PETSC_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
  PETSC_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  PETSC_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/conf )
  PETSC_TSLIB = $(call searchdirs, $(EXTERNAL_CM_DIR)/lib/, libpetscts.a )
  ifeq (,$(PETSC_TSLIB))
    PETSC_LIBRARIES = -lpetsc -lMrm -lXmu -lXm -lXt -lX11 -lstdc++ -ldl #Prob newer version of PETSc which has rolled all the libraries into one.
  else
    PETSC_LIBRARIES = -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lpetsccontrib -lMrm -lXmu -lXm -lXt -lX11 -lstdc++ -ldl
  endif
  PETSC_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
  ifeq ($(ABI),64)
    PETSC_LIB_PATH += $(addprefix -L, /usr/X11R6/lib64/ )
    PETSC_LIB_PATH += $(addprefix -L, /usr/lib64/ )
    PETSC_LIB_PATH += $(addprefix -L, /usr/lib64/openmotif/ )
  else
    PETSC_LIB_PATH += $(addprefix -L, /usr/X11R6/lib/ )
    PETSC_LIB_PATH += $(addprefix -L, /usr/lib/ )
    PETSC_LIB_PATH += $(addprefix -L, /usr/lib/openmotif/ )
  endif
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    PETSC_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
    PETSC_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
    PETSC_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/conf/ )
    PETSC_TSLIB = $(call searchdirs, $(EXTERNAL_CM_DIR)/lib/, libpetscts.a )
    ifeq (,$(PETSC_TSLIB))
      PETSC_LIBRARIES = -lpetsc #Prob newer version of PETSc which has rolled all the libraries into one.
    else
      PETSC_LIBRARIES = -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lpetsccontrib -lMrm -lXmu -lXm -lXt -lX11 -lstdc++
    endif
    PETSC_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
    PETSC_LIB_PATH += $(addprefix -L, /usr/X11R6/lib/ )
  else# windows
    PETSC_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
    PETSC_INCLUDE_PATH = $(addprefix -I, /home/users/local/ )
  endif
endif

#SUNDIALS
SUNDIALS_INCLUDE_PATH =#
SUNDIALS_LIBRARIES =#
SUNDIALS_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  SUNDIALS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
  SUNDIALS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  SUNDIALS_LIBRARIES = -lsundials_cvode -lsundials_nvecparallel -lsundials_nvecserial
  SUNDIALS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    SUNDIALS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
    SUNDIALS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
    SUNDIALS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    SUNDIALS_INCLUDE_PATH = $(addprefix -I, /home/users/local/ )
    SUNDIALS_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#HYPRE
HYPRE_INCLUDE_PATH =#
HYPRE_LIBRARIES =#
HYPRE_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  HYPRE_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
  HYPRE_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  HYPRE_LIBRARIES = -lHYPRE
  HYPRE_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    HYPRE_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
    HYPRE_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
    HYPRE_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    HYPRE_INCLUDE_PATH = $(addprefix -I, /home/users/local/ )
    HYPRE_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#PaStiX
PASTIX_LIBRARIES =#
PASTIX_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  PASTIX_LIBRARIES = -lpastix
  PASTIX_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  PASTIX_LIBRARIES = -lpastix
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    PASTIX_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    PASTIX_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#Scotch
SCOTCH_LIBRARIES =#
SCOTCH_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  SCOTCH_LIBRARIES = -lptscotch -lptscotcherr -lptscotcherrexit -lscotch -lscotcherr -lscotcherrexit
  SCOTCH_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  SCOTCH_LIBRARIES = -lptscotch -lptscotchexit -lptscotcherrexit -lptscotchparmetis
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    SCOTCH_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    SCOTCH_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#SuperLU_DIST
SUPERLU_LIBRARIES =#
SUPERLU_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  SUPERLU_LIBRARIES = -lsuperlu_dist
  SUPERLU_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  SUPERLU_LIBRARIES = -lsuperlu_dist
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    SUPERLU_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    SUPERLU_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#MUMPS
MUMPS_INCLUDE_PATH =#
MUMPS_LIBRARIES =#
MUMPS_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  MUMPS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
  MUMPS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  MUMPS_LIBRARIES = -ldmumps -lmumps_common -lpord
  MUMPS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    MUMPS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
    MUMPS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
    MUMPS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    MUMPS_INCLUDE_PATH = $(addprefix -I, /home/users/local/ )
    MUMPS_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#ScaLAPACK
SCALAPACK_INCLUDE_PATH =#
SCALAPACK_LIBRARIES =#
SCALAPACK_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  SCALAPACK_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
  SCALAPACK_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  SCALAPACK_LIBRARIES = -lscalapack
  SCALAPACK_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    SCALAPACK_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
    SCALAPACK_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
    SCALAPACK_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    SCALAPACK_INCLUDE_PATH = $(addprefix -I, /home/users/local/ )
    SCALAPACK_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#BLACS
BLACS_INCLUDE_PATH =#
BLACS_LIBRARIES =#
BLACS_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  BLACS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
  BLACS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  BLACS_LIBRARIES = -lblacs -lblacsC -lblacsF77 -lblacs -lblacsC -lblacsF77
  BLACS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    BLACS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/ )
    BLACS_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
    BLACS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib )
  else# windows
    BLACS_INCLUDE_PATH = $(addprefix -I, /home/users/local/ )
    BLACS_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#ParMETIS
PARMETIS_INCLUDE_PATH =#
PARMETIS_LIBRARIES = -lparmetis -lmetis
PARMETIS_LIB_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  PARMETIS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    PARMETIS_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
  else# windows
    PARMETIS_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
  endif
endif

#TAU/PAPI
TAU_LIBRARIES =#
TAU_LIB_PATH =#
ifeq ($(TAUPROF),true)
  TAU_LIBRARIES += -ltau-memory-phase-papi-mpi -l:libpapi.a -l:libpfm.a
  ifeq ($(ABI),32)
    TAU_LIB_PATH += -L$(EXTERNAL_CM_DIR)/lib
  else
    TAU_LIB_PATH = -L$(EXTERNAL_CM_DIR)/lib -L$(EXTERNAL_CM_DIR)/lib64
  endif
endif

#CELLML
CELLML_INCLUDE_PATH =#
CELLML_LIBRARIES =#
CELLML_LIB_PATH =#
ifeq ($(USECELLML),true)
  CELLML_LIBRARIES += -loccellml -lcellml -lccgs -lcevas -lcuses -lmalaes -lvacss -lannotools -lxml2 -lz
  ifeq ($(OPERATING_SYSTEM),linux)# Linux
    ifdef COMPILER_VERSION
      CELLML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)_$(COMPILER_VERSION)/include/ )
     CELLML_LIB_PATH += $(addprefix -L, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)_$(COMPILER_VERSION)/lib/ )
    else
      CELLML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/include/ )
      CELLML_LIB_PATH += $(addprefix -L, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/lib/ )
    endif
    CELLML_LIB_PATH += $(addprefix -L, $(EXTERNAL_CELLML_DIR)/lib/ )
    CELLML_LIB_PATH += $(addprefix -L, $(EXTERNAL_COMMON_DIR)/lib/ )
  else
    ifeq ($(OPERATING_SYSTEM),aix)# AIX
         CELLML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/include/ )
         CELLML_LIB_PATH += $(addprefix -L, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/lib/ )
         CELLML_LIB_PATH += $(addprefix -L, $(EXTERNAL_CELLML_DIR)/lib/ )
         CELLML_LIB_PATH += $(addprefix -L, $(EXTERNAL_COMMON_DIR)/lib/ )
    else# windows
         CELLML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/include/ )
         CELLML_LIB_PATH += $(addprefix -L, $(GLOBAL_CELLML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/lib/ )
         CELLML_LIB_PATH += $(addprefix -L, $(EXTERNAL_CELLML_DIR)/lib/ )
         CELLML_LIB_PATH += $(addprefix -L, $(EXTERNAL_COMMON_DIR)/lib/ )
    endif
  endif
endif

#FIELDML
FIELDML_INCLUDE_PATH =#
MOD_FIELDML_TARGET = #
FIELDML_LIBRARIES =#
FIELDML_LIB_PATH =#
ifeq ($(USEFIELDML),true)
  MOD_FIELDML_TARGET = MOD_FIELDML
  FIELDML_LIBRARIES += -lfieldml -lxml2 -lz
  ifeq ($(OPERATING_SYSTEM),linux)# Linux
    ifdef COMPILER_VERSION
      FIELDML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)_$(COMPILER_VERSION)/include/ )
    else
      FIELDML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/include/ )
    endif
    FIELDML_LIB_PATH += $(addprefix -L, $(EXTERNAL_FIELDML_DIR)/lib/ )
    FIELDML_LIB_PATH += $(addprefix -L, $(EXTERNAL_COMMON_DIR)/lib/ )
    FIELDML_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_FIELDML_DIR)/include/ )
  else
    ifeq ($(OPERATING_SYSTEM),aix)# AIX
       FIELDML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/include/ )
      FIELDML_LIB_PATH += $(addprefix -L, $(EXTERNAL_FIELDML_DIR)/lib/ )
      FIELDML_LIB_PATH += $(addprefix -L, $(EXTERNAL_COMMON_DIR)/lib/ )
      FIELDML_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_FIELDML_DIR)/include/ )
    else# windows
       FIELDML_INCLUDE_PATH += $(addprefix -I, $(GLOBAL_FIELDML_ROOT)/$(LIB_ARCH_DIR)$(MT_SUFFIX)$(DEBUG_SUFFIX)$(PROF_SUFFIX)/$(COMPILER)/include/ )
      FIELDML_LIB_PATH += $(addprefix -L, $(EXTERNAL_FIELDML_DIR)/lib/ )
      FIELDML_LIB_PATH += $(addprefix -L, $(EXTERNAL_COMMON_DIR)/lib/ )
      FIELDML_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_FIELDML_DIR)/include/ )
    endif
  endif
endif

#TAU/PAPI
ifeq ($(TAUPROF),true)
  CPPFLAGS += -DTAUPROF
  FPPFLAGS += -DTAUPROF
endif

#MPI
MPI_LIB_PATH =#
MPI_LIBRARIES =#
MPI_INCLUDE_PATH =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  ifeq ($(MPI),mpich2)
    MPI_LIBRARIES += -lmpichf90
    ifeq ($(MPIPROF),true)
      MPI_LIBRARIES += -lmpe_f2cmpi -llmpe -lmpe
    endif
    MPL_LIB = $(call searchdirs, $(EXTERNAL_CM_DIR)/lib/, libmpl* )
    MPI_LIBRARIES += -lmpich
    ifneq (,$(MPL_LIB))
      MPI_LIBRARIES += -lmpl
    endif
    MPI_LIBRARIES += -lpthread -lrt
    MPI_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
    MPI_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
  else
    ifeq ($(MPI),intel)
      ifeq ($(MPIPROF),true)
        MPI_INCLUDE_PATH += $(addprefix -I, $(VT_ROOT)/include/ )
        MPI_LIB_PATH += $(addprefix -L, $(VT_LIB_DIR) )
        MPI_LIBRARIES += -lVT $(VT_ADD_LIBS)
      endif
      MPI_LIBRARIES += -lmpi -lmpiif -lmpigi -lrt -lpthread -ldl
      ifeq ($(ABI),64)
        MPI_INCLUDE_PATH += $(addprefix -I, $(I_MPI_ROOT)/include64/ )
        MPI_LIB_PATH += $(addprefix -L, $(I_MPI_ROOT)/lib64/ )
      else
        MPI_INCLUDE_PATH += $(addprefix -I, $(I_MPI_ROOT)/include/ )
        MPI_LIB_PATH += $(addprefix -L, $(I_MPI_ROOT)/lib/ )
      endif
    else
      MPI_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/lib/ )
      ifeq ($(MPI),mvapich2)
        MPI_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/include/ )
        MPI_LIBRARIES += -lmpich -lmpichf90 -lpthread  -lrt -luuid -lrt
        MPI_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
      else
        MPI_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/lib/ )
        MPI_LIBRARIES += -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil
        MPI_LIB_PATH += $(addprefix -L, $(EXTERNAL_CM_DIR)/lib/ )
      endif
    endif
  endif
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    MPI_LIBRARIES = -lmpi
    MPI_LIB_PATH += $(addprefix -L, /usr/lpp/ppe.poe/lib/ )
    MPI_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/lib/ )
  else# windows
    MPI_LIBRARIES = -lmpichf90 -lmpich -lpthread -lrt
    MPI_LIB_PATH += $(addprefix -L, /home/users/local/lib/ )
    MPI_INCLUDE_PATH += $(addprefix -I, $(EXTERNAL_CM_DIR)/lib/ )
  endif
endif

#BLAS/lapack
BLAS_INCLUDE_PATH =#
BLAS_LIB_PATH =#
BLAS_LIBRARIES =#
ifeq ($(OPERATING_SYSTEM),linux)# Linux
  BLAS_LIBRARIES = -llapack -lblas
else
  ifeq ($(OPERATING_SYSTEM),aix)# AIX
    BLAS_LIBRARIES = -lessl
  else# windows
    BLAS_LIBRARIES = -llapack -lblas
  endif
endif

EXTERNAL_LIBRARIES = #
EXTERNAL_LIB_PATH = #
EXTERNAL_INCLUDE_PATH = #

ifeq ($(USECELLML),true)
  EXTERNAL_INCLUDE_PATH += $(CELLML_INCLUDE_PATH)
  FPPFLAGS += -DUSECELLML
  CPPFLAGS += -DUSECELLML
  EXTERNAL_LIBRARIES += $(CELLML_LIBRARIES)
  EXTERNAL_LIB_PATH += $(CELLML_LIB_PATH)
endif

ifeq ($(USEFIELDML),true)
  EXTERNAL_LIBRARIES += $(FIELDML_LIBRARIES)
  EXTERNAL_LIB_PATH += $(FIELDML_LIB_PATH)
  EXTERNAL_INCLUDE_PATH += $(FIELDML_INCLUDE_PATH)
  FPPFLAGS += -DUSEFIELDML
  CPPFLAGS += -DUSEFIELDML
endif

EXTERNAL_LIBRARIES += $(strip $(TAO_LIBRARIES) $(PETSC_LIBRARIES) $(SUNDIALS_LIBRARIES) $(HYPRE_LIBRARIES) $(SUPERLU_LIBRARIES) $(PASTIX_LIBRARIES) $(SCOTCH_LIBRARIES) $(MUMPS_LIBRARIES) $(SCALAPACK_LIBRARIES) $(BLACS_LIBRARIES) $(PARMETIS_LIBRARIES) $(TAU_LIBRARIES) $(MPI_LIBRARIES) $(BLAS_LIBRARIES))
EXTERNAL_LIB_PATH += $(strip $(TAO_LIB_PATH) $(PETSC_LIB_PATH) $(SUNDIALS_LIB_PATH) $(HYPRE_LIB_PATH) $(SUPERLU_LIB_PATH) $(PASTIX_LIB_PATH) $(SCOTCH_LIB_PATH) $(MUMPS_LIB_PATH) $(SCALAPACK_LIB_PATH) $(BLACS_LIB_PATH) $(PARMETIS_LIB_PATH) $(TAU_LIB_PATH) $(MPI_LIB_PATH) $(BLAS_LIB_PATH))
EXTERNAL_INCLUDE_PATH += $(strip $(TAO_INCLUDE_PATH) $(PETSC_INCLUDE_PATH) $(SUNDIALS_INCLUDE_PATH) $(HYPRE_INCLUDE_PATH) $(MUMPS_INCLUDE_PATH) $(SCALAPACK_INCLUDE_PATH) $(BLACS_INCLUDE_PATH) $(PARMETIS_INCLUDE_PATH) $(MPI_INCLUDE_PATH) $(BLAS_INCLUDE_PATH))

# vim:filetype=make
