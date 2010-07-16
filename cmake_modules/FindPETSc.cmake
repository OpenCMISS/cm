# - Try to find PETSc
# Once done this will define
#
# PETSC_FOUND - system has PETSc
# PETSC_INCLUDES - the PETSc include directories
# PETSC_LIBRARIES - Link these to use PETSc
# PETSC_COMPILER - Compiler used by PETSc, helpful to find a compatible MPI
# PETSC_DEFINITIONS - Compiler switches for using PETSc
# PETSC_MPIEXEC - Executable for running MPI programs
# PETSC_VERSION - Version string (MAJOR.MINOR.SUBMINOR)
#
# Hack: PETSC_VERSION currently decides on the version based on the
# layout. Otherwise we need to run C code to determine the version.
#
# Setting these changes the behavior of the search
# PETSC_DIR - directory in which PETSc resides
# PETSC_ARCH - build architecture
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

#PETSc
IF (${OPERATING_SYSTEM} MATCHES linux)# Linux
  INCLUDE_DIRECTORIES(${EXTERNAL_CM_ROOT}/x86_64-linux-debug/mpich2/gnu/)
  INCLUDE_DIRECTORIES(${EXTERNAL_CM_ROOT}/x86_64-linux-debug/mpich2/gnu/include/)
  INCLUDE_DIRECTORIES(${EXTERNAL_CM_ROOT}/x86_64-linux-debug/mpich2/gnu/conf/)
ELSEIF (${OPERATING_SYSTEM} MATCHES aix)# AIX
  SET(PETSC_INCLUDE_PATH
    -I${EXTERNAL_CM_ROOT}/x86_64-linux-debug/mpich2/gnu
    -I${EXTERNAL_CM_ROOT}/x86_64-linux-debug/mpich2/gnu/include
    -I${EXTERNAL_CM_ROOT}/x86_64-linux-debug/mpich2/gnu/include/finclude)
ELSE (${OPERATING_SYSTEM} MATCHES linux)# windows
    SET(PETSC_INCLUDE_PATH
      -L/home/users/local/lib
      -I/home/users/local)
ENDIF(${OPERATING_SYSTEM} MATCHES linux)# Linux