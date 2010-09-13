# - Try to find PETSc
# Once done this will define
#
# PETSC_FOUND - system has PETSc
# PETSC_INCLUDES - the PETSc include directories
# PETSC_LIBRARIES - Link these to use PETSc
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

find_path (PETSC_DIR include/petsc.h
  HINTS ENV PETSC_DIR
  PATHS
  ${CMAKE_SYSTEM_PREFIX_PATH}
  DOC "PETSc Directory")
  
#TODO This is to be removed (for Petsc version 2.3.3 in windows only) 
#SET(PETSC_ARCH cygwin-c-debug)

# Determine whether the PETSc layout is old-style (through 2.3.3) or
# new-style (3.0.0)
if (EXISTS ${PETSC_DIR}/include/petscconf.h) # > 2.3.3
  set (PETSC_VERSION "3.0.0")
  find_path (PETSC_CONF_DIR rules HINTS "${PETSC_DIR}" PATH_SUFFIXES conf NO_DEFAULT_PATH)
  set (PETSC_FOUND Yes)
elseif (EXISTS ${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf.h) # <= 2.3.3
  set (PETSC_VERSION "2.3.3")
  find_path (PETSC_CONF_DIR petscrules HINTS "${PETSC_DIR}" PATH_SUFFIXES bmake/${PETSC_ARCH} NO_DEFAULT_PATH)
  set (PETSC_FOUND Yes)
else (EXISTS ${PETSC_DIR}/include/petscconf.h)
  set (PETSC_FOUND No)
endif (EXISTS ${PETSC_DIR}/include/petscconf.h)

IF(PETSC_FOUND)
  SET(PETSC_INCLUDES ${PETSC_DIR})
  find_path (PETSC_INCLUDE_DIR petscts.h HINTS "${PETSC_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND PETSC_INCLUDES ${PETSC_INCLUDE_DIR})
  list(APPEND PETSC_INCLUDES ${PETSC_CONF_DIR})
  FILE(GLOB PETSC_LIBRARIES RELATIVE "${PETSC_DIR}/lib" "${PETSC_DIR}/lib/libpetsc*.a")
ENDIF(PETSC_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDES)