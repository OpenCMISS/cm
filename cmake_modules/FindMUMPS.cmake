# - Try to find MUMPS
#

find_path (MUMPS_DIR include/mumps_compat.h HINTS ENV MUMPS_DIR PATHS $ENV{HOME}/mumps DOC "Mumps Directory")

IF(EXISTS ${MUMPS_DIR}/include/mumps_compat.h)
  SET(MUMPS_FOUND YES)
  SET(MUMPS_INCLUDES ${MUMPS_DIR})
  find_path (MUMPS_INCLUDE_DIR mumps_compat.h HINTS "${MUMPS_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND MUMPS_INCLUDES ${MUMPS_INCLUDE_DIR})
  FILE(GLOB MUMPS_LIBRARIES RELATIVE "${MUMPS_DIR}/lib" "${MUMPS_DIR}/lib/libmumps*.a")
ELSE(EXISTS ${MUMPS_DIR}/include/mumps_compat.h)
  SET(MUMPS_FOUND NO)
  message(FATAL_ERROR "Cannot find MUMPS!")
ENDIF(EXISTS ${MUMPS_DIR}/include/mumps_compat.h)