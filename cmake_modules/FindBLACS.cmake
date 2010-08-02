# - Try to find BLACS
#

find_path (BLACS_DIR lib/libblacs.a HINTS ENV BLACS_DIR PATHS $ENV{HOME}/blacs DOC "Blacs Directory")


IF(EXISTS ${BLACS_DIR}/lib/libblacs.a)
  SET(BLACS_FOUND YES)
  SET(BLACS_INCLUDES ${BLACS_DIR})
  find_path (BLACS_INCLUDE_DIR mumps_compat.h HINTS "${BLACS_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND BLACS_INCLUDES ${BLACS_INCLUDE_DIR})
  FILE(GLOB BLACS_LIBRARIES RELATIVE "${BLACS_DIR}/lib" "${BLACS_DIR}/lib/libblacs*.a")
ELSE(EXISTS ${BLACS_DIR}/lib/libblacs.a)
  SET(BLACS_FOUND NO)
  message(FATAL_ERROR "Cannot find BLACS!")
ENDIF(EXISTS ${BLACS_DIR}/lib/libblacs.a)