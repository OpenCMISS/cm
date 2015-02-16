# - Try to find BLACS
#

find_path (BLACS_DIR lib/libblacs.a HINTS ENV BLACS_DIR PATHS $ENV{HOME}/blacs DOC "Blacs Directory")

IF(EXISTS ${BLACS_DIR}/lib/libblacs.a)
  SET(BLACS_FOUND YES)
  SET(BLACS_INCLUDES ${BLACS_DIR})
  SET(BLACS_INCLUDE_DIR ${BLACS_DIR}/include)
  list(APPEND BLACS_INCLUDES ${BLACS_INCLUDE_DIR})
  FILE(GLOB BLACS_LIBRARIES RELATIVE "${BLACS_DIR}/lib" "${BLACS_DIR}/lib/libblacs*.a")
ELSE(EXISTS ${BLACS_DIR}/lib/libblacs.a)
  SET(BLACS_FOUND NO)
ENDIF(EXISTS ${BLACS_DIR}/lib/libblacs.a)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLACS DEFAULT_MSG BLACS_LIBRARIES BLACS_INCLUDES)