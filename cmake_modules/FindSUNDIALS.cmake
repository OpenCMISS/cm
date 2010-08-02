# - Try to find SUNDIALS
#

find_path (SUNDIALS_DIR include/sundials/sundials_config.h HINTS ENV SUNDIALS_DIR PATHS $ENV{HOME}/sundials DOC "Sundials Directory")

IF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)
  SET(SUNDIALS_FOUND YES)
  SET(SUNDIALS_INCLUDES ${SUNDIALS_DIR})
  find_path (SUNDIALS_INCLUDE_DIR sundials_config.h HINTS "${SUNDIALS_DIR}" PATH_SUFFIXES include/sundials NO_DEFAULT_PATH)
  list(APPEND SUNDIALS_INCLUDES ${SUNDIALS_INCLUDE_DIR})
  FILE(GLOB SUNDIALS_LIBRARIES RELATIVE "${SUNDIALS_DIR}/lib" "${SUNDIALS_DIR}/lib/libsundials*.a")
ELSE(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)
  SET(SUNDIALS_FOUND NO)
  message(FATAL_ERROR "Cannot find SUNDIALS!")
ENDIF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)