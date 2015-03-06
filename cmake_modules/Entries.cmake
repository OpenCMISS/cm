IF(${OPERATING_SYSTEM} MATCHES linux)# Linux
    SET(MACHINE_ENTRY ${SOURCE_DIR}/machine_constants_linux.f90)
ELSEIF(${OPERATING_SYSTEM} MATCHES aix)#AIX
    SET(MACHINE_ENTRY ${SOURCE_DIR}/machine_constants_aix.f90)
ELSE(${OPERATING_SYSTEM} MATCHES linux)# windows
    SET(MACHINE_ENTRY ${SOURCE_DIR}/machine_constants_windows.f90)
ENDIF()
FILE(GLOB MACHINE_ENTRY ${MACHINE_ENTRY}) # necessary to get an exact match of the filename in the list (for out-of-source builds)

FILE(GLOB MACHINE_EXCLUDES "${SOURCE_DIR}" "${SOURCE_DIR}/machine_constants_*.f90")
LIST(REMOVE_ITEM MACHINE_EXCLUDES ${MACHINE_ENTRY})

FILE(GLOB FIELDML_EXCLUDES "${SOURCE_DIR}" "${SOURCE_DIR}/fieldml_*.f90")

SET(EXCLUDED_ROUTINES 
  ${SOURCE_DIR}/Helmholtz_TEMPLATE_equations_routines.f90
  ${SOURCE_DIR}/binary_file_f.f90
  ${SOURCE_DIR}/finite_element_routines.f90
  ${SOURCE_DIR}/binary_file_c.c
)
LIST(APPEND EXCLUDED_ROUTINES ${MACHINE_EXCLUDES})
# Kick out fieldml files if support is disabled
IF(NOT WITH_FIELDML)
  LIST(APPEND EXCLUDED_ROUTINES ${FIELDML_EXCLUDES})
ENDIF()
