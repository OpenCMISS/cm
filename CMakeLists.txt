# License

cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project(Iron VERSION 0.4.0 LANGUAGES C Fortran)

option(WITH_CELLML "Build iron with CellML support" ON)
option(WITH_FIELDML "Build iron with FieldML support" ON)
#option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
set(PACKAGE_CONFIG_DIR "cmake" CACHE STRING "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)")

# Config
# See http://cmake.3232098.n2.nabble.com/How-to-install-Fortran-module-files-td7584399.html
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules)
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)
SET(IRON_MODULES ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
LIST(APPEND CMAKE_MODULE_PATH
    ${IRON_MODULES}
)

INCLUDE(MacroDefinitions)
INCLUDE(CompileFlags)

# External components
# MPI
find_package(MPI REQUIRED)

SET(LINK_LIBS )
# CellML
message(STATUS "===${WITH_CELLML}===")
message(STATUS "===${WITH_FIELDML}===")
if (WITH_CELLML)
    find_package(LIBCELLML CONFIG REQUIRED)
    find_package(CELLML-API CONFIG REQUIRED)
    LIST(APPEND LINK_LIBS cellml cellml_model_definition)
endif()
# FieldML
if (WITH_FIELDML)
    find_package(FIELDML REQUIRED)
    LIST(APPEND LINK_LIBS fieldml_api fieldml_io_api)
endif()

# Other Dependencies
include(Variables)
SET(DEPS PETSC SUNDIALS HYPRE MUMPS SCALAPACK) #BLACS
foreach(DEP ${DEPS})
    # Look for config-aware implementations first
    find_package(${DEP} CONFIG QUIET)
    if (NOT ${DEP}_FOUND)
        message(STATUS "===${DEP}===")
        FIND_PACKAGE(${DEP} REQUIRED)
        INCLUDE_DIRECTORIES(${${DEP}_INCLUDES})
    else()
        message(STATUS "===CONFIG=${DEP}=${${DEP}_FOUND}= ===${${DEP}_TARGETS}===")
        LIST(APPEND LINK_LIBS ${${DEP}_TARGETS})
    endif()
endforeach()
message("Iron link libaries: ${LINK_LIBS}")
set(_expected_target)
set(_found_target)
set(_missing_target)
foreach(_link_lib ${LINK_LIBS})
    list(APPEND _expected_target ${_link_lib})
    if(TARGET ${_link_lib})
list(APPEND _found_target ${_link_lib})
else()
list(APPEND _missing_target ${_link_lib})
endif()
endforeach()
message(STATUS "==${_found_target}=")
message(STATUS "==${_missing_target}=")
unset(_expected_target)
#unset(_found_target)
unset(_missing_target)

# Config
SET(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
#SET(HEADER_INC_NAME opencmiss.h)
#SET(HEADER_SOURCE_INC ${SOURCE_DIR}/${HEADER_INC_NAME})
SET(LIBRARY_INSTALL_TYPE ARCHIVE)

# Remove some source files from compilation
FILE(GLOB IRON_SOURCES ${SOURCE_DIR}/*.*)
INCLUDE(Entries)
FOREACH(infileName ${EXCLUDED_ROUTINES})
    FILE(GLOB infileName ${infileName}) # necessary to get an exact match of the filename in the list (for out-of-source builds)
    IF (infileName) # if file was found
        LIST(REMOVE_ITEM IRON_SOURCES ${infileName})
    ENDIF (infileName)
ENDFOREACH(infileName)

# how portable is dumpversion?
execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion OUTPUT_VARIABLE Fortran_COMPILER_VERSION)
string(STRIP ${Fortran_COMPILER_VERSION} Fortran_COMPILER_VERSION)

# Libraries / targets
IF(WIN32)
  # Need separate C/Fortran code due to Project natures in VS
  FILE(GLOB OPENCMISS_FORTRAN_SRCS ${SOURCE_DIR}/*.f90)
  FILE(GLOB OPENCMISS_C_SRCS ${SOURCE_DIR}/*.c)
  FILE(GLOB OPENCMISS_C_HEADER ${SOURCE_DIR}/*.h)
  LIST(APPEND OPENCMISS_C_SRCS ${OPENCMISS_C_HEADER})

  INCLUDE(Entries)
  FOREACH(infileName ${EXCLUDED_ROUTINES})
    LIST(REMOVE_ITEM OPENCMISS_FORTRAN_SRCS ${infileName})
    LIST(REMOVE_ITEM OPENCMISS_C_SRCS ${infileName})
  ENDFOREACH(infileName)

  # ? maybe use object target an then combine to one target..
  ADD_LIBRARY(iron_Fortran ${OPENCMISS_FORTRAN_SRCS})
  target_include_directories(iron_Fortran PUBLIC ${MPI_Fortran_INCLUDE_PATH})
  if (WITH_CELLML)
      target_compile_definitions(iron_Fortran PUBLIC USECELLML)
  endif()
  ADD_LIBRARY(iron_C ${OPENCMISS_C_SRCS})
  if (WITH_CELLML)
      target_compile_definitions(iron_C PUBLIC USECELLML)
  endif()
  target_include_directories(iron_C PUBLIC ${MPI_C_INCLUDE_PATH})

ELSE()

  ADD_LIBRARY(iron ${IRON_SOURCES})
  if (WITH_CELLML)
      target_compile_definitions(iron PUBLIC USECELLML)
  endif()
  if (WITH_FIELDML)
      target_compile_definitions(iron PUBLIC USEFIELDML)
  endif()
  target_include_directories(iron PUBLIC ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
  target_include_directories(iron INTERFACE $<INSTALL_INTERFACE:fortran/${Fortran_COMPILER_VERSION}/opencmiss/iron>)
  target_link_libraries(iron PRIVATE ${_found_target}) # ${LINK_LIBS})
ENDIF()

# Install steps
INSTALL(TARGETS iron
      EXPORT iron-config
      DESTINATION lib/opencmiss
      INCLUDES DESTINATION include/opencmiss/iron)
#INSTALL(FILES ${HEADER_SOURCE_INC} DESTINATION include)
install(EXPORT iron-config DESTINATION ${PACKAGE_CONFIG_DIR})
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake DESTINATION ${PACKAGE_CONFIG_DIR})

FILE(GLOB OPENCMISS_MODS ${CMAKE_Fortran_MODULE_DIRECTORY}/*.mod)
INSTALL(FILES ${OPENCMISS_MODS} DESTINATION fortran/${Fortran_COMPILER_VERSION}/opencmiss/iron)

