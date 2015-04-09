cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project(Iron VERSION 1.0 LANGUAGES C Fortran)

option(WITH_CELLML "Build iron with CellML support" YES)
option(WITH_FIELDML "Build iron with FieldML support" NO)
#option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
option(PACKAGE_CONFIG_DIR "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)" "lib/cmake")

# Config
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED YES)
SET(IRON_MODULES ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
INCLUDE(${IRON_MODULES}/MacroDefinitions.cmake)
INCLUDE(${IRON_MODULES}/CompileFlags.cmake)

# External components
# MPI
find_package(MPI REQUIRED)

SET(LINK_LIBS )
# CellML
if (WITH_CELLML)
    find_package(LIBCELLML CONFIG REQUIRED)
    find_package(CELLML-API CONFIG REQUIRED)
    LIST(APPEND LINK_LIBS cellml cellml_api cellml_model_definition)
endif()
# FieldML
if (WITH_FIELDML)
    find_package(FIELDML REQUIRED)
endif()

# Other Dependencies
include(Variables)
SET(DEPS PETSC SUNDIALS HYPRE MUMPS SCALAPACK) #BLACS
foreach(DEP ${DEPS})
    # Look for config-aware implementations first
    find_package(${DEP} CONFIG QUIET)
    if (NOT ${DEP}_FOUND)
        FIND_PACKAGE(${DEP} REQUIRED)
        INCLUDE_DIRECTORIES(${${DEP}_INCLUDES})
    else()
        LIST(APPEND LINK_LIBS ${${DEP}_TARGETS})
    endif()
endforeach()
message("Iron link libaries: ${LINK_LIBS}")

# Config
SET(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
#SET(HEADER_INC_NAME opencmiss.h)
#SET(HEADER_SOURCE_INC ${SOURCE_DIR}/${HEADER_INC_NAME})
SET(LIBRARY_INSTALL_TYPE ARCHIVE)

##############################
# Libraries / targets

# C-Bindings extra target
find_package(PythonInterp QUIET)
set(HAVE_BINDINGS FALSE)
if (PYTHONINTERP_FOUND)
    set(OPENCMISS_H ${CMAKE_CURRENT_BINARY_DIR}/opencmiss.h)
    set(OPENCMISS_C_F90 ${CMAKE_CURRENT_BINARY_DIR}/opencmiss_c.f90)
    set_source_files_properties(${OPENCMISS_C_F90} PROPERTIES GENERATED TRUE)
    add_custom_target(bindings 
        COMMAND ${PYTHON_EXECUTABLE} generate_bindings ${CMAKE_CURRENT_SOURCE_DIR} C ${OPENCMISS_H} ${OPENCMISS_C_F90}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bindings)
    set(HAVE_BINDINGS TRUE)
else()
    message(WARNING "No Python interpreter found. Unable to generate C bindings for Iron.")
endif()

IF(WIN32)
  # Need separate C/Fortran code due to Project natures in VS
  FILE(GLOB OPENCMISS_FORTRAN_SRCS ${SOURCE_DIR}/*.f90)
  FILE(GLOB OPENCMISS_C_SRCS ${SOURCE_DIR}/*.c)
  FILE(GLOB OPENCMISS_C_HEADER ${SOURCE_DIR}/*.h)
  LIST(APPEND OPENCMISS_C_SRCS ${OPENCMISS_C_HEADER})
  
  INCLUDE(${IRON_MODULES}/Entries.cmake ) 
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

  # Remove some source files from compilation
  FILE(GLOB IRON_SOURCES ${SOURCE_DIR}/*.*)
  INCLUDE(${IRON_MODULES}/Entries.cmake)
  FOREACH(infileName ${EXCLUDED_ROUTINES})
    FILE(GLOB infileName ${infileName}) # necessary to get an exact match of the filename in the list (for out-of-source builds)
    IF (infileName) # if file was found
        LIST(REMOVE_ITEM IRON_SOURCES ${infileName})
    ENDIF ()
  ENDFOREACH()
  
  if (HAVE_BINDINGS)
      list(APPEND IRON_SOURCES ${OPENCMISS_C_F90})
  endif()

  ADD_LIBRARY(iron ${IRON_SOURCES})
  if (WITH_CELLML)
      target_compile_definitions(iron PRIVATE USECELLML)
  endif()
  if (WITH_FIELDML)
      target_compile_definitions(iron PRIVATE USEFIELDML)
  endif()
  target_include_directories(iron PUBLIC ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(iron ${LINK_LIBS})
  
  if (HAVE_BINDINGS)
      add_dependencies(iron bindings)
  endif()
  
ENDIF()

# Install steps
INSTALL(TARGETS iron 
      EXPORT iron-config
      DESTINATION iron
      INCLUDES DESTINATION include)     
install(EXPORT iron-config DESTINATION ${PACKAGE_CONFIG_DIR})
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake DESTINATION ${PACKAGE_CONFIG_DIR})      

FILE(GLOB OPENCMISS_MODS ${CMAKE_BINARY_DIR}/*.mod)
INSTALL(FILES ${OPENCMISS_MODS} DESTINATION include)

if (HAVE_BINDINGS)
    INSTALL(FILES ${OPENCMISS_H} DESTINATION include)
endif()