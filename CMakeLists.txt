# License

cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project(Iron VERSION 0.4.0 LANGUAGES C Fortran)

option(WITH_CELLML "Build iron with CellML support" ON)
option(WITH_FIELDML "Build iron with FieldML support" OFF)
option(WITH_HYPRE "Build iron with HYPRE integration" ON)
option(WITH_SUNDIALS "Build iron with SUNDIALS integration" ON)
option(WITH_MUMPS "Build iron with MUMPS integration" ON)
option(WITH_SCALAPACK "Build iron with SCALAPACK integration" ON)
option(WITH_PETSC "Build iron with PETSC integration" ON)
option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
set(PACKAGE_CONFIG_DIR "cmake" CACHE STRING "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)")

# Config
# See http://cmake.3232098.n2.nabble.com/How-to-install-Fortran-module-files-td7584399.html
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules)
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)
# Have cmake find our own scripts
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

INCLUDE(MacroDefinitions)
INCLUDE(CompileFlags)

# External components
# MPI
find_package(MPI REQUIRED)

# Collect libraries to configure for/link against
SET(LINK_LIBS )
if (WITH_CELLML)
    find_package(LIBCELLML CONFIG REQUIRED)
    find_package(CELLML CONFIG REQUIRED)
    message(STATUS "Building with CELLML (${LIBCELLML_DIR})")
    LIST(APPEND LINK_LIBS cellml cellml_model_definition)
endif()
if (WITH_FIELDML)
    message(STATUS "Building with FIELDML (${FIELDML_DIR})")
    find_package(FIELDML-API REQUIRED)
    LIST(APPEND LINK_LIBS fieldml_api fieldml_io_api)
endif()
# This include file is in the main "manage" directory.
# We can include this here as iron is an integrated part of the overall build system.
include(Variables)
SET(DEPS PETSC SUNDIALS HYPRE MUMPS SCALAPACK) #BLACS
foreach(DEP ${DEPS})
    if (WITH_${DEP})
        # Look for config-aware implementations first
        find_package(${DEP} CONFIG QUIET)
        if (NOT ${DEP}_FOUND)
            #message(STATUS "=== Looking for ${DEP} - (wrapper) module mode ===")
            find_package(${DEP} REQUIRED)
            # This is probably not a complete list of variables that might be set by module-mode package finds :-(
            INCLUDE_DIRECTORIES(${${DEP}_INCLUDES} ${${DEP}_INCLUDE_DIRS} ${${DEP}_INCLUDE_DIRECTORY})
            message(STATUS "Building with ${DEP} (${${DEP}_LIBRARIES}})")
        else()
            message(STATUS "Building with ${DEP} (${${DEP}_DIR})")
            #message(STATUS "${DEP} targets===CONFIG=${DEP}=${${DEP}_FOUND}= ===${${DEP}_TARGETS}===")
            LIST(APPEND LINK_LIBS ${${DEP}_TARGETS})
        endif()
    endif()
endforeach()
#message("Iron link libaries: ${LINK_LIBS}")
#set(_expected_target)
#set(_found_target)
#set(_missing_target)
#foreach(_link_lib ${LINK_LIBS})
#    list(APPEND _expected_target ${_link_lib})
#    if(TARGET ${_link_lib})
#list(APPEND _found_target ${_link_lib})
#else()
#list(APPEND _missing_target ${_link_lib})
#endif()
#endforeach()
#message(STATUS "==${_found_target}=")
#message(STATUS "==${_missing_target}=")
#unset(_expected_target)
#unset(_found_target)
#unset(_missing_target)

# Config
#SET(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
#SET(HEADER_INC_NAME opencmiss.h)
#SET(HEADER_SOURCE_INC ${SOURCE_DIR}/${HEADER_INC_NAME})
SET(LIBRARY_INSTALL_TYPE ARCHIVE)

# Populates the IRON_XX_SRC variables
include(Sources)
IF(WIN32)
    # Need separate C/Fortran code due to Project natures in VS
    add_library(iron_Fortran ${IRON_Fortran_SRC})
    add_library(iron_C ${IRON_C_SRC})
    set(IRON_TARGETS iron_Fortran iron_C) 
    # Add joining library to have the "iron" target in each case 
    add_library(iron <$TARGET_OBJECTS:iron_Fortran> <$TARGET_OBJECTS:iron_C>)
ELSE()
    add_library(iron ${IRON_SRC})
    set(IRON_TARGETS iron)
ENDIF()
# Set the same stuff for either the separate targets (windows) or one target (else)
foreach(irontarget ${IRON_TARGETS})
    target_include_directories(${irontarget} 
        PUBLIC ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
        PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
    if (WITH_CELLML)
        target_compile_definitions(${irontarget} PUBLIC USECELLML)
    endif()
    if (WITH_FIELDML)
        target_compile_definitions(${irontarget} PUBLIC USEFIELDML)
    endif()
    target_link_libraries(${irontarget} PRIVATE ${LINK_LIBS}) # ${LINK_LIBS})
endforeach()

# Install steps
INSTALL(TARGETS iron
      EXPORT iron-config
      DESTINATION lib
      INCLUDES DESTINATION include/iron)
#INSTALL(FILES ${HEADER_SOURCE_INC} DESTINATION include)

install(EXPORT iron-config DESTINATION ${PACKAGE_CONFIG_DIR})
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake DESTINATION ${PACKAGE_CONFIG_DIR})

# how portable is dumpversion? NOT portable - and we have a toolchain architecture path segment
#execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion OUTPUT_VARIABLE Fortran_COMPILER_VERSION)
#string(STRIP ${Fortran_COMPILER_VERSION} Fortran_COMPILER_VERSION)
INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/iron
    PATTERN *.mod) #fortran/${Fortran_COMPILER_VERSION}/opencmiss/iron
