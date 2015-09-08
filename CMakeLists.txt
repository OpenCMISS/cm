cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project(Iron VERSION 0.4.0 LANGUAGES C CXX Fortran)# Need all languages as iron's deps cover all (link phase) 

option(WITH_CELLML "Build iron with CellML support" ON)
option(WITH_FIELDML "Build iron with FieldML support" OFF)
option(WITH_HYPRE "Build iron with HYPRE integration" ON)
option(WITH_SUNDIALS "Build iron with SUNDIALS integration" ON)
option(WITH_MUMPS "Build iron with MUMPS integration" ON)
option(WITH_SCALAPACK "Build iron with SCALAPACK integration" ON)
option(WITH_PETSC "Build iron with PETSC integration" ON)
option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
option(WITH_PROFILING "${PROJECT_NAME} - Build with profiling flags" OFF)
option(WITH_C_BINDINGS "Build iron C bindings" ON)
option(WITH_Python_BINDINGS "Build iron Python bindings" ON)
set(PACKAGE_CONFIG_DIR "cmake" CACHE STRING "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)")

# Config
# See http://cmake.3232098.n2.nabble.com/How-to-install-Fortran-module-files-td7584399.html
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules)
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)
# Have cmake find our own scripts
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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
    message(STATUS "Building with FIELDML (${FIELDML-API_DIR})")
    find_package(FIELDML-API REQUIRED)
endif()

# This include file is in the main "manage" directory.
# We can include this here as iron is an integrated part of the overall build system.
SET(IRON_UPSTREAM_DEPS PETSC SUNDIALS HYPRE MUMPS SCALAPACK PARMETIS)
foreach(DEP ${IRON_UPSTREAM_DEPS})
    if (WITH_${DEP})
        message(STATUS "Looking for ${DEP}..")
        find_package(${DEP} REQUIRED)
        if (${DEP}_DIR)
            # This is set if CONFIG mode was used
            set(DEP_LOC ${${DEP}_DIR})
        else()
            # This is set if MODULE mode was used
            if (${DEP}_LIBRARIES)
                set(_ALTERNATIVE ${DEP}_LIBRARIES)
            else()
                set(_ALTERNATIVE ${DEP}_INCLUDE_PATH)
            endif()
            list(GET ${_ALTERNATIVE} 0 DEP_LOC)
            get_filename_component(DEP_LOC ${DEP_LOC} PATH)
        endif()
        message(STATUS "Looking for ${DEP}.. found at ${DEP_LOC}")
        string(TOLOWER ${DEP} targetname)
        LIST(APPEND LINK_LIBS ${targetname})
    endif()
endforeach()

# Populates the IRON_XX_SRC variables
include(Sources)
IF(WIN32 AND NOT MINGW)
    # Need separate C/Fortran code due to Project natures in VS
    add_library(iron_Fortran OBJECT ${IRON_Fortran_SRC})
    add_library(iron_C OBJECT ${IRON_C_SRC})
    target_compile_options(iron_Fortran PRIVATE -cpp)
    set(IRON_TARGETS iron_Fortran iron_C) 
    # Add joining library to have the "iron" target in each case 
    add_library(iron $<TARGET_OBJECTS:iron_Fortran> $<TARGET_OBJECTS:iron_C>)
ELSE()
    add_library(iron ${IRON_SRC})
    target_compile_options(iron PRIVATE -cpp)
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
endforeach()
target_link_libraries(iron PRIVATE ${LINK_LIBS})
# Dont use "USE MPI" instructions if the percompiled mpi module files are incompatible.
# This happens when the fortran compiler used to build mpi differs from the currently used one.
# This holds on windows (only pre-compiled mpi implementations available; MPICH2 and MSMPI)
# and on e.g. clusters with inappropriately configured toolchain environments.
if (MINGW OR NOT MPI_Fortran_MODULE_COMPATIBLE)
    target_compile_definitions(iron PRIVATE NOMPIMOD)
endif()

# Fieldml extra library - some source of iron include the module files,
# but the cmake dependency parser doesnt get the USE declarations if wrapped in preprocessor #ifdefs
# hence the build breaks if those files are compiled before the fieldml-module files.
if (WITH_FIELDML)
    add_library(iron-fieldml ${IRON_FIELDML_SRC})
    target_link_libraries(iron-fieldml fieldml-api)
    target_compile_definitions(iron-fieldml PUBLIC USEFIELDML)
    install(TARGETS iron-fieldml EXPORT iron-config DESTINATION lib)
    
    # Let the iron lib know what's going on
    add_dependencies(iron iron-fieldml)
    target_link_libraries(iron PUBLIC iron-fieldml)
endif()

# Build bindings as selected
add_subdirectory(bindings)

# Install steps
INSTALL(TARGETS iron
      EXPORT iron-config
      DESTINATION lib
      INCLUDES DESTINATION include/iron)

install(EXPORT iron-config DESTINATION ${PACKAGE_CONFIG_DIR})
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake COMPATIBILITY AnyNewerVersion)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake 
    "include(CMakeFindDependencyMacro)\r\n"
    "set(${PROJECT_NAME}_IMPORT_PREFIX \${_IMPORT_PREFIX})\r\n"
)
foreach(IRON_DEP ${IRON_UPSTREAM_DEPS})
    if (WITH_${IRON_DEP})
        file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(${IRON_DEP})\r\n")
    endif()
endforeach()
if (WITH_CELLML)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(CELLML)\r\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(LIBCELLML)\r\n")
endif()
if (WITH_FIELDML)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(FIELDML-API)\r\n")
endif()
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "set(_IMPORT_PREFIX \${${PROJECT_NAME}_IMPORT_PREFIX})")
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake
    DESTINATION ${PACKAGE_CONFIG_DIR})  

INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/iron
    PATTERN *.mod)

if (HAVE_C_BINDINGS)
    INSTALL(FILES ${OPENCMISS_H} DESTINATION include/iron)
endif()
