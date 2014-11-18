cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project(Iron VERSION 1.0 LANGUAGES C Fortran)

#################Define OPENCMISS_ROOT & OPENCMISSEXTRAS_ROOT##################

# Set OPENCMISS_ROOT
IF(NOT DEFINED OPENCMISS_ROOT)
  IF(DEFINED ENV{OPENCMISS_ROOT})
    SET(OPENCMISS_ROOT $ENV{OPENCMISS_ROOT})
  ELSE()
    SET(OPENCMISS_ROOT ..)
  ENDIF()
ENDIF() 
SET(GLOBAL_CELLML_ROOT ${OPENCMISS_ROOT}/cellml)

# Need some way of getting OPENCMISSEXTRAS_ROOT defined.
IF(NOT DEFINED OPENCMISSEXTRAS_ROOT)
  IF(DEFINED ENV{OPENCMISSEXTRAS_ROOT})
    set (OPENCMISSEXTRAS_ROOT $ENV{OPENCMISSEXTRAS_ROOT})
  ELSE()
    set (OPENCMISSEXTRAS_ROOT ${OPENCMISS_ROOT}/../opencmissextras) 
  ENDIF()
ENDIF()
SET(EXTERNAL_CM_ROOT ${OPENCMISSEXTRAS_ROOT}/cm/external)
SET(GLOBAL_FIELDML_ROOT ${OPENCMISSEXTRAS_ROOT}/fieldml)

#LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")
SET(IRON_MODULES ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
INCLUDE(${IRON_MODULES}/MacroDefinitions.cmake)
INCLUDE(${IRON_MODULES}/CompileFlags.cmake)

########################Include opencmissextras libraries#############################################

#IF(${OPERATING_SYSTEM} STREQUAL win32)
#  #TODO This is to be removed while new opencmissextras is available.
#  list(APPEND CMAKE_SYSTEM_PREFIX_PATH ${EXTERNAL_CM_ROOT}/i386-win32-debug)
#ELSE()
#  list(APPEND CMAKE_SYSTEM_PREFIX_PATH ${EXTERNAL_CM_ROOT}/${LIB_ARCH_DIR}${DEBUG_SUFFIX}${PROF_SUFFIX}/${MPI}/${COMPILER})
#ENDIF()

# MPI
IF(NOT ${OPERATING_SYSTEM} STREQUAL aix)
  FIND_PACKAGE(MPI REQUIRED)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Todo: write a script that parses the config files for added targets
SET(PARMETIS_TARGETS parmetis metis)
SET(PTSCOTCH_TARGETS scotch ptscotch)
SET(SUITESPARSE_TARGETS suitesparseconfig amd btf camd cholmod colamd ccolamd klu umfpack)
SET(PASTIX_TARGETS pastix)
SET(SCALAPACK_TARGETS scalapack)
SET(MUMPS_TARGETS mumps)
SET(SUPERLU_DIST_TARGETS superlu_dist)
SET(SUPERLU_TARGETS superlu)
SET(SUNDIALS_TARGETS sundials_cvode sundials_fcvode sundials_cvodes
    sundials_ida sundials_fida sundials_idas
    sundials_kinsol sundials_fkinsol
    sundials_nvecparallel sundials_nvecserial
    )
SET(HYPRE_TARGETS hypre)
SET(PETSC_TARGETS petsc)
SET(DEPS PETSC SUNDIALS HYPRE MUMPS SCALAPACK) #BLACS
SET(LINK_LIBS )
foreach(DEP ${DEPS})
    # Look for config-aware implementations first
    FIND_PACKAGE(${DEP} CONFIG QUIET)
    if (NOT ${DEP}_FOUND)
        FIND_PACKAGE(${DEP} REQUIRED)
        INCLUDE_DIRECTORIES(${${DEP}_INCLUDES})
    else()
        LIST(APPEND LINK_LIBS ${${DEP}_TARGETS})
    endif()
endforeach()
message(${LINK_LIBS})

IF(USECELLML)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSECELLML ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSECELLML ")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DUSECELLML ")
  INCLUDE_DIRECTORIES(${GLOBAL_CELLML_ROOT}/${LIB_ARCH_DIR}${DEBUG_SUFFIX}${PROF_SUFFIX}/${COMPILER}/include)
ENDIF(USECELLML)

IF(USEFIELDML)
  INCLUDE_DIRECTORIES(${GLOBAL_FIELDML_ROOT}/${LIB_ARCH_DIR}${DEBUG_SUFFIX}${PROF_SUFFIX}/include)
  INCLUDE_DIRECTORIES(${GLOBAL_FIELDML_ROOT}/${LIB_ARCH_DIR}${DEBUG_SUFFIX}${PROF_SUFFIX}/${COMPILER}/include)
ENDIF(USEFIELDML)

SET(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
SET(OBJECT_DIR object/${LIB_ARCH_DIR}${MY_SUFFIX}${DEBUG_SUFFIX}${PROF_SUFFIX}/${MPI}/${COMPILER})
SET(MODULE_DIR ${OBJECT_DIR})
SET(INC_DIR include/${BIN_ARCH_DIR}/${MPI}/${COMPILER})
SET(MOD_INC_NAME opencmiss.mod)
SET(MOD_SOURCE_INC ${CMAKE_BINARY_DIR}/${MOD_INC_NAME})
SET(HEADER_INC_NAME opencmiss.h)
SET(HEADER_SOURCE_INC ${SOURCE_DIR}/${HEADER_INC_NAME})
SET(LIB_DIR lib/${BIN_ARCH_DIR}/${MPI}/${COMPILER})
SET(LIBRARY_TARGET_NAME ${CMAKE_PROJECT_NAME}${EXE_ABI_SUFFIX}${MT_SUFFIX}${DEBUG_SUFFIX}${PRO_SUFFIX})
SET(LIBRARY_BUILD_TYPE STATIC)
SET(LIBRARY_INSTALL_TYPE ARCHIVE)

FILE(GLOB OPENCMISS_SRCS ${SOURCE_DIR}/*.*)
INCLUDE(${IRON_MODULES}/Entries.cmake)
FOREACH(infileName ${EXCLUDED_ROUTINES})
    FILE(GLOB infileName ${infileName}) # necessary to get an exact match of the filename in the list (for out-of-source builds)
    IF (infileName) # if file was found
        LIST(REMOVE_ITEM OPENCMISS_SRCS ${infileName})
    ENDIF (infileName)
ENDFOREACH(infileName)

# Preprocess files mentioned in "PreprocessedSources.cmake"
#[[INCLUDE(${IRON_MODULES}/PreprocessedSources.cmake)
FOREACH(filenameSrc ${PREPROCESSED_SRCS})
    
    # Get full filename of raw, unexpanded src file
    SET(fullFilenameSrc "${SOURCE_DIR}/${filenameSrc}")
    
    # Generate filename for expanded (intermediate) file
    GET_FILENAME_COMPONENT(filenameBase ${filenameSrc} NAME_WE)
    GET_FILENAME_COMPONENT(filenameExtension ${filenameSrc} EXT)
    SET(fullFilenameExpanded "${OBJECT_DIR}/${filenameBase}-expanded${filenameExtension}")
    
    # Don't compile the raw, unexpanded source file...
    FILE(GLOB fullFilenameSrc ${fullFilenameSrc}) # necessary to get an exact match of the filename in the list (for out-of-source builds)
    LIST(REMOVE_ITEM OPENCMISS_SRCS ${fullFilenameSrc})
    
    # ...but compile the preprocessed (expanded) intermediate file instead!
    LIST(APPEND OPENCMISS_SRCS ${fullFilenameExpanded})
    
    # Add rule that tells cmake how to create the intermediate file
    ADD_CUSTOM_COMMAND(OUTPUT      ${fullFilenameExpanded}
                       COMMAND     mkdir
                         ARGS      -p ${OBJECT_DIR}
                         VERBATIM
                       COMMAND     m4
                         ARGS      --prefix-builtins ${fullFilenameSrc} > ${fullFilenameExpanded}
                         VERBATIM
                       DEPENDS     ${fullFilenameSrc}
    )
    
ENDFOREACH(filenameSrc)]]

IF(${OPERATING_SYSTEM} STREQUAL win32)
  FILE(GLOB OPENCMISS_FORTRAN_SRCS ${SOURCE_DIR}/*.f90)
  FILE(GLOB OPENCMISS_C_SRCS ${SOURCE_DIR}/*.c)
  FILE(GLOB OPENCMISS_C_HEADER ${SOURCE_DIR}/*.h)
  LIST(APPEND OPENCMISS_C_SRCS ${OPENCMISS_C_HEADER})
  INCLUDE(${IRON_MODULES}/Entries.cmake ) 
  FOREACH(infileName ${EXCLUDED_ROUTINES})
    LIST(REMOVE_ITEM OPENCMISS_FORTRAN_SRCS ${infileName})
    LIST(REMOVE_ITEM OPENCMISS_C_SRCS ${infileName})
  ENDFOREACH(infileName)
  ADD_LIBRARY(${LIBRARY_TARGET_NAME}_Fortran ${LIBRARY_BUILD_TYPE} ${OPENCMISS_FORTRAN_SRCS} ${PUBLIC_HDRS} ${PRIVATE_HDRS})
  ADD_LIBRARY(${LIBRARY_TARGET_NAME}_C ${LIBRARY_BUILD_TYPE} ${OPENCMISS_C_SRCS} ${PUBLIC_HDRS} ${PRIVATE_HDRS})
ELSE()
  ADD_LIBRARY(${LIBRARY_TARGET_NAME} ${LIBRARY_BUILD_TYPE} ${OPENCMISS_SRCS} ${PUBLIC_HDRS} ${PRIVATE_HDRS})
  target_link_libraries(${LIBRARY_TARGET_NAME} ${LINK_LIBS})
ENDIF()

IF(${OPERATING_SYSTEM} STREQUAL linux)
  INSTALL(TARGETS ${LIBRARY_TARGET_NAME} ${LIBRARY_INSTALL_TYPE} DESTINATION lib/x86_64-linux/mpich2/gnu)
  INSTALL(FILES ${MOD_SOURCE_INC} ${HEADER_SOURCE_INC} DESTINATION ${INC_DIR})

  FILE(GLOB OPENCMISS_MODS ${CMAKE_BINARY_DIR}/*.mod)
  INSTALL(FILES ${OPENCMISS_MODS} DESTINATION ${MODULE_DIR})
  FILE(GLOB OPENCMISS_OBJS ${CMAKE_BINARY_DIR}/CMakeFiles/${CMAKE_PROJECT_NAME}${EXE_ABI_SUFFIX}${MT_SUFFIX}${DEBUG_SUFFIX}${PRO_SUFFIX}.dir/src/*.o)
  INSTALL(FILES ${OPENCMISS_OBJS} DESTINATION ${OBJECT_DIR})
ENDIF(${OPERATING_SYSTEM} STREQUAL linux)