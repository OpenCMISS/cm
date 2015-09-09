#########################
# I know there's a UseSWIG.cmake file.
# This did not seem to be very flexible at the current stage.

find_package(PythonInterp QUIET)
if (NOT PYTHONINTERP_FOUND)
    message(WARNING "No Python interpreter found. Unable to generate any Iron bindings")
else()
    set(_CLEANUP)
    if (WITH_C_BINDINGS OR WITH_Python_BINDINGS)
        # Generate the c files in each case, as swig wraps again around the C interface
        set(OPENCMISS_H "${CMAKE_CURRENT_BINARY_DIR}/opencmiss.h")
        set(OPENCMISS_C_F90 "${CMAKE_CURRENT_BINARY_DIR}/opencmiss_c.f90")
        add_custom_command(OUTPUT ${OPENCMISS_C_F90} ${OPENCMISS_H}
                COMMAND "${PYTHON_EXECUTABLE}" generate_bindings "${Iron_SOURCE_DIR}" C "${OPENCMISS_H}" "${OPENCMISS_C_F90}"
                WORKING_DIRECTORY ${Iron_SOURCE_DIR}/bindings)
        list(APPEND _CLEANUP ${OPENCMISS_H} ${OPENCMISS_C_F90}) 
    endif()
    
    # C-Bindings extra target
    set(HAVE_C_BINDINGS FALSE)
    if (WITH_C_BINDINGS)
        message(STATUS "Creating C bindings")
        add_library(iron_c SHARED ${OPENCMISS_C_F90})
        target_link_libraries(iron_c iron)
        install(TARGETS iron_c
                    EXPORT iron-config
                    DESTINATION lib
                    INCLUDES DESTINATION include/iron)
        install(FILES ${OPENCMISS_H}
            DESTINATION include/iron)
        set(HAVE_C_BINDINGS TRUE)
    endif()
    
    # Python-Bindings extra target
    set(HAVE_Python_BINDINGS FALSE)
    if (WITH_Python_BINDINGS)
        message(STATUS "Creating Python bindings")
        find_package(PythonLibs QUIET)
        find_package(SWIG QUIET)
        # I know there's also a UseSWIG.cmake file (http://www.cmake.org/cmake/help/v3.3/module/UseSWIG.html)
        # This did not seem to be very flexible at the current stage.
        if (SWIG_FOUND)
            if(PYTHONLIBS_FOUND)
            
                set(IRON_PYTHON_MODULE iron_python)
                set(IRON_OPENCMISS_DIR ${CMAKE_CURRENT_BINARY_DIR}/opencmiss)
            
                # Copy interface files to binary dir
                set(SWIG_INTERFACE_FILES ${CMAKE_CURRENT_BINARY_DIR}/numpy.i ${CMAKE_CURRENT_BINARY_DIR}/numpy_extra.i ${CMAKE_CURRENT_BINARY_DIR}/iron.i)
                add_custom_command(OUTPUT ${SWIG_INTERFACE_FILES}
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/numpy.i ${CMAKE_CURRENT_BINARY_DIR}
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/numpy_extra.i ${CMAKE_CURRENT_BINARY_DIR}
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/iron.i ${CMAKE_CURRENT_BINARY_DIR}
                )
                list(APPEND _CLEANUP ${SWIG_INTERFACE_FILES})
     
                # iron_generated.i is included by iron.i later
                set(SWIG_IFACE "${CMAKE_CURRENT_BINARY_DIR}/iron_generated.i")
                add_custom_command(OUTPUT ${SWIG_IFACE}
                    COMMAND "${PYTHON_EXECUTABLE}" generate_bindings "${Iron_SOURCE_DIR}" SWIG "${SWIG_IFACE}"
                    COMMENT "Generating swig interface file for Iron"
                    WORKING_DIRECTORY ${Iron_SOURCE_DIR}/bindings
                )
    
                set(SWIG_OUTDIR "${IRON_OPENCMISS_DIR}/iron")
                # Generate C wrapper
                set(PYTHON_WRAPPER "${CMAKE_CURRENT_BINARY_DIR}/iron_python_wrapper.c")
                add_custom_command(OUTPUT ${PYTHON_WRAPPER}
                    DEPENDS ${SWIG_IFACE}
                    DEPENDS ${SWIG_INTERFACE_FILES}
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${SWIG_OUTDIR}
                    COMMAND ${SWIG_EXECUTABLE} -python -o "${PYTHON_WRAPPER}"
                        -module ${IRON_PYTHON_MODULE} -outdir ${SWIG_OUTDIR} iron.i
                    COMMENT "Generating python wrapper for Iron with swig"
                )
                list(APPEND _CLEANUP ${SWIG_IFACE} ${PYTHON_WRAPPER})
    
                #Generate wrapper object
                execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${Iron_SOURCE_DIR}/utils/numpy_include.py"
                    OUTPUT_VARIABLE NUMPY_INCLUDES
                    RESULT_VARIABLE RES_NUMPY_INC
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                )
                if (NOT RES_NUMPY_INC)
                    # MODULE: See http://www.cmake.org/cmake/help/v3.3/command/add_library.html
                    add_library(${IRON_PYTHON_MODULE} MODULE ${PYTHON_WRAPPER} ${OPENCMISS_H} ${OPENCMISS_C_F90})
                    target_link_libraries(${IRON_PYTHON_MODULE} PUBLIC iron ${PYTHON_LIBRARIES})
                    target_include_directories(${IRON_PYTHON_MODULE} PRIVATE
                        ${CMAKE_CURRENT_BINARY_DIR}
                        ${PYTHON_INCLUDE_DIRS}
                        ${NUMPY_INCLUDES}
                    )
                    set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES PREFIX ""
                        LIBRARY_OUTPUT_DIRECTORY ${SWIG_OUTDIR}
                        RUNTIME_OUTPUT_DIRECTORY ${SWIG_OUTDIR})
    
                    # We dont want any libXX or .dll - see UseSWIG.cmake
                    set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES PREFIX "_")
                    if(WIN32 AND NOT CYGWIN)
                        set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES SUFFIX ".pyd")
                    endif()
                        
                    # Generate Python module
                    set(IRON_PY ${CMAKE_CURRENT_BINARY_DIR}/iron.py)
    
                    set(PKGUTIL_INIT_PY ${IRON_OPENCMISS_DIR}/__init__.py)
                    set(INIT_PY ${SWIG_OUTDIR}/__init__.py)
                    set(UTILS_PY ${SWIG_OUTDIR}/_utils.py)
                    set(SETUP_PY ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
                    add_custom_command(TARGET ${IRON_PYTHON_MODULE} POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -DIRON_TARGET_FILE="$<TARGET_FILE:iron>" -DIRON_PYTHON_MODULE=${IRON_PYTHON_MODULE} -DUTILS_PY=${UTILS_PY} -DSETUP_PY=${SETUP_PY} -P ${CMAKE_CURRENT_SOURCE_DIR}/configurefiles.cmake
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/opencmiss/__init__.py ${PKGUTIL_INIT_PY}
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/opencmiss/iron/__init__.py ${INIT_PY}
                        COMMAND "${PYTHON_EXECUTABLE}" generate_bindings "${Iron_SOURCE_DIR}" Python ${IRON_PYTHON_MODULE} "${SWIG_OUTDIR}"
                        COMMENT "Generating Python top-level module"
                        WORKING_DIRECTORY "${Iron_SOURCE_DIR}/bindings"
                    )
                    list(APPEND _CLEANUP ${UTILS_PY} ${SETUP_PY} ${PKGUTIL_INIT_PY} ${INIT_PY})
                    install(TARGETS ${IRON_PYTHON_MODULE}
                        EXPORT iron-config
                        DESTINATION python/opencmiss
                    )
                    list(APPEND _CLEANUP ${SWIG_IFACE} ${PYTHON_WRAPPER})
                    
                    #Generate wrapper object
                    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${Iron_SOURCE_DIR}/utils/numpy_include.py"
                        OUTPUT_VARIABLE NUMPY_INCLUDES
                        RESULT_VARIABLE RES_NUMPY_INC
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
                    # "Bugfix": If we use the same source file name, there will be two processes messing with
                    # the same fortran module file which can lead to errors using -j parallel builds.
                    # By this we simply force to have iron_c build before the python module.
                    if (HAVE_C_BINDINGS)
                        add_dependencies(${IRON_PYTHON_MODULE} iron_c)
                    endif()
                        
                    list(APPEND _CLEANUP ${IRON_PY} ${SWIG_OUTDIR} ${IRON_OPENCMISS_DIR})
                else()
                    message(FATAL_ERROR "Could not detect NumPy includes")
                endif()
                
                set(HAVE_Python_BINDINGS TRUE)
            else()
                message(WARNING "No Python libraries found. Unable to generate Python bindings for Iron.")
            endif()
        else()
            message(WARNING "SWIG could not be found. Unable to generate Python bindings for Iron.")            
        endif()
    endif()
    if (_CLEANUP)
        set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${_CLEANUP}")
    endif()
endif()
