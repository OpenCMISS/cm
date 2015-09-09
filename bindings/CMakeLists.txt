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
            
                # Copy interface files to binary dir
                file(GLOB SWIG_INTERFACES "${CMAKE_CURRENT_SOURCE_DIR}/python/*.i")
                file(COPY ${SWIG_INTERFACES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
                
                # iron_generated.i is included by iron.i later
                set(SWIG_IFACE "${CMAKE_CURRENT_BINARY_DIR}/iron_generated.i")
                add_custom_command(OUTPUT ${SWIG_IFACE}
                    COMMAND "${PYTHON_EXECUTABLE}" generate_bindings "${Iron_SOURCE_DIR}" SWIG "${SWIG_IFACE}"
                    COMMENT "Generating swig interface file for Iron"
                    WORKING_DIRECTORY ${Iron_SOURCE_DIR}/bindings
                )
                                        
                # Generate C wrapper
                set(PYTHON_WRAPPER "${CMAKE_CURRENT_BINARY_DIR}/iron_python_wrapper.c")
                add_custom_command(OUTPUT ${PYTHON_WRAPPER}
                    DEPENDS ${SWIG_IFACE}
                    COMMAND ${SWIG_EXECUTABLE} -python -o "${PYTHON_WRAPPER}"
                        -module ${IRON_PYTHON_MODULE} -outdir . iron.i
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
                    if (HAVE_C_BINDINGS)
                        set(C_BINDIGS_SRC $<TARGET_OBJECTS:iron_c>)
                    else()
                        set(C_BINDIGS_SRC ${OPENCMISS_C_F90})
                    endif()
                    add_library(${IRON_PYTHON_MODULE} MODULE ${PYTHON_WRAPPER} ${OPENCMISS_H} )
                    target_link_libraries(${IRON_PYTHON_MODULE} PUBLIC iron)
                    target_include_directories(${IRON_PYTHON_MODULE} PRIVATE
                        ${CMAKE_CURRENT_BINARY_DIR}
                        ${PYTHON_INCLUDE_DIRS}
                        ${NUMPY_INCLUDES}
                    )
                    
                    # We dont want any libXX or .dll - see UseSWIG.cmake
                    set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES PREFIX "_")
                    if(WIN32 AND NOT CYGWIN)
                        set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES SUFFIX ".pyd")
                    endif()
                        
                    # Generate Python module
                    set(IRON_PY ${CMAKE_CURRENT_BINARY_DIR}/Iron.py)
                    add_custom_command(TARGET ${IRON_PYTHON_MODULE} POST_BUILD
                        COMMAND "${PYTHON_EXECUTABLE}" generate_bindings "${Iron_SOURCE_DIR}" Python "${CMAKE_CURRENT_BINARY_DIR}"
                        COMMENT "Generating Python top-level module"
                        WORKING_DIRECTORY "${Iron_SOURCE_DIR}/bindings"
                    )
                    install(TARGETS ${IRON_PYTHON_MODULE}
                        EXPORT iron-config
                        DESTINATION python/opencmiss
                    )
                    
                    # Install the Python module
                    file(GLOB additional_opencmiss_python_files "${CMAKE_CURRENT_SOURCE_DIR}/python/opencmiss/*.py")
                    install(FILES ${IRON_PY} 
                        ${additional_opencmiss_python_files} 
                        ${CMAKE_CURRENT_BINARY_DIR}/${IRON_PYTHON_MODULE}.py
                        DESTINATION python/opencmiss
                    )
                    list(APPEND _CLEANUP ${IRON_PY})
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
        set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${_CLEANUP})
    endif()
endif()