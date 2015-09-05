find_package(PythonInterp QUIET)
if (NOT PYTHONINTERP_FOUND)
    message(WARNING "No Python interpreter found. Unable to generate C bindings for Iron.")
else()
    set(BINDINGS_DIR ${CMAKE_CURRENT_BINARY_DIR}/bindings)
    file(MAKE_DIRECTORY ${BINDINGS_DIR})
    set(OPENCMISS_H ${BINDINGS_DIR}/opencmiss.h)
    
    # C-Bindings extra target
    set(HAVE_C_BINDINGS FALSE)
    if (WITH_C_BINDINGS)
        message(STATUS "Configuring to create C bindings")
        set(OPENCMISS_C_F90 ${BINDINGS_DIR}/opencmiss_c.f90)
        #set_source_files_properties(${OPENCMISS_C_F90} PROPERTIES GENERATED TRUE)
        add_custom_command(OUTPUT ${OPENCMISS_C_F90} ${OPENCMISS_H}
            COMMAND ${PYTHON_EXECUTABLE} generate_bindings ${CMAKE_CURRENT_SOURCE_DIR} C ${OPENCMISS_H} ${OPENCMISS_C_F90}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bindings)
        set(HAVE_C_BINDINGS TRUE)
    endif()
    
    # Python-Bindings extra target
    set(HAVE_Python_BINDINGS FALSE)
    if (WITH_Python_BINDINGS)
        message(STATUS "Configuring to create Python bindings")
        find_package(PythonLibs QUIET)
        find_package(SWIG QUIET)
        if (SWIG_FOUND AND PYTHONLIBS_FOUND AND HAVE_C_BINDINGS)
            # Copy interface files
            file(GLOB SWIG_INTERFACES ${CMAKE_CURRENT_SOURCE_DIR}/bindings/python/*.i)
            file(COPY ${SWIG_INTERFACES} DESTINATION ${BINDINGS_DIR})
            
            #Generate interface
            set(SWIG_IFACE ${BINDINGS_DIR}/opencmiss.i)
            add_custom_command(OUTPUT ${SWIG_IFACE}
                COMMAND ${PYTHON_EXECUTABLE}2 generate_bindings ${CMAKE_CURRENT_SOURCE_DIR} SWIG ${SWIG_IFACE}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bindings)
            
            #Generate module
            set(CMISS_PY ${BINDINGS_DIR}/CMISS.py)
            add_custom_command(OUTPUT ${CMISS_PY}
                COMMAND ${PYTHON_EXECUTABLE} generate_bindings ${CMAKE_CURRENT_BINARY_DIR} Python ${CMAKE_CURRENT_BINARY_DIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bindings)
            
            # Generate C wrapper
            set(PYTHON_WRAPPER ${BINDINGS_DIR}/python_wrapper.c)
            add_custom_command(OUTPUT ${PYTHON_WRAPPER}
                DEPENDS ${SWIG_IFACE}
                COMMAND ${SWIG_EXECUTABLE} -python -o ${PYTHON_WRAPPER}
                    -module opencmiss_swig -outdir . opencmiss_py.i
                WORKING_DIRECTORY ${BINDINGS_DIR})
            
            #Generate wrapper object
            execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/utils/numpy_include.py
                OUTPUT_VARIABLE NUMPY_INCLUDES
                RESULT_VARIABLE RES_NUMPY_INC
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            if (NOT RES_NUMPY_INC)
                add_library(iron_python SHARED ${PYTHON_WRAPPER} ${OPENCMISS_H})
                target_link_libraries(iron_python PUBLIC iron)
                target_include_directories(iron_python PRIVATE ${BINDINGS_DIR} ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
                install(TARGETS iron_python
                    EXPORT iron-config
                    DESTINATION lib
                    INCLUDES DESTINATION include/iron)
            else()
                message(FATAL_ERROR "boo")
            endif()
            
            set(HAVE_Python_BINDINGS TRUE)
        else()
            message(WARNING "No SWIG or Python libraries found or no C bindings built. Unable to generate Python bindings for Iron.")            
        endif()
    endif()   
endif()