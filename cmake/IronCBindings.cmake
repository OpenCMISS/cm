# C-Bindings extra target
find_package(PythonInterp QUIET)
set(HAVE_C_BINDINGS FALSE)
if (PYTHONINTERP_FOUND)
    set(OPENCMISS_H ${CMAKE_CURRENT_BINARY_DIR}/opencmiss.h)
    set(OPENCMISS_C_F90 ${CMAKE_CURRENT_BINARY_DIR}/opencmiss_c.f90)
    set_source_files_properties(${OPENCMISS_C_F90} PROPERTIES GENERATED TRUE)
    add_custom_target(cbindings 
        COMMAND ${PYTHON_EXECUTABLE} generate_bindings ${CMAKE_CURRENT_SOURCE_DIR} C ${OPENCMISS_H} ${OPENCMISS_C_F90}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bindings)
    set(HAVE_C_BINDINGS TRUE)
else()
    message(WARNING "No Python interpreter found. Unable to generate C bindings for Iron.")
endif()