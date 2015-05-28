set(IRON_C_SRC
    #binary_file_c.c
    cmiss_c.c
    external_dae_solver_routines.c
    FieldExport.c
    timer_c.c
)
set(IRON_HEADERS
    external_dae_solver_routines.h
    FieldExport.h
    FieldExportConstants.h
)
set(IRON_Fortran_SRC
    advection_diffusion_equation_routines.f90
    analytic_analysis_routines.f90
    base_routines.f90
    basis_routines.f90
    #binary_file_f.f90
    biodomain_equation_routines.f90
    bioelectric_finite_elasticity_routines.f90
    bioelectric_routines.f90
    blas.f90
    boundary_condition_routines.f90
    Burgers_equation_routines.f90
    characteristic_equation_routines.f90
    classical_field_routines.f90
    cmiss_fortran_c.f90
    cmiss_cellml.f90
    cmiss_mpi.f90
    cmiss_parmetis.f90
    cmiss_petsc_types.f90
    cmiss_petsc.f90
    cmiss.f90
    computational_environment.f90
    constants.f90
    control_loop_routines.f90
    coordinate_routines.f90
    Darcy_equations_routines.f90
    Darcy_pressure_equations_routines.f90
    data_point_routines.f90
    data_projection_routines.f90
    diffusion_advection_diffusion_routines.f90
    diffusion_diffusion_routines.f90
    diffusion_equation_routines.f90
    distributed_matrix_vector_IO.f90
    distributed_matrix_vector.f90
    domain_mappings.f90
    elasticity_routines.f90
    electromechanics_routines.f90
    electrophysiology_cell_routines.f90
    equations_mapping_routines.f90
    equations_matrices_routines.f90
    equations_routines.f90
    equations_set_constants.f90
    equations_set_routines.f90
    field_IO_routines.f90
    field_routines.f90
    finite_elasticity_Darcy_routines.f90
    finite_elasticity_fluid_pressure_routines.f90
    finite_elasticity_routines.f90
    fitting_routines.f90
    fluid_mechanics_IO_routines.f90
    fluid_mechanics_routines.f90
    fsi_routines.f90
    generated_mesh_routines.f90
    Hamilton_Jacobi_equations_routines.f90
    Helmholtz_equations_routines.f90
    #Helmholtz_TEMPLATE_equations_routines.f90
    history_routines.f90
    input_output.f90
    interface_conditions_constants.f90
    interface_conditions_routines.f90
    interface_equations_routines.f90
    interface_mapping_routines.f90
    interface_matrices_constants.f90
    interface_matrices_routines.f90
    interface_operators_routines.f90
    interface_routines.f90
    iso_varying_string.f90
    kinds.f90
    lapack.f90
    Laplace_equations_routines.f90
    linear_elasticity_routines.f90
    linkedlist_routines.f90
    lists.f90
    maths.f90
    matrix_vector.f90
    mesh_routines.f90
    monodomain_equations_routines.f90
    multi_compartment_transport_routines.f90
    multi_physics_routines.f90
    Navier_Stokes_equations_routines.f90
    node_routines.f90
    opencmiss.f90
    Poiseuille_equations_routines.f90
    Poisson_equations_routines.f90
    problem_constants.f90
    problem_routines.f90
    reaction_diffusion_equation_routines.f90
    reaction_diffusion_IO_routines.f90
    region_routines.f90
    solver_mapping_routines.f90
    solver_matrices_routines.f90
    solver_routines.f90
    sorting.f90
    Stokes_equations_routines.f90
    strings.f90
    test_framework_routines.f90
    timer_f.f90
    trees.f90
    types.f90
    util_array.f90
)
# Add platform dependent files
IF(${OPERATING_SYSTEM} MATCHES linux)
    list(APPEND IRON_Fortran_SRC machine_constants_linux.f90)
    #list(INSERT IRON_Fortran_SRC 0 machine_constants_linux.f90)
ELSEIF(${OPERATING_SYSTEM} MATCHES darwin)
    list(APPEND IRON_Fortran_SRC machine_constants_linux.f90)
ELSEIF(${OPERATING_SYSTEM} MATCHES aix)
    list(APPEND IRON_Fortran_SRC machine_constants_aix.f90)
ELSE(${OPERATING_SYSTEM} MATCHES linux)
    list(APPEND IRON_Fortran_SRC machine_constants_windows.f90)
ENDIF()
# Add optional files for interop
if (WITH_FIELDML)
    list(APPEND IRON_Fortran_SRC
        fieldml_input_routines.f90
        fieldml_output_routines.f90
        fieldml_types.f90
        fieldml_util_routines.f90
    )
endif()

# Fix paths to files
set(FIXPATH_VARS IRON_C_SRC IRON_Fortran_SRC)
foreach(varname ${FIXPATH_VARS})
    set(_TMP )
    foreach(filename ${${varname}})
        list(APPEND _TMP ${CMAKE_CURRENT_SOURCE_DIR}/src/${filename}) 
    endforeach()
    set(${varname} ${_TMP})
endforeach()
# Set combined sources variable
set(IRON_SRC ${IRON_C_SRC} ${IRON_Fortran_SRC})