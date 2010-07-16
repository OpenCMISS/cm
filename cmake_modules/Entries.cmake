IF(${OPERATING_SYSTEM} MATCHES linux)# Linux
  SET(MACHINE_ENTRY machine_constants_linux)
ELSEIF(${OPERATING_SYSTEM} MATCHES aix)#AIX
    SET(MACHINE_ENTRY machine_constants_aix)
ELSE(${OPERATING_SYSTEM} MATCHES linux)# windows
    SET(MACHINE_ENTRY machine_constants_win32)
ENDIF(${OPERATING_SYSTEM} MATCHES linux)

IF(${USECELLML} MATCHES true)
     SET(CELLML_ENTRY  cmiss_cellml)
ELSE(${USECELLML} MATCHES true)
     SET(CELLML_ENTRY  cmiss_cellml_dummy)
ENDIF(${USECELLML} MATCHES true)



SET(OPENCMISS_ENTRIES 
  advection_diffusion_equation_routines 
  analytic_analysis_routines
  base_routines
  basis_routines
  bioelectric_routines
  biodomain_equation_routines
  boundary_condition_routines
  blas
  classical_field_routines
  cmiss
  cmiss_mpi
  cmiss_parmetis
  cmiss_petsc
  cmiss_petsc_types
  computational_environment
  constants
  control_loop_routines
  coordinate_routines
  Darcy_equations_routines
  data_point_routines
  data_projection_routines
  diffusion_advection_diffusion_routines
  diffusion_diffusion_routines
  diffusion_equation_routines
  distributed_matrix_vector
  distributed_matrix_vector_IO
  domain_mappings
  elasticity_routines
  electromechanics_routines
  equations_routines
  equations_mapping_routines
  equations_matrices_routines
  equations_set_constants
  equations_set_routines
  field_routines
  field_IO_routines
  finite_elasticity_Darcy_routines
  finite_elasticity_routines
  fluid_mechanics_routines
  fluid_mechanics_IO_routines
  FieldExport.c
  Galerkin_projection_routines
  generated_mesh_routines
  Helmholtz_equations_routines
  history_routines
  input_output
  interface_routines
  interface_conditions_constants
  interface_conditions_routines
  interface_equations_routines
  interface_mapping_routines
  interface_matrices_routines
  iso_varying_string
  kinds
  Laplace_equations_routines
  linear_elasticity_routines
  lists
  maths
  matrix_vector
  mesh_routines
  multi_compartment_transport_routines
  multi_physics_routines
  Navier_Stokes_equations_routines
  node_routines
  opencmiss
  opencmiss_c
  Poisson_equations_routines
  problem_constants
  problem_routines
  region_routines
  Stokes_equations_routines
  solver_routines
  solver_mapping_routines
  solver_matrices_routines
  sorting
  Stokes_equations_routines
  strings
  test_framework_routines
  timer_c.c
  timer_f
  trees
  types
  ${MACHINE_ENTRY} 
  ${CELLML_ENTRY}
)