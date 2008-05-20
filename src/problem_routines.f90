!> \file
!> $Id: problem_routines.f90 28 2007-07-27 08:35:14Z cpb $
!> \author Chris Bradley
!> \brief This module handles all problem routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module handles all problem routines.
MODULE PROBLEM_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE FIELD_ROUTINES
  USE GLOBAL_MATRICES_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE MPI
  USE SOLUTION_MAPPING_ROUTINES
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Problem Classes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_CLASS=0
  
  INTEGER(INTG), PARAMETER :: PROBLEM_ELASTICITY_CLASS=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FLUID_MECHANICS_CLASS=2
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROMAGNETICS_CLASS=3
  INTEGER(INTG), PARAMETER :: PROBLEM_CLASSICAL_FIELD_CLASS=4
  
  INTEGER(INTG), PARAMETER :: PROBLEM_MODAL_CLASS=5
  INTEGER(INTG), PARAMETER :: PROBLEM_FITTING_CLASS=6
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISATION_CLASS=7

  !Problem types
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_TYPE=0
  !Elasticity class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTICITY_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_TYPE=2
  !Fluid mechanics class
  INTEGER(INTG), PARAMETER :: PROBLEM_STOKES_FLUID_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NAVIER_STOKES_FLUID_TYPE=2
  !Electromagnetics class
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROSTATIC_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_MAGNETOSTATIC_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MAXWELLS_EQUATIONS_TYPE=3
  !Classical field class
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_POISSON_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_HELMHOLTZ_EQUATION_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_WAVE_EQUATION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_EQUATION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_BIHARMONIC_EQUATION_TYPE=8
  !Modal class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTIC_MODAL_TYPE=1
  
  !Problem subtypes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SUBTYPE=0
  !Elasticity class
  !Fluid mechanics class
  !Electromagnetics class
  !Classical field class
  !  Laplace equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_LAPLACE_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_LAPLACE_SUBTYPE=2
  !  Poisson equation
  !  Helmholtz equation
  !  Wave equation
  !  Diffusion equation
  !Modal class

  !> \addtogroup PROBLEM_ROUTINES_LinearityTypes PROBLEM_ROUTINES::LinearityTypes
  !> \brief The linearity type parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_PROBLEM_LINEARITIES=3 !<The number of problem linearity types defined. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR=1 !<The problemis linear. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR=2 !<The problem is non-linear. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_BCS=3 !<The problem has non-linear boundary conditions. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  !>@}
  
  !> \addtogroup PROBLEM_ROUTINES_TimeDepedenceTypes PROBLEM_ROUTINES::TimeDepedenceTypes
  !> \brief The time dependence type parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_PROBLEM_TIME_TYPES=3 !<The number of problem time dependence types defined. \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC=1 !<The problem is static and has no time dependence \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC=2 !<The problem is dynamic \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC=3 !<The problem is quasi-static \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  !>@}
  
  !> \addtogroup PROBLEM_ROUTINES_SolutionMethods PROBLEM_ROUTINES::SolutionMethods
  !> \brief The solution method parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_PROBLEM_SOLUTION_METHODS=6 !<The number of solution methods defined. \see PROBLEM_ROUTINES_SolutionMethods,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_FEM_SOLUTION_METHOD=1 !<Finite Element Method solution method \see PROBLEM_ROUTINES_SolutionMethods,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_BEM_SOLUTION_METHOD=2 !<Boundary Element Method solution method \see PROBLEM_ROUTINES_SolutionMethods,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_FD_SOLUTION_METHOD=3 !<Finite Difference solution method \see PROBLEM_ROUTINES_SolutionMethods,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_FV_SOLUTION_METHOD=4 !<Finite Volume solution method \see PROBLEM_ROUTINES_SolutionMethods,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_GFEM_SOLUTION_METHOD=5 !<Grid-based Finite Element Method solution method \see PROBLEM_ROUTINES_SolutionMethods,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_GFV_SOLUTION_METHOD=6 !<Grid-based Finite Volume solution method \see PROBLEM_ROUTINES_SolutionMethods,PROBLEM_ROUTINES
  !>@}

  !> \addtogroup PROBLEM_ROUTINES_FixedConditions PROBLEM_ROUTINES::FixedConditions
  !> \brief Fixed conditons parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_NOT_FIXED=0 !<The dof is not fixed. \see PROBLEM_ROUTINES_FixedConditions,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_FIXED_BOUNDARY_CONDITION=1 !<The dof is fixed as a boundary condition. \see PROBLEM_ROUTINES_FixedConditions,PROBLEM_ROUTINES
  !>@}
  
  !> \addtogroup PROBLEM_ROUTINES_SetupTypes PROBLEM_ROUTINES::SetupTypes
  !> \brief Setup type parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_INITIAL_SETUP_TYPE=1 !<Initial setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_GEOMETRY_SETUP_TYPE=2 !<Geometry setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_DEPENDENT_SETUP_TYPE=3 !<Dependent variables setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_MATERIALS_SETUP_TYPE=4 !<Materials setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SOURCE_SETUP_TYPE=5 !<Source setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_ANALYTIC_SETUP_TYPE=6 !<Analytic setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_FIXED_CONDITIONS_SETUP_TYPE=7 !<Fixed conditions setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_SETUP_TYPE=8 !<Solver setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_SETUP_TYPE=9 !<Solution parameters setup for a problem. \see PROBLEM_ROUTINES_SetupTypes,PROBLEM_ROUTINES
  !>@}
  
  !> \addtogroup PROBLEM_ROUTINES_SetupActionTypes PROBLEM_ROUTINES::SetupActionTypes
  !> \brief Setup action type parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_START_ACTION=1 !<Start setup action. \see PROBLEM_ROUTINES_SetupActionTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_FINISH_ACTION=2 !<Finish setup action. \see PROBLEM_ROUTINES_SetupActionTypes,PROBLEM_ROUTINES
  !>@}

  !> \addtogroup PROBLEM_ROUTINES_GlobalMatrixStructureTypes PROBLEM_ROUTINES::GlobalMatrixStructureTypes
  !> \brief Global matrices structure (sparsity) types
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_GLOBAL_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see PROBLEM_ROUTINES_GlobalMatrixStructureTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_GLOBAL_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see PROBLEM_ROUTINES_GlobalMatrixStructureTypes,PROBLEM_ROUTINES
 !>@}

  !> \addtogroup PROBLEM_ROUTINES_SolutionOutputTypes PROBLEM_ROUTINES::SolutionOutputTypes
  !> \brief Solution output types
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_NO_OUTPUT=0 !<No output \see PROBLEM_ROUTINES_SolutionOutputTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_TIMING_OUTPUT=1 !<Timing information output \see PROBLEM_ROUTINES_SolutionOutputTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_GLOBAL_MATRIX_OUTPUT=2 !<All below and Global matrices output \see PROBLEM_ROUTINES_SolutionOutputTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_ELEMENT_MATRIX_OUTPUT=3 !<All below and Element matrices output \see PROBLEM_ROUTINES_SolutionOutputTypes,PROBLEM_ROUTINES
  !>@}

  !> \addtogroup PROBLEM_ROUTINES_SolutionGlobalSparsityTypes PROBLEM_ROUTINES::SolutionGlobalSparsityTypes
  !> \brief Solution output types
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_SPARSE_GLOBAL_MATRICES=1 !<Use sparse global matrices \see PROBLEM_ROUTINES_SolutionGlobalSparsityTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_FULL_GLOBAL_MATRICES=2 !<Use fully populated global matrices \see PROBLEM_ROUTINES_SolutionGlobalSparsityTypes,PROBLEM_ROUTINES
 !>@}

  !Module types

  !Module variables

  !Interfaces

  INTERFACE PROBLEM_SPECIFICATION_SET
    MODULE PROCEDURE PROBLEM_SPECIFICATION_SET_NUMBER
    MODULE PROCEDURE PROBLEM_SPECIFICATION_SET_PTR
  END INTERFACE !PROBLEM_TYPE_SET

  INTERFACE PROBLEM_FIXED_CONDITIONS_SET_DOF
    MODULE PROCEDURE PROBLEM_FIXED_CONDITIONS_SET_DOFS
    MODULE PROCEDURE PROBLEM_FIXED_CONDITIONS_SET_DOF1
  END INTERFACE !PROBLEM_FIXED_CONDITIONS_SET_DOF

  PUBLIC PROBLEM_NO_CLASS,PROBLEM_ELASTICITY_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
    & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_MODAL_CLASS

  PUBLIC PROBLEM_NO_TYPE,PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE,PROBLEM_STOKES_FLUID_TYPE, &
    & PROBLEM_NAVIER_STOKES_FLUID_TYPE,PROBLEM_ELECTROSTATIC_TYPE,PROBLEM_MAGNETOSTATIC_TYPE, &
    & PROBLEM_MAXWELLS_EQUATIONS_TYPE,PROBLEM_LAPLACE_EQUATION_TYPE,PROBLEM_POISSON_EQUATION_TYPE, &
    & PROBLEM_HELMHOLTZ_EQUATION_TYPE,PROBLEM_WAVE_EQUATION_TYPE,PROBLEM_DIFFUSION_EQUATION_TYPE, &
    & PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE,PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,PROBLEM_BIHARMONIC_EQUATION_TYPE, &
    & PROBLEM_LINEAR_ELASTIC_MODAL_TYPE

  PUBLIC PROBLEM_NO_SUBTYPE,PROBLEM_STANDARD_LAPLACE_SUBTYPE,PROBLEM_GENERALISED_LAPLACE_SUBTYPE

  PUBLIC PROBLEM_LINEAR,PROBLEM_NONLINEAR,PROBLEM_NONLINEAR_BCS

  PUBLIC PROBLEM_STATIC,PROBLEM_DYNAMIC,PROBLEM_QUASISTATIC

  PUBLIC PROBLEM_FEM_SOLUTION_METHOD,PROBLEM_BEM_SOLUTION_METHOD,PROBLEM_FD_SOLUTION_METHOD,PROBLEM_FV_SOLUTION_METHOD, &
    & PROBLEM_GFEM_SOLUTION_METHOD,PROBLEM_GFV_SOLUTION_METHOD

  PUBLIC PROBLEM_NOT_FIXED,PROBLEM_FIXED_BOUNDARY_CONDITION

  PUBLIC PROBLEM_SOLUTION_NO_OUTPUT,PROBLEM_SOLUTION_TIMING_OUTPUT,PROBLEM_SOLUTION_GLOBAL_MATRIX_OUTPUT, &
    & PROBLEM_SOLUTION_ELEMENT_MATRIX_OUTPUT

  PUBLIC PROBLEM_SOLUTION_SPARSE_GLOBAL_MATRICES,PROBLEM_SOLUTION_FULL_GLOBAL_MATRICES

  PUBLIC PROBLEM_CREATE_START,PROBLEM_CREATE_FINISH,PROBLEM_DESTROY,PROBLEMS_INITIALISE,PROBLEMS_FINALISE

  PUBLIC PROBLEM_FIXED_CONDITIONS_CREATE_START,PROBLEM_FIXED_CONDITIONS_CREATE_FINISH,PROBLEM_FIXED_CONDITIONS_SET_DOF

  PUBLIC PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET,PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET, &
    & PROBLEM_MATERIALS_CREATE_START,PROBLEM_MATERIALS_CREATE_FINISH,PROBLEM_MATERIALS_SCALING_SET

  PUBLIC PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET,PROBLEM_DEPENDENT_CREATE_START,PROBLEM_DEPENDENT_CREATE_FINISH, &
    & PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET,PROBLEM_DEPENDENT_SCALING_SET
  
  PUBLIC PROBLEM_SOLUTION_CREATE_START,PROBLEM_SOLUTION_CREATE_FINISH,PROBLEM_SOLUTION_GLOBAL_SPARSITY_TYPE_SET, &
    & PROBLEM_SOLUTION_OUTPUT_TYPE_SET
  
  PUBLIC PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET,PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET, &
    & PROBLEM_SOURCE_CREATE_START,PROBLEM_SOURCE_CREATE_FINISH,PROBLEM_SOURCE_SCALING_SET
  
  PUBLIC PROBLEM_SPECIFICATION_SET

  PUBLIC PROBLEM_SOLVE

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_ANALYTIC_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_ANALYTIC_CREATE_FINISH
    !###  Description:
    !###    Finish the creation of a analytic solution for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_ANALYTIC_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%ANALYTIC)) THEN
        !Finish the problem specific source setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_ANALYTIC_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
        !Finish the source creation
        PROBLEM%ANALYTIC%ANALYTIC_FINISHED=.FALSE.
      ELSE
        CALL FLAG_ERROR("The problem analytic is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_ANALYTIC_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_ANALYTIC_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_ANALYTIC_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_ANALYTIC_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_ANALYTIC_CREATE_START(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_ANALYTIC_CREATE_START
    !###  Description:
    !###    Start the creation of a analytic solution for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_ANALYTIC_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%ANALYTIC)) THEN
        CALL FLAG_ERROR("The problem analytic is already associated",ERR,ERROR,*999)        
      ELSE
        !Initialise the problem analytic
        CALL PROBLEM_ANALYTIC_INITIALISE(PROBLEM,ERR,ERROR,*999)
        !Start the problem specific analytic setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_ANALYTIC_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_ANALYTIC_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_ANALYTIC_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_ANALYTIC_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_ANALYTIC_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_ANALYTIC_FINALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_ANALYTIC_FINALISE
    !###  Description:
    !###    Finalise the analytic solution for a problem and deallocate all memory.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_ANALYTIC_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%ANALYTIC)) THEN        
        DEALLOCATE(PROBLEM%ANALYTIC)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_ANALYTIC_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_ANALYTIC_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_ANALYTIC_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_ANALYTIC_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_ANALYTIC_INITIALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_ANALYTIC_INITIALISE
    !###  Description:
    !###    Initialises the analytic solution for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_ANALYTIC_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%ANALYTIC)) THEN
        CALL FLAG_ERROR("Analytic is already associated for this problem",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM%ANALYTIC,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem analytic",ERR,ERROR,*999)
        PROBLEM%ANALYTIC%PROBLEM=>PROBLEM
        PROBLEM%ANALYTIC%ANALYTIC_FINISHED=.FALSE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_ANALYTIC_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_ANALYTIC_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_ANALYTIC_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_ANALYTIC_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Assembles the global stiffness matrix and rhs for a linear static problem using the finite element method.
  SUBROUTINE PROBLEM_ASSEMBLE_LINEAR_STATIC_FEM(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    
    CALL ENTERS("PROBLEM_ASSEMBLE_LINEAR_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      DEPENDENT_FIELD=>PROBLEM%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        PROBLEM_SOLUTION=>PROBLEM%SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          GLOBAL_MATRICES=>PROBLEM%SOLUTION%GLOBAL_MATRICES
          IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
            IF(PROBLEM_SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Start the transfer of the solution values that have been set as part of the boundary conditions
            CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Problem interpolation setup
            CALL PROBLEM_INTERPOLATION_INITIALISE(PROBLEM%SOLUTION,ERR,ERROR,*999)
            !Initialise the matrices and rhs vector
            CALL GLOBAL_MATRICES_VALUES_INITIALISE(GLOBAL_MATRICES,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL GLOBAL_MATRICES_ELEMENT_INITIALISE(GLOBAL_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(PROBLEM_SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for Global setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for Global setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL
              ne=ELEMENTS_MAPPING%INTERNAL_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL GLOBAL_MATRICES_ELEMENT_CALCULATE(GLOBAL_MATRICES,ne,ERR,ERROR,*999)
              CALL PROBLEM_FINITE_ELEMENT_CALCULATE(PROBLEM,ne,ERR,ERROR,*999)
              CALL GLOBAL_MATRICES_ELEMENT_ADD(GLOBAL_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(PROBLEM_SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal Global assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal Global assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
            !Finish the transfer of the solution values.
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Output timing information if required
            IF(PROBLEM_SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
!!TODO: sort out combining boundary and ghost list
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY
              ne=ELEMENTS_MAPPING%BOUNDARY_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL GLOBAL_MATRICES_ELEMENT_CALCULATE(GLOBAL_MATRICES,ne,ERR,ERROR,*999)
              CALL PROBLEM_FINITE_ELEMENT_CALCULATE(PROBLEM,ne,ERR,ERROR,*999)
              CALL GLOBAL_MATRICES_ELEMENT_ADD(GLOBAL_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_GHOST
              ne=ELEMENTS_MAPPING%GHOST_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL GLOBAL_MATRICES_ELEMENT_CALCULATE(GLOBAL_MATRICES,ne,ERR,ERROR,*999)
              CALL PROBLEM_FINITE_ELEMENT_CALCULATE(PROBLEM,ne,ERR,ERROR,*999)
              CALL GLOBAL_MATRICES_ELEMENT_ADD(GLOBAL_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx          
            !Output timing information if required
            IF(PROBLEM_SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost Global assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost Global assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for Global assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for Global assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL GLOBAL_MATRICES_ELEMENT_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*999)
            !Finalise the problem interpolation
            CALL PROBLEM_INTERPOLATION_FINALISE(PROBLEM%SOLUTION,ERR,ERROR,*999)
            !Output global matrices and RHS vector if required
            IF(PROBLEM_SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_GLOBAL_MATRIX_OUTPUT) THEN
              CALL GLOBAL_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,GLOBAL_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(PROBLEM_SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for Global assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for Global assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_ASSEMBLE_LINEAR_STATIC_FEM")
    RETURN
999 CALL ERRORS("PROBLEM_ASSEMBLE_LINEAR_STATIC_FEM",ERR,ERROR)
    CALL EXITS("PROBLEM_ASSEMBLE_LINEAR_STATIC_FEM")
    RETURN 1
  END SUBROUTINE PROBLEM_ASSEMBLE_LINEAR_STATIC_FEM

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a clasical field class finite element problem.
  SUBROUTINE PROBLEM_CLASSICAL_FIELD_CLASS_FINITE_ELEMENT_CALCULATE(PROBLEM,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_CLASSICAL_FIELD_CLASS_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%TYPE)
      CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
        CALL PROBLEM_LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE(PROBLEM,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(PROBLEM_POISSON_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_WAVE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a classical field problem class"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CLASSICAL_FIELD_CLASS_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("PROBLEM_CLASSICAL_FIELD_CLASS_FINITE_ELEMENT_CALCUALTE",ERR,ERROR)
    CALL EXITS("PROBLEM_CLASSICAL_FIELD_CLASS_FINITE_ELEMENT_CALCUALTE")
    RETURN 1
  END SUBROUTINE PROBLEM_CLASSICAL_FIELD_CLASS_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_CLASSICAL_FIELD_CLASS_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_TYPE_CLASSICAL_FIELD_CLASS_SETUP
    !###  Description:
    !###    Sets up the solution for a classical field problem class.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_CLASSICAL_FIELD_CLASS_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%TYPE)
      CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
        CALL PROBLEM_LAPLACE_EQUATION_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(PROBLEM_POISSON_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_WAVE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a classical field problem class"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CLASSICAL_FIELD_CLASS_SETUP")
    RETURN
999 CALL ERRORS("PROBLEM_CLASSICAL_FIELD_CLASS_SETUP",ERR,ERROR)
    CALL EXITS("PROBLEM_CLASSICAL_FIELD_CLASS_SETUP")
    RETURN 1
  END SUBROUTINE PROBLEM_CLASSICAL_FIELD_CLASS_SETUP

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_CLASSICAL_FIELD_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_TYPE_CLASSICAL_FIELD_CLASS_SET
    !###  Description:
    !###    Sets/changes the problem type and subtype for a classical field problem class.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_TYPE_CLASSICAL_FIELD_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_EQUATION_TYPE)
      CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
        CALL PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*999)
      CASE(PROBLEM_POISSON_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_WAVE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM_EQUATION_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a classical field problem class"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CLASSICAL_FIELD_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("PROBLEM_CLASSICAL_FIELD_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_CLASSICAL_FIELD_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_CLASSICAL_FIELD_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_CREATE_FINISH(REGION,PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_CREATE_FINISH
    !###  Description:
    !###    Finishes the process of creating a problem on a region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(PROBLEM)) THEN
        IF(ASSOCIATED(PROBLEM%REGION)) THEN
          IF(PROBLEM%REGION%USER_NUMBER==REGION%USER_NUMBER) THEN            
            !Finish the problem specific setup
            CALL PROBLEM_SETUP(PROBLEM,PROBLEM_INITIAL_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
            !Finish the problem specific geometry setup
            CALL PROBLEM_SETUP(PROBLEM,PROBLEM_GEOMETRY_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
            !Finish the problem creation
            PROBLEM%PROBLEM_FINISHED=.TRUE.
          ELSE
            LOCAL_ERROR="The region user number of the specified problem ("// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%REGION%USER_NUMBER,"*",ERR,ERROR))// &
              & ") does not match the user number of the specified region ("// &
              & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//")"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem region is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE        
        CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Region = ",REGION%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of problems = ",REGION%PROBLEMS%NUMBER_OF_PROBLEMS,ERR,ERROR,*999)
      DO problem_idx=1,REGION%PROBLEMS%NUMBER_OF_PROBLEMS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Problem number = ",problem_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number   = ", &
          & REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number = ", &
          & REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER,ERR,ERROR,*999)
      ENDDO !problem_idx    
    ENDIF
    
    CALL EXITS("PROBLEM_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("PROBLEM_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE PROBLEM_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_CREATE_START(USER_NUMBER,REGION,GEOM_FIBRE_FIELD,PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_CREATE_START
    !###  Description:
    !###    Starts the process of creating a problem defined by USER_NUMBER in the region identified by REGION.
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(FIELD_TYPE), POINTER :: GEOM_FIBRE_FIELD
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx
    TYPE(PROBLEM_TYPE), POINTER :: NEW_PROBLEM
    TYPE(PROBLEM_PTR_TYPE), POINTER :: NEW_PROBLEMS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    NULLIFY(NEW_PROBLEM)
    NULLIFY(NEW_PROBLEMS)

    CALL ENTERS("PROBLEM_CREATE_START",ERR,ERROR,*999)

    NULLIFY(PROBLEM)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%PROBLEMS)) THEN
        CALL PROBLEM_USER_NUMBER_FIND(USER_NUMBER,REGION,PROBLEM,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM)) THEN
          LOCAL_ERROR="Problem number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
            IF(GEOM_FIBRE_FIELD%FIELD_FINISHED) THEN
              IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR.GEOM_FIBRE_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
                !Allocate the new problem
                ALLOCATE(NEW_PROBLEM,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problem",ERR,ERROR,*999)
                !Initalise problem
                CALL PROBLEM_INITIALISE(NEW_PROBLEM,ERR,ERROR,*999)
                !Set default problem values
                NEW_PROBLEM%USER_NUMBER=USER_NUMBER
                NEW_PROBLEM%GLOBAL_NUMBER=REGION%PROBLEMS%NUMBER_OF_PROBLEMS+1
                NEW_PROBLEM%PROBLEMS=>REGION%PROBLEMS
                NEW_PROBLEM%REGION=>REGION
                !Default to a standardised Laplace.
                NEW_PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
                NEW_PROBLEM%TYPE=PROBLEM_LAPLACE_EQUATION_TYPE
                NEW_PROBLEM%SUBTYPE=PROBLEM_STANDARD_LAPLACE_SUBTYPE
                !Start problem specific setup
                CALL PROBLEM_SETUP(NEW_PROBLEM,PROBLEM_INITIAL_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)
                !Set up the geometry
                CALL PROBLEM_GEOMETRY_INITIALISE(NEW_PROBLEM,ERR,ERROR,*999)
                IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
                  NEW_PROBLEM%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD
                  NULLIFY(NEW_PROBLEM%GEOMETRY%FIBRE_FIELD)
                ELSE
                  NEW_PROBLEM%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD%GEOMETRIC_FIELD
                  NEW_PROBLEM%GEOMETRY%FIBRE_FIELD=>GEOM_FIBRE_FIELD
                ENDIF
                !Set up problem specific geometry
                CALL PROBLEM_SETUP(NEW_PROBLEM,PROBLEM_GEOMETRY_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)                
                
                !Add new problem into list of problems in the region
                ALLOCATE(NEW_PROBLEMS(REGION%PROBLEMS%NUMBER_OF_PROBLEMS+1),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problems",ERR,ERROR,*999)
                DO problem_idx=1,REGION%PROBLEMS%NUMBER_OF_PROBLEMS
                  NEW_PROBLEMS(problem_idx)%PTR=>REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR
                ENDDO !problem_idx
                NEW_PROBLEMS(REGION%PROBLEMS%NUMBER_OF_PROBLEMS+1)%PTR=>NEW_PROBLEM
                IF(ASSOCIATED(REGION%PROBLEMS%PROBLEMS)) DEALLOCATE(REGION%PROBLEMS%PROBLEMS)
                REGION%PROBLEMS%PROBLEMS=>NEW_PROBLEMS
                REGION%PROBLEMS%NUMBER_OF_PROBLEMS=REGION%PROBLEMS%NUMBER_OF_PROBLEMS+1
                PROBLEM=>NEW_PROBLEM
              ELSE
                CALL FLAG_ERROR("The specified geometric field is not a geometric or fibre field",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified geometric field is not finished",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The specified geometric field is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The problems on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated. Initialise problems first."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_CREATE_START")
    RETURN 1   
  END SUBROUTINE PROBLEM_CREATE_START
  
  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_DESTROY(USER_NUMBER,REGION,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_DESTROY
    !###  Description:
    !###    Destroys a problem identified by a user number on the give region and deallocates all memory.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx,problem_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_PTR_TYPE), POINTER :: NEW_PROBLEMS(:)

    NULLIFY(NEW_PROBLEMS)

    CALL ENTERS("PROBLEM_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%PROBLEMS)) THEN
        
        !Find the problem identified by the user number
        FOUND=.FALSE.
        problem_position=0
        DO WHILE(problem_position<REGION%PROBLEMS%NUMBER_OF_PROBLEMS.AND..NOT.FOUND)
          problem_position=problem_position+1
          IF(REGION%PROBLEMS%PROBLEMS(problem_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          PROBLEM=>REGION%PROBLEMS%PROBLEMS(problem_position)%PTR
          
          !Destroy all the problem components
          CALL PROBLEM_FINALISE(PROBLEM,ERR,ERROR,*999)
          
          !Remove the problem from the list of problems
          IF(REGION%PROBLEMS%NUMBER_OF_PROBLEMS>1) THEN
            ALLOCATE(NEW_PROBLEMS(REGION%PROBLEMS%NUMBER_OF_PROBLEMS-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problems",ERR,ERROR,*999)
            DO problem_idx=1,REGION%PROBLEMS%NUMBER_OF_PROBLEMS
              IF(problem_idx<problem_position) THEN
                NEW_PROBLEMS(problem_idx)%PTR=>REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR
              ELSE IF(problem_idx>problem_position) THEN
                REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER=REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER-1
                NEW_PROBLEMS(problem_idx-1)%PTR=>REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR
              ENDIF
            ENDDO !problem_idx
            DEALLOCATE(REGION%PROBLEMS%PROBLEMS)
            REGION%PROBLEMS%PROBLEMS=>NEW_PROBLEMS
            REGION%PROBLEMS%NUMBER_OF_PROBLEMS=REGION%PROBLEMS%NUMBER_OF_PROBLEMS-1
          ELSE
            DEALLOCATE(REGION%PROBLEMS%PROBLEMS)
            REGION%PROBLEMS%NUMBER_OF_PROBLEMS=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Problem number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The problems on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("PROBLEM_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_PROBLEMS)) DEALLOCATE(NEW_PROBLEMS)
998 CALL ERRORS("PROBLEM_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_DESTROY")
    RETURN 1   
  END SUBROUTINE PROBLEM_DESTROY
  
  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_FINALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_FINALISE
    !###  Description:
    !###    Finalise the problem and deallocate all memory.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_GEOMETRY_FINALISE(PROBLEM,ERR,ERROR,*999)
      CALL PROBLEM_DEPENDENT_FINALISE(PROBLEM,ERR,ERROR,*999)
      CALL PROBLEM_MATERIALS_FINALISE(PROBLEM,ERR,ERROR,*999)
      CALL PROBLEM_SOURCE_FINALISE(PROBLEM,ERR,ERROR,*999)
      CALL PROBLEM_ANALYTIC_FINALISE(PROBLEM,ERR,ERROR,*999)
      CALL PROBLEM_FIXED_CONDITIONS_FINALISE(PROBLEM,ERR,ERROR,*999)
      CALL PROBLEM_SOLUTION_FINALISE(PROBLEM,ERR,ERROR,*999)
      DEALLOCATE(PROBLEM)
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_FINALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a finite element problem.
  SUBROUTINE PROBLEM_FINITE_ELEMENT_CALCULATE(PROBLEM,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(ELEMENT_VECTOR_TYPE), POINTER :: ELEMENT_VECTOR
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%CLASS)
      CASE(PROBLEM_ELASTICITY_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_FLUID_MECHANICS_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
        CALL PROBLEM_CLASSICAL_FIELD_CLASS_FINITE_ELEMENT_CALCULATE(PROBLEM,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(PROBLEM_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM%CLASS,"*",ERR,ERROR))//" is not valid"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT      
      IF(PROBLEM%SOLUTION%OUTPUT_TYPE>=PROBLEM_SOLUTION_ELEMENT_MATRIX_OUTPUT) THEN
        GLOBAL_MATRICES=>PROBLEM%SOLUTION%GLOBAL_MATRICES
        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element stiffness matrices:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",GLOBAL_MATRICES%NUMBER_OF_MATRICES, &
          & ERR,ERROR,*999)
        DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",GLOBAL_MATRICES%MATRICES(matrix_idx)%UPDATE_MATRIX, &
            & ERR,ERROR,*999)
          ELEMENT_MATRIX=>GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
            & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%COLUMN_DOFS, &
            & '("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
          CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
            & 8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
            & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))','(16X,8(X,E13.6))', &
            & ERR,ERROR,*999)        
        ENDDO !matrix_idx
        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element RHS vector :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",GLOBAL_MATRICES%UPDATE_VECTOR,ERR,ERROR,*999)
        ELEMENT_VECTOR=>GLOBAL_MATRICES%ELEMENT_VECTOR
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
          & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
          & '("  Vector(:)    :",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)      
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF

       
    CALL EXITS("PROBLEM_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("PROBLEM_FINITE_ELEMENT_CALCUALTE",ERR,ERROR)
    CALL EXITS("PROBLEM_FINITE_ELEMENT_CALCUALTE")
    RETURN 1
  END SUBROUTINE PROBLEM_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !
  
  !> Finalises the interpolation information for a problem and deallocates all memory
  SUBROUTINE PROBLEM_INTERPOLATION_FINALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_INTERPOLATION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%INTERPOLATION)) THEN
        CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_INTERP_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_INTERP_POINT,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_POINT,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_INTERP_POINT,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_INTERP_POINT,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS, &
          & ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_POINT_METRICS,ERR,ERROR,*999)
        DEALLOCATE(PROBLEM_SOLUTION%INTERPOLATION)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_INTERPOLATION_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_INTERPOLATION_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_INTERPOLATION_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_INTERPOLATION_FINALISE

  !
  !================================================================================================================================
  !

  !> Initialises the interpolation information for a problem solution
  SUBROUTINE PROBLEM_INTERPOLATION_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<The pointer to the problem solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_INTERPOLATION_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%INTERPOLATION)) THEN
        CALL FLAG_ERROR("Interpolation is already associated for this problem solution",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLUTION%INTERPOLATION,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution interpolation",ERR,ERROR,*999)
        PROBLEM_SOLUTION%INTERPOLATION%PROBLEM_SOLUTION=>PROBLEM_SOLUTION
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_PARAMETERS)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_INTERP_PARAMETERS)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_INTERP_PARAMETERS)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_POINT)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_INTERP_POINT)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_INTERP_POINT)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_INTERP_POINT)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS)
        NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_POINT_METRICS)
        
        PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_FIELD=>PROBLEM_SOLUTION%PROBLEM%GEOMETRY%GEOMETRIC_FIELD
        PROBLEM_SOLUTION%INTERPOLATION%FIBRE_FIELD=>PROBLEM_SOLUTION%PROBLEM%GEOMETRY%FIBRE_FIELD
        PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_FIELD=>PROBLEM_SOLUTION%PROBLEM%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(PROBLEM_SOLUTION%PROBLEM%MATERIALS)) THEN
          PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_FIELD=>PROBLEM_SOLUTION%PROBLEM%MATERIALS%MATERIAL_FIELD
        ELSE
          NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_FIELD)
        ENDIF
        IF(ASSOCIATED(PROBLEM_SOLUTION%PROBLEM%SOURCE)) THEN
          PROBLEM_SOLUTION%INTERPOLATION%SOURCE_FIELD=>PROBLEM_SOLUTION%PROBLEM%SOURCE%SOURCE_FIELD
        ELSE
          NULLIFY(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_FIELD)
        ENDIF
 
        CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_FIELD, &
          & FIELD_STANDARD_VARIABLE_TYPE,PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS, &
          & PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT, &
          & PROBLEM_SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_FIELD, &
          & FIELD_STANDARD_VARIABLE_TYPE,PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS, &
          & PROBLEM_SOLUTION%INTERPOLATION%DEPENDENT_INTERP_POINT,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_FIELD)) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_FIELD, &
            & FIELD_STANDARD_VARIABLE_TYPE,PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_PARAMETERS,  &
            &  PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_POINT,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_POINT,  &
            &  PROBLEM_SOLUTION%INTERPOLATION%FIBRE_INTERP_POINT_METRICS,ERR,ERROR,*999)
        ENDIF
        IF(ASSOCIATED(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_FIELD)) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_FIELD, &
            & FIELD_STANDARD_VARIABLE_TYPE,PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_INTERP_PARAMETERS,  &
            & PROBLEM_SOLUTION%INTERPOLATION%MATERIAL_INTERP_POINT,ERR,ERROR,*999)
        ENDIF
        IF(ASSOCIATED(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_FIELD)) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_FIELD, &
            & FIELD_STANDARD_VARIABLE_TYPE,PROBLEM_SOLUTION%INTERPOLATION%SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_INITIALISE(PROBLEM_SOLUTION%INTERPOLATION%SOURCE_INTERP_PARAMETERS, &
            & PROBLEM_SOLUTION%INTERPOLATION%SOURCE_INTERP_POINT,ERR,ERROR,*999)
        ENDIF
        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_INTERPOLATION_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_INTERPOLATION_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_INTERPOLATION_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_INTERPOLATION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and rhs for a laplace equation finite element problem
  SUBROUTINE PROBLEM_LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE(PROBLEM,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) ng,mh,mhs,mi,ms,nh,nhs,ni,ns
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)
        !Store all these in global matrices/somewhere else?????
        DEPENDENT_FIELD=>PROBLEM%SOLUTION%INTERPOLATION%DEPENDENT_FIELD
        GEOMETRIC_FIELD=>PROBLEM%SOLUTION%INTERPOLATION%GEOMETRIC_FIELD
        GLOBAL_MATRICES=>PROBLEM%SOLUTION%GLOBAL_MATRICES
        FIELD_VARIABLE=>GLOBAL_MATRICES%MATRICES(1)%VARIABLE
        DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,PROBLEM%SOLUTION%INTERPOLATION% &
          & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
        !Loop over gauss points
        DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
          CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,PROBLEM%SOLUTION%INTERPOLATION% &
            & GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,PROBLEM%SOLUTION%INTERPOLATION% &
            & GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
          !Calculate RWG.
!!TODO: Think about symmetric problems. 
          RWG=PROBLEM%SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN*QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
          !Loop over field components
          mhs=0          
          DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns.
            DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
              mhs=mhs+1
              nhs=0
              IF(GLOBAL_MATRICES%MATRICES(1)%UPDATE_MATRIX) THEN
                !Loop over element columns
                DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    nhs=nhs+1
                    DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                      PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                      PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                    ENDDO !ni
                    SUM=0.0_DP
                    DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        SUM=SUM+PGMSI(mi)*PGNSI(ni)*PROBLEM%SOLUTION%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%GU(mi,ni)
                      ENDDO !ni
                    ENDDO !mi
                    GLOBAL_MATRICES%MATRICES(1)%ELEMENT_MATRIX%MATRIX(mhs,nhs)= &
                      & GLOBAL_MATRICES%MATRICES(1)%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                  ENDDO !ns
                ENDDO !nh
              ENDIF
              IF(GLOBAL_MATRICES%UPDATE_VECTOR) THEN
                GLOBAL_MATRICES%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDIF
            ENDDO !ms
          ENDDO !mh
        ENDDO !ng
        
      CASE(PROBLEM_GENERALISED_LAPLACE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a laplace equation type of a classical field problem class"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_LAPLACE_EQUATION_FINITE_ELEMENT_CALCUALTE")
    RETURN
999 CALL ERRORS("PROBLEM_LAPLACE_EQUATION_FINITE_ELEMENT_CALCUALTE",ERR,ERROR)
    CALL EXITS("PROBLEM_LAPLACE_EQUATION_FINITE_ELEMENT_CALCUALTE")
    RETURN 1
  END SUBROUTINE PROBLEM_LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_LAPLACE_EQUATION_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_LAPLACE_EQUATION_SETUP
    !###  Description:
    !###    Sets up the Laplace equation type of a classical field problem class.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_LAPLACE_EQUATION_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)
        CALL PROBLEM_STANDARD_LAPLACE_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_LAPLACE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a laplace equation type of a classical field problem class"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_LAPLACE_EQUATION_SETUP")
    RETURN
999 CALL ERRORS("PROBLEM_LAPLACE_EQUATION_SETUP",ERR,ERROR)
    CALL EXITS("PROBLEM_LAPLACE_EQUATION_SETUP")
    RETURN 1
  END SUBROUTINE PROBLEM_LAPLACE_EQUATION_SETUP

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET
    !###  Description:
    !###    Sets/changes the problem subtype for a Laplace equation type of a classical field problem class.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)        
        CALL PROBLEM_STANDARD_LAPLACE_SETUP(PROBLEM,PROBLEM_INITIAL_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_LAPLACE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a laplace equation type of a classical field problem class"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_LAPLACE_EQUATION_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_STANDARD_LAPLACE_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_STANDARD_LAPLACE_SETUP
    !###  Description:
    !###    Sets up the standard laplace equation.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEXT_NUMBER
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_STANDARD_LAPLACE_EQUATION_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%CLASS==PROBLEM_CLASSICAL_FIELD_CLASS) THEN
        IF(PROBLEM%TYPE==PROBLEM_LAPLACE_EQUATION_TYPE) THEN
          IF(PROBLEM%SUBTYPE==PROBLEM_STANDARD_LAPLACE_SUBTYPE) THEN
            SELECT CASE(SETUP_TYPE)
            CASE(PROBLEM_INITIAL_SETUP_TYPE)
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                PROBLEM%LINEARITY=PROBLEM_LINEAR
                PROBLEM%TIME_TYPE=PROBLEM_STATIC
                PROBLEM%SOLUTION_METHOD=PROBLEM_FEM_SOLUTION_METHOD
              CASE(PROBLEM_FINISH_ACTION)
                !!TODO: Check valid setup
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_GEOMETRY_SETUP_TYPE)
              !Do nothing???
            CASE(PROBLEM_DEPENDENT_SETUP_TYPE)
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                CALL FIELD_NEXT_NUMBER_FIND(PROBLEM%REGION,NEXT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_CREATE_START(NEXT_NUMBER,PROBLEM%REGION,PROBLEM%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,PROBLEM%GEOMETRY%GEOMETRIC_FIELD% &
                  & DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,PROBLEM%GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,1,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,FIELD_STANDARD_VARIABLE_TYPE,1, &
                  & PROBLEM%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(FIELD_STANDARD_VARIABLE_TYPE)%COMPONENTS(1)%MESH_COMPONENT_NUMBER, &
                  & ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,FIELD_NORMAL_VARIABLE_TYPE,1, &
                  & PROBLEM%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(FIELD_STANDARD_VARIABLE_TYPE)%COMPONENTS(1)%MESH_COMPONENT_NUMBER, &
                  & ERR,ERROR,*999)
                SELECT CASE(PROBLEM%SOLUTION_METHOD)
                CASE(PROBLEM_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,FIELD_STANDARD_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,FIELD_NORMAL_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,PROBLEM%GEOMETRY%GEOMETRIC_FIELD%SCALINGS% &
                    & SCALING_TYPE,ERR,ERROR,*999)
                CASE(PROBLEM_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE(PROBLEM_FINISH_ACTION)
                CALL FIELD_CREATE_FINISH(PROBLEM%REGION,PROBLEM%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_MATERIALS_SETUP_TYPE)
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                !Do nothing
              CASE(PROBLEM_FINISH_ACTION)
                !Do nothing
                !? Maybe set finished flag????
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SOURCE_SETUP_TYPE)
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                !Do nothing
              CASE(PROBLEM_FINISH_ACTION)
                !Do nothing
                !? Maybe set finished flag????
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_ANALYTIC_SETUP_TYPE)
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                IF(PROBLEM%DEPENDENT%DEPENDENT_FINISHED) THEN
                  !Do nothing
                ELSE
                  CALL FLAG_ERROR("Problem dependent has not been finished",ERR,ERROR,*999)
                ENDIF
              CASE(PROBLEM_FINISH_ACTION)
                !Do nothing
                !? Maybe set finished flag????
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_FIXED_CONDITIONS_SETUP_TYPE)
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                IF(PROBLEM%DEPENDENT%DEPENDENT_FINISHED) THEN
                  !Do nothing
                ELSE
                  CALL FLAG_ERROR("Problem dependent has not been finished",ERR,ERROR,*999)
                ENDIF
              CASE(PROBLEM_FINISH_ACTION)
                !Do nothing
                !? Maybe set finished flag????
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SOLUTION_SETUP_TYPE)
              SOLUTION=>PROBLEM%SOLUTION
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                IF(ASSOCIATED(PROBLEM%FIXED_CONDITIONS)) THEN
                  IF(PROBLEM%FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
                    !Do nothing
                    !?Initialise problem solution???
                  ELSE
                    CALL FLAG_ERROR("Problem fixed conditions has not been finished",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Problem fixed conditions is not associated",ERR,ERROR,*999)
                ENDIF
              CASE(PROBLEM_FINISH_ACTION)
                SELECT CASE(PROBLEM%SOLUTION_METHOD)
                CASE(PROBLEM_FEM_SOLUTION_METHOD)
                  !Create the global problem.
                  CALL GLOBAL_MATRICES_CREATE_START(SOLUTION,GLOBAL_MATRICES,ERR,ERROR,*999)
                  CALL GLOBAL_MATRICES_NUMBER_SET(GLOBAL_MATRICES,1,ERR,ERROR,*999)
                  CALL GLOBAL_MATRICES_VARIABLE_TYPES_SET(GLOBAL_MATRICES,(/FIELD_STANDARD_VARIABLE_TYPE/), &
                    & ERR,ERROR,*999)
                  CALL GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET(GLOBAL_MATRICES,FIELD_NORMAL_VARIABLE_TYPE,ERR,ERROR,*999)
                  IF(SOLUTION%GLOBAL_SPARSITY_TYPE==PROBLEM_SOLUTION_FULL_GLOBAL_MATRICES) THEN
                    CALL GLOBAL_MATRICES_STORAGE_TYPE_SET(GLOBAL_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                      & ERR,ERROR,*999)
                  ELSE IF(SOLUTION%GLOBAL_SPARSITY_TYPE==PROBLEM_SOLUTION_SPARSE_GLOBAL_MATRICES) THEN
                    CALL GLOBAL_MATRICES_STORAGE_TYPE_SET(GLOBAL_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                      & ERR,ERROR,*999)
                    CALL GLOBAL_MATRICES_STRUCTURE_TYPE_SET(GLOBAL_MATRICES,(/PROBLEM_GLOBAL_MATRIX_FEM_STRUCTURE/), &
                      & ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The problem solution global matrices type of "// &
                      & TRIM(NUMBER_TO_VSTRING(SOLUTION%GLOBAL_SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid"
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  CALL GLOBAL_MATRICES_CREATE_FINISH(GLOBAL_MATRICES,ERR,ERROR,*999)
                  !Create the solution mapping                  
                  CALL SOLUTION_MAPPING_CREATE_START(GLOBAL_MATRICES,SOLUTION_MAPPING,ERR,ERROR,*999)
                  CALL SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLUTION_MAPPING,1,ERR,ERROR,*999)
                  CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET(SOLUTION_MAPPING,1,(/FIELD_STANDARD_VARIABLE_TYPE/), &
                    & ERR,ERROR,*999)
                  CALL SOLUTION_MAPPING_CREATE_FINISH(SOLUTION_MAPPING,ERR,ERROR,*999)
                  !Create the solver
                  CALL SOLVER_CREATE_START(SOLUTION%SOLUTION_MAPPING,SOLVER_LINEAR_TYPE,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_LIBRARY_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                  CALL SOLVER_CREATE_FINISH(SOLVER,ERR,ERROR,*999)
                CASE(PROBLEM_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT                  
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SOLVER_SETUP_TYPE)
              SOLUTION=>PROBLEM%SOLUTION
              SELECT CASE(ACTION_TYPE)
              CASE(PROBLEM_START_ACTION)
                GLOBAL_MATRICES=>SOLUTION%GLOBAL_MATRICES
                IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
                  !Start the creation of a linear solver
                ELSE
                  CALL FLAG_ERROR("Problem global matrices is not associated",ERR,ERROR,*999)
                ENDIF
              CASE(PROBLEM_FINISH_ACTION)
                SELECT CASE(PROBLEM%SOLUTION_METHOD)
                CASE(PROBLEM_FEM_SOLUTION_METHOD)
                CASE(PROBLEM_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE(PROBLEM_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT                  
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
               LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a standard Laplace equation"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
             END SELECT            
          ELSE
            LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " does not equal a standard Laplace equation subtype"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The problem type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
            & " does not equal a Laplace equation type"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The problem class of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%CLASS,"*",ERR,ERROR))// &
          & " does not equal a classical field problem class"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_STANDARD_LAPLACE_SETUP")
    RETURN
999 CALL ERRORS("PROBLEM_STANDARD_LAPLACE_SETUP",ERR,ERROR)
    CALL EXITS("PROBLEM_STANDARD_LAPLACE_SETUP")
    RETURN 1
  END SUBROUTINE PROBLEM_STANDARD_LAPLACE_SETUP

  !
  !================================================================================================================================
  !

  !>Initialises a problem.
  SUBROUTINE PROBLEM_INITIALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<The pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The errror code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      PROBLEM%USER_NUMBER=0
      PROBLEM%GLOBAL_NUMBER=0
      NULLIFY(PROBLEM%PROBLEMS)
      NULLIFY(PROBLEM%REGION)
      PROBLEM%CLASS=PROBLEM_NO_CLASS
      PROBLEM%TYPE=PROBLEM_NO_TYPE
      PROBLEM%SUBTYPE=PROBLEM_NO_SUBTYPE
      PROBLEM%LINEARITY=0
      PROBLEM%TIME_TYPE=0
      PROBLEM%SOLUTION_METHOD=0
      CALL PROBLEM_GEOMETRY_INITIALISE(PROBLEM,ERR,ERROR,*999)
      CALL PROBLEM_DEPENDENT_INITIALISE(PROBLEM,ERR,ERROR,*999)
      NULLIFY(PROBLEM%MATERIALS)
      NULLIFY(PROBLEM%SOURCE)
      NULLIFY(PROBLEM%ANALYTIC)
      NULLIFY(PROBLEM%FIXED_CONDITIONS)
      NULLIFY(PROBLEM%SOLUTION)
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_FIXED_CONDITIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_FIXED_CONDITIONS_CREATE_FINISH
    !###  Description:
    !###    Finish the creation of fixed conditions for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_computational_nodes,MPI_IERROR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_DEPENDENT_TYPE), POINTER :: PROBLEM_DEPENDENT
    TYPE(PROBLEM_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS

    CALL ENTERS("PROBLEM_FIXED_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      FIXED_CONDITIONS=>PROBLEM%FIXED_CONDITIONS
      IF(ASSOCIATED(PROBLEM%FIXED_CONDITIONS)) THEN
        PROBLEM_DEPENDENT=>PROBLEM%DEPENDENT
        IF(ASSOCIATED(PROBLEM_DEPENDENT)) THEN
          IF(PROBLEM_DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>PROBLEM_DEPENDENT%DEPENDENT_FIELD
            DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
            !Start the transfer of the boundary conditions array. Note that the acutal boundary condition values will be
            !transferred in the assemble routines.
            CALL DISTRIBUTED_VECTOR_UPDATE_START(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
            IF(ERR/=0) GOTO 999
            IF(number_computational_nodes>1) THEN
              !Transfer all the fixed conditions to all the computational nodes. At the moment just use an MPI_ALLREDUCE as the
              !dofs belonging to each computational node are not continuous (different field components) which prevents the
              !straightforward use of MPI_ALLGATHERV. The ALLREDUCE will have more transfers but we will see how bad it is later.
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS,DEPENDENT_DOFS_MAPPING% &
                & NUMBER_OF_GLOBAL,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
            ENDIF
            !Finish the transfer of the boundary conditions
            CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            !Finish problem specific setting up 
            CALL PROBLEM_SETUP(PROBLEM,PROBLEM_FIXED_CONDITIONS_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
            !Finish the fixed conditions
            FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Problem dependent has not been finished",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem dependent is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The problem fixed conditions is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Fixed conditions:",ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DEPENDENT_DOFS_MAPPING%NUMBER_OF_GLOBAL,8,8, &
        & FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS,'("  Global BCs:",8(X,I8))','(13X,8(X,I8))',ERR,ERROR,*999)      
    ENDIF
       
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_FIXED_CONDITIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_FIXED_CONDITIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_FIXED_CONDITIONS_CREATE_START(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_FIXED_CONDITIONS_CREATE_START
    !###  Description:
    !###    Start the creation of fixed conditions for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_DEPENDENT_TYPE), POINTER :: PROBLEM_DEPENDENT

    CALL ENTERS("PROBLEM_FIXED_CONDITIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%FIXED_CONDITIONS)) THEN
        CALL FLAG_ERROR("The problem fixed conditions is already associated",ERR,ERROR,*999)        
      ELSE
        PROBLEM_DEPENDENT=>PROBLEM%DEPENDENT
        IF(ASSOCIATED(PROBLEM_DEPENDENT)) THEN
          IF(PROBLEM_DEPENDENT%DEPENDENT_FINISHED) THEN
            CALL PROBLEM_FIXED_CONDITIONS_INITIALISE(PROBLEM,ERR,ERROR,*999)
            DEPENDENT_FIELD=>PROBLEM_DEPENDENT%DEPENDENT_FIELD
            DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
            ALLOCATE(PROBLEM%FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(DEPENDENT_DOFS_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global fixed conditions",ERR,ERROR,*999)
            PROBLEM%FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS=PROBLEM_NOT_FIXED
            CALL DISTRIBUTED_VECTOR_CREATE_START(DEPENDENT_DOFS_MAPPING,PROBLEM%FIXED_CONDITIONS%BOUNDARY_CONDITIONS, &
              & ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(PROBLEM%FIXED_CONDITIONS%BOUNDARY_CONDITIONS,MATRIX_VECTOR_INTG_TYPE, &
              & ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_CREATE_FINISH(PROBLEM%FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            !Initialise boundary conditions
            CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(PROBLEM%FIXED_CONDITIONS%BOUNDARY_CONDITIONS,PROBLEM_NOT_FIXED,ERR,ERROR,*999)
            !Perform problem specific setup
            CALL PROBLEM_SETUP(PROBLEM,PROBLEM_FIXED_CONDITIONS_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)            
          ELSE
            CALL FLAG_ERROR("The problem dependent has not been finished",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem dependent is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_FIXED_CONDITIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_FIXED_CONDITIONS_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_FIXED_CONDITIONS_FINALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_FIXED_CONDITIONS_FINALISE
    !###  Description:
    !###    Finalise the fixed conditions for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_FIXED_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%FIXED_CONDITIONS)) THEN
        IF(ALLOCATED(PROBLEM%FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS))  &
          & DEALLOCATE(PROBLEM%FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS)
        CALL DISTRIBUTED_VECTOR_DESTROY(PROBLEM%FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
        DEALLOCATE(PROBLEM%FIXED_CONDITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_FIXED_CONDITIONS_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_FIXED_CONDITIONS_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_FIXED_CONDITIONS_INITIALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_FIXED_CONDITIONS_INITIALISE
    !###  Description:
    !###    Initialises the fixed conditions for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_FIXED_CONDITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%FIXED_CONDITIONS)) THEN
        CALL FLAG_ERROR("Fixed conditions is already associated for this problem",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM%FIXED_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem fixed conditions",ERR,ERROR,*999)
        PROBLEM%FIXED_CONDITIONS%PROBLEM=>PROBLEM
        PROBLEM%FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED=.FALSE.
        NULLIFY(PROBLEM%FIXED_CONDITIONS%BOUNDARY_CONDITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_FIXED_CONDITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_FIXED_CONDITIONS_INITIALISE

  !
  !================================================================================================================================
  !
  
  !#### Generic-Subroutine: PROBLEM_FIXED_CONDITIONS_SET_DOF
  !###  Description:
  !###    Sets the fixed conditions for the specified dofs.
  !###  Child-subroutines: PROBLEM_FIXED_CONDITIONS_SET_DOFS,PROBLEM_FIXED_CONDITIONS_SET_DOF1

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_FIXED_CONDITIONS_SET_DOFS(PROBLEM,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_FIXED_CONDITIONS_SET_DOFS
    !###  Description:
    !###    Sets fixed conditions for the problem on the specified dofs.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: DOF_INDICES(:)
    INTEGER(INTG), INTENT(IN) :: CONDITIONS(:)
    REAL(DP), INTENT(IN) :: VALUES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,local_ny,global_ny
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_DEPENDENT_TYPE), POINTER :: PROBLEM_DEPENDENT
    TYPE(PROBLEM_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_FIXED_CONDITIONS_SET_DOFS",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      FIXED_CONDITIONS=>PROBLEM%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
          CALL FLAG_ERROR("Fixed conditions have been finished for this problem",ERR,ERROR,*999)
        ELSE
          PROBLEM_DEPENDENT=>PROBLEM%DEPENDENT
          IF(ASSOCIATED(PROBLEM_DEPENDENT)) THEN
            IF(PROBLEM_DEPENDENT%DEPENDENT_FINISHED) THEN
              IF(SIZE(DOF_INDICES,1)==SIZE(CONDITIONS,1)) THEN
                IF(SIZE(DOF_INDICES,1)==SIZE(VALUES,1)) THEN
                  DEPENDENT_FIELD=>PROBLEM_DEPENDENT%DEPENDENT_FIELD
                  DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
                  DO i=1,SIZE(DOF_INDICES,1)
                    local_ny=DOF_INDICES(i)
                    IF(local_ny>0.AND.local_ny<=DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                      global_ny=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      IF(DEPENDENT_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                        SELECT CASE(CONDITIONS(i))
                        CASE(PROBLEM_NOT_FIXED)
                          FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=PROBLEM_NOT_FIXED                     
                        CASE(PROBLEM_FIXED_BOUNDARY_CONDITION)
!!TODO: need to think how initial conditions and increments for non-linear problems are set.
                          FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=PROBLEM_FIXED_BOUNDARY_CONDITION
                          CALL DISTRIBUTED_VECTOR_VALUES_SET(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,local_ny, &
                            & PROBLEM_FIXED_BOUNDARY_CONDITION,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_UPDATE_DOF(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,local_ny, &
                            & VALUES(i),ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The condition for index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))// &
                          & " is "//TRIM(NUMBER_TO_VSTRING(CONDITIONS(i),"*",ERR,ERROR))//" which is invalid"
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        !?Error
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Invalid dof indices. The dof number for index number "// &
                        & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                        & ". The allowed range is 1 to "// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_DOFS_MAPPING%NUMBER_OF_LOCAL,"*",ERR,ERROR))
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !i
                ELSE
                  LOCAL_ERROR="The size of the dof indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                    & ") does not match the size of the values array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The size of the dof indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                  & ") does not match the size of the fixed conditions array ("// &
                  & TRIM(NUMBER_TO_VSTRING(SIZE(CONDITIONS,1),"*",ERR,ERROR))//")"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem dependent has not been finished",ERR,ERROR,*999)              
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem dependent is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem fixed conditions are not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_SET_DOFS")
    RETURN
999 CALL ERRORS("PROBLEM_FIXED_CONDITIONS_SET_DOFS",ERR,ERROR)
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_SET_DOFS")
    RETURN 1
  END SUBROUTINE PROBLEM_FIXED_CONDITIONS_SET_DOFS

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_FIXED_CONDITIONS_SET_DOF1(PROBLEM,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_FIXED_CONDITIONS_SET_DOF1
    !###  Description:
    !###    Sets a fixed condition for the problem on the specified dof.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX
    INTEGER(INTG), INTENT(IN) :: CONDITION
    REAL(DP), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(PROBLEM_DEPENDENT_TYPE), POINTER :: DEPENDENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_FIXED_CONDITIONS_SET_DOF1",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      FIXED_CONDITIONS=>PROBLEM%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
          CALL FLAG_ERROR("Fixed conditions have been finished for this problem",ERR,ERROR,*999)
        ELSE
          DEPENDENT=>PROBLEM%DEPENDENT
          IF(ASSOCIATED(DEPENDENT)) THEN
            IF(DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>DEPENDENT%DEPENDENT_FIELD
              DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
              local_ny=DOF_INDEX
              IF(local_ny>0.AND.local_ny<=DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                global_ny=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                SELECT CASE(CONDITION)
                CASE(PROBLEM_NOT_FIXED)
                  FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=PROBLEM_NOT_FIXED                     
                CASE(PROBLEM_FIXED_BOUNDARY_CONDITION)
                  FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=PROBLEM_FIXED_BOUNDARY_CONDITION
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,local_ny, &
                    & PROBLEM_FIXED_BOUNDARY_CONDITION,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_DOF(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,local_ny, &
                    & VALUE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The specified condition of "//TRIM(NUMBER_TO_VSTRING(CONDITION,"*",ERR,ERROR))// &
                    & " is invalid"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                LOCAL_ERROR="The specified dof of "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                  & "is invalid. The allowed range is 1 to "// &
                  & TRIM(NUMBER_TO_VSTRING(DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem dependent has not been finished",ERR,ERROR,*999)              
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem dependent is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem fixed conditions are not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_SET_DOF1")
    RETURN
999 CALL ERRORS("PROBLEM_FIXED_CONDITIONS_SET_DOF1",ERR,ERROR)
    CALL EXITS("PROBLEM_FIXED_CONDITIONS_SET_DOF1")
    RETURN 1
  END SUBROUTINE PROBLEM_FIXED_CONDITIONS_SET_DOF1

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_GEOMETRY_FINALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_GEOMETRY_FINALISE
    !###  Description:
    !###    Finalise the geometry for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_GEOMETRY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Do nothing
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_GEOMETRY_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_GEOMETRY_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_GEOMETRY_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_GEOMETRY_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_GEOMETRY_INITIALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_GEOMETRY_INITIALISE
    !###  Description:
    !###    Initialises the geometry for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_GEOMETRY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      PROBLEM%GEOMETRY%PROBLEM=>PROBLEM
      NULLIFY(PROBLEM%GEOMETRY%GEOMETRIC_FIELD)
      NULLIFY(PROBLEM%GEOMETRY%FIBRE_FIELD)
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_GEOMETRY_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_GEOMETRY_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_GEOMETRY_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_GEOMETRY_INITIALISE

  !
  !================================================================================================================================
  !
  
  !> Finalises the linear data information for a problem and deallocates all memory
  SUBROUTINE PROBLEM_LINEAR_DATA_FINALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_LINEAR_DATA_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%LINEAR_DATA)) THEN
        DEALLOCATE(PROBLEM_SOLUTION%LINEAR_DATA)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_LINEAR_DATA_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_LINEAR_DATA_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_LINEAR_DATA_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_LINEAR_DATA_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the linear data information for a problem solution
  SUBROUTINE PROBLEM_LINEAR_DATA_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<The pointer to the problem solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_LINEAR_DATA_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%LINEAR_DATA)) THEN
        CALL FLAG_ERROR("Linear data is already associated for this problem solution",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLUTION%LINEAR_DATA,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution linear data",ERR,ERROR,*999)
        PROBLEM_SOLUTION%LINEAR_DATA%PROBLEM_SOLUTION=>PROBLEM_SOLUTION        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_LINEAR_DATA_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_LINEAR_DATA_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_LINEAR_DATA_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_LINEAR_DATA_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET(PROBLEM,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET
    !###  Description:
    !###    Sets/changes the field component interpolation for a materials field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        IF(PROBLEM%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Problem materials has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_COMPONENT_INTERPOLATION_SET(PROBLEM%MATERIALS%MATERIAL_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
            & COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem materials is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_COMPONENT_INTERPOLATION_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET(PROBLEM,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET
    !###  Description:
    !###    Sets/changes the field component mesh component for a materials field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        IF(PROBLEM%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Problem materials has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(PROBLEM%MATERIALS%MATERIAL_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
            & COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem materials is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_COMPONENT_MESH_COMPONENT_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_CREATE_FINISH
    !###  Description:
    !###    Finish the creation of materials for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_MATERIALS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        !Finish problem specific startup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_MATERIALS_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
        !Finish materials creation
        PROBLEM%MATERIALS%MATERIALS_FINISHED=.FALSE.
      ELSE
        CALL FLAG_ERROR("The problem materials is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_CREATE_START(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_CREATE_START
    !###  Description:
    !###    Start the creation of materials for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_MATERIALS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        CALL FLAG_ERROR("The problem materials is already associated",ERR,ERROR,*999)        
      ELSE
        CALL PROBLEM_MATERIALS_INITIALISE(PROBLEM,ERR,ERROR,*999)
        !Start problem specific startup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_MATERIALS_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_FINALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_FINALISE
    !###  Description:
    !###    Finalise the materials for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_MATERIALS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        IF(ASSOCIATED(PROBLEM%MATERIALS%MATERIAL_FIELD))  &
          & CALL FIELD_DESTROY(PROBLEM%MATERIALS%MATERIAL_FIELD,PROBLEM%REGION,ERR,ERROR,*999)
        DEALLOCATE(PROBLEM%MATERIALS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_INITIALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_INITIALISE
    !###  Description:
    !###    Initialises the materials for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_MATERIALS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        CALL FLAG_ERROR("Materials is already associated for this problem",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM%MATERIALS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem materials",ERR,ERROR,*999)
        PROBLEM%MATERIALS%PROBLEM=>PROBLEM
        PROBLEM%MATERIALS%MATERIALS_FINISHED=.FALSE.
        NULLIFY(PROBLEM%MATERIALS%MATERIAL_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_MATERIAL_FIELD_GET(PROBLEM,MATERIAL_FIELD,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_MATERIAL_FIELD_GET
    !###  Description:
    !###    Returns a pointer to the material field of the materials for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(FIELD_TYPE), POINTER :: MATERIAL_FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_MATERIALS_MATERIAL_FIELD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        IF(PROBLEM%MATERIALS%MATERIALS_FINISHED) THEN
          MATERIAL_FIELD=>PROBLEM%MATERIALS%MATERIAL_FIELD
        ELSE
          CALL FLAG_ERROR("Materials has not been finished for this problem",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Materials is not associated for this problem",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_MATERIAL_FIELD_GET")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_MATERIAL_FIELD_GET",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_MATERIAL_FIELD_GET")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_MATERIAL_FIELD_GET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_MATERIALS_SCALING_SET(PROBLEM,SCALING_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_MATERIALS_SCALING_SET
    !###  Description:
    !###    Sets/changes the field scaling for a materials field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_MATERIALS_SCALING_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%MATERIALS)) THEN
        IF(PROBLEM%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Problem materials has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_SCALING_TYPE_SET(PROBLEM%MATERIALS%MATERIAL_FIELD,SCALING_TYPE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem materials is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_MATERIALS_SCALING_SET")
    RETURN
999 CALL ERRORS("PROBLEM_MATERIALS_SCALING_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_MATERIALS_SCALING_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_MATERIALS_SCALING_SET

  !
  !================================================================================================================================
  !
  
  !> Finalises the nonlinear data information for a problem and deallocates all memory
  SUBROUTINE PROBLEM_NONLINEAR_DATA_FINALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_NONLINEAR_DATA_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%NONLINEAR_DATA)) THEN
        DEALLOCATE(PROBLEM_SOLUTION%NONLINEAR_DATA)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_NONLINEAR_DATA_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_NONLINEAR_DATA_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_NONLINEAR_DATA_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_NONLINEAR_DATA_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nonlinear data information for a problem solution
  SUBROUTINE PROBLEM_NONLINEAR_DATA_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<The pointer to the problem solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_NONLINEAR_DATA_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%NONLINEAR_DATA)) THEN
        CALL FLAG_ERROR("Nonlinear data is already associated for this problem solution",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLUTION%NONLINEAR_DATA,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution nonlinear data",ERR,ERROR,*999)
        PROBLEM_SOLUTION%NONLINEAR_DATA%PROBLEM_SOLUTION=>PROBLEM_SOLUTION        
        PROBLEM_SOLUTION%NONLINEAR_DATA%NUMBER_OF_ITERATIONS=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_NONLINEAR_DATA_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_NONLINEAR_DATA_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_NONLINEAR_DATA_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_NONLINEAR_DATA_INITIALISE

   !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET(PROBLEM,VARIABLE_NUMBER,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER, &
    & ERR,ERROR,*)

    !#### Subroutine: PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET
    !###  Description:
    !###    Sets/changes the field component mesh component for a dependent field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("Problem dependent has been finished",ERR,ERROR,*999)
      ELSE
        CALL FIELD_COMPONENT_MESH_COMPONENT_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER, &
          & MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET")
    RETURN
999 CALL ERRORS("PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_DEPENDENT_COMPONENT_MESH_COMPONENT_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_DEPENDENT_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_DEPENDENT_CREATE_FINISH
    !###  Description:
    !###    Finish the creation of a dependent variables for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_DEPENDENT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Finish problem specific setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_DEPENDENT_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
      !Finish the problem creation
      PROBLEM%DEPENDENT%DEPENDENT_FINISHED=.TRUE.
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_DEPENDENT_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_DEPENDENT_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_DEPENDENT_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_DEPENDENT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_DEPENDENT_CREATE_START(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_DEPENDENT_CREATE_START
    !###  Description:
    !###    Start the creation of dependent variables for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_DEPENDENT_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Start the problem specfic solution setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_DEPENDENT_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_DEPENDENT_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_DEPENDENT_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_DEPENDENT_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_DEPENDENT_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET(PROBLEM,DEPENDENT_FIELD,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET
    !###  Description:
    !###    Returns a pointer to the dependent field of the dependent variables for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%DEPENDENT%DEPENDENT_FINISHED) THEN
        DEPENDENT_FIELD=>PROBLEM%DEPENDENT%DEPENDENT_FIELD
      ELSE
        CALL FLAG_ERROR("Dependent variables have not been finished for this problem",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET")
    RETURN
999 CALL ERRORS("PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET",ERR,ERROR)
    CALL EXITS("PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET")
    RETURN 1
  END SUBROUTINE PROBLEM_DEPENDENT_DEPENDENT_FIELD_GET

  !
  !================================================================================================================================
  !
  
  !>Finalises the dependent variables for a problem and deallocates all memory.
  SUBROUTINE PROBLEM_DEPENDENT_FINALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<The pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_DEPENDENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%DEPENDENT%DEPENDENT_FIELD))  &
        & CALL FIELD_DESTROY(PROBLEM%DEPENDENT%DEPENDENT_FIELD,PROBLEM%REGION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_DEPENDENT_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_DEPENDENT_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_DEPENDENT_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_DEPENDENT_FINALISE
  
  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_DEPENDENT_INITIALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_DEPENDENT_INITIALISE
    !###  Description:
    !###    Initialises the dependent variables for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_DEPENDENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      PROBLEM%DEPENDENT%PROBLEM=>PROBLEM
      PROBLEM%DEPENDENT%DEPENDENT_FINISHED=.FALSE.
      NULLIFY(PROBLEM%DEPENDENT%DEPENDENT_FIELD)
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_DEPENDENT_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_DEPENDENT_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_DEPENDENT_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_DEPENDENT_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_DEPENDENT_SCALING_SET(PROBLEM,SCALING_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_DEPENDENT_SCALING_SET
    !###  Description:
    !###    Sets/changes the field scaling for a dependent field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_DEPENDENT_SCALING_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("Problem dependent has been finished",ERR,ERROR,*999)
      ELSE
        CALL FIELD_SCALING_TYPE_SET(PROBLEM%DEPENDENT%DEPENDENT_FIELD,SCALING_TYPE,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_DEPENDENT_SCALING_SET")
    RETURN
999 CALL ERRORS("PROBLEM_DEPENDENT_SCALING_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_DEPENDENT_SCALING_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_DEPENDENT_SCALING_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SETUP
    !###  Description:
    !###    Sets up the specifices for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%CLASS)
      CASE(PROBLEM_ELASTICITY_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_FLUID_MECHANICS_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
        CALL PROBLEM_CLASSICAL_FIELD_CLASS_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(PROBLEM_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM%CLASS,"*",ERR,ERROR))//" is not valid"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOLVE
    !###  Description:
    !###    Solves a problem.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOLUTION)) THEN
        IF(PROBLEM%SOLUTION%SOLUTION_FINISHED) THEN
          SELECT CASE(PROBLEM%LINEARITY)
          CASE(PROBLEM_LINEAR)
            SELECT CASE(PROBLEM%TIME_TYPE)
            CASE(PROBLEM_STATIC)
              SELECT CASE(PROBLEM%SOLUTION_METHOD)
              CASE(PROBLEM_FEM_SOLUTION_METHOD)
                CALL PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM(PROBLEM,ERR,ERROR,*999)
                CALL PROBLEM_SOLVE_LINEAR_STATIC_FEM(PROBLEM,ERR,ERROR,*999)
              CASE(PROBLEM_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
              CASE(PROBLEM_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
              CASE(PROBLEM_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
              CASE(PROBLEM_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
              CASE(PROBLEM_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The problem solution method  of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_DYNAMIC,PROBLEM_QUASISTATIC)
              CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The problem time type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TIME_TYPE,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_NONLINEAR,PROBLEM_NONLINEAR_BCS)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The problem linearity of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%LINEARITY,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem solution has not been finished",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOLVE_LINEAR_STATIC_FEM(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOLVE_LINEAR_STATIC_FEM
    !###  Description:
    !###    Solves a linear static problem using the finite element method.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVE_LINEAR_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Set up the FEM problem and allocate solution arrays
      CALL PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM(PROBLEM,ERR,ERROR,*999)
      !Assemble the finite element equations into the global stiffness matrix and rhs vector
      CALL PROBLEM_ASSEMBLE_LINEAR_STATIC_FEM(PROBLEM,ERR,ERROR,*999)
      !Calculate the solution matrix and right hand side
      !Solve the system
      !Back substitute to find flux values.
      !Update flux values back to the field vector.
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVE_LINEAR_STATIC_FEM")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVE_LINEAR_STATIC_FEM",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVE_LINEAR_STATIC_FEM")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVE_LINEAR_STATIC_FEM

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM
    !###  Description:
    !###    Sets up a linear static problem using the finite element method.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Do nothing???
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVE_SETUP_LINEAR_STATIC_FEM

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a global matrix.
  SUBROUTINE PROBLEM_GLOBAL_MATRIX_STRUCTURE_CALCULATE(GLOBAL_MATRICES,matrix_idx,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(IN) :: matrix_idx !<The matrix number to calculate the structure of
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return a pointer to row location indices. The calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return a pointer to the column location indices. The calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) ::  column_idx,DUMMY_ERR,elem_idx,global_column,local_column,local_ny,mk,mp,ne,nh,nn,nnk,np, &
      & NUMBER_OF_COLUMNS,nyy
    INTEGER(INTG), POINTER :: COLUMNS(:)
    REAL(DP) :: SPARSITY
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_DOMAIN_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES    
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: DEPENDENT_DOFS_PARAM_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(COLUMNS)
    
    CALL ENTERS("PROBLEM_GLOBAL_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*998)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(matrix_idx>=1.AND.matrix_idx<=GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN            
        IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
          IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
            SELECT CASE(GLOBAL_MATRICES%MATRICES(matrix_idx)%STRUCTURE_TYPE)
            CASE(PROBLEM_GLOBAL_MATRIX_NO_STRUCTURE)
              CALL FLAG_ERROR("Not implemented",ERR,ERROR,*998)
            CASE(PROBLEM_GLOBAL_MATRIX_FEM_STRUCTURE)
              SELECT CASE(GLOBAL_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE)
              CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                DEPENDENT_FIELD=>GLOBAL_MATRICES%PROBLEM_SOLUTION%PROBLEM%DEPENDENT%DEPENDENT_FIELD
                FIELD_VARIABLE=>GLOBAL_MATRICES%MATRICES(matrix_idx)%VARIABLE
                DEPENDENT_DOFS_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                DEPENDENT_DOFS_PARAM_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOF_TO_PARAM_MAP
                !Allocate lists
                ALLOCATE(COLUMN_INDICES_LISTS(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists",ERR,ERROR,*999)
                !Allocate row indices
                ALLOCATE(ROW_INDICES(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices",ERR,ERROR,*999)
                ROW_INDICES(1)=1
                !First, loop over the rows and calculate the number of non-zeros
                NUMBER_OF_NON_ZEROS=0
                DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                  IF(DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(1,local_ny)==FIELD_NODE_DOF_TYPE) THEN
                    nyy=DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(2,local_ny)
                    np=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(2,nyy)
                    nh=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(3,nyy)
                    DOMAIN_NODES=>FIELD_VARIABLE%COMPONENTS(nh)%DOMAIN%TOPOLOGY%NODES
                    DOMAIN_ELEMENTS=>FIELD_VARIABLE%COMPONENTS(nh)%DOMAIN%TOPOLOGY%ELEMENTS
                    !Set up list
                    NULLIFY(COLUMN_INDICES_LISTS(local_ny)%PTR)
                    CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR,DOMAIN_NODES%NODES(np)% &
                      & NUMBER_OF_SURROUNDING_ELEMENTS*FIELD_VARIABLE%COMPONENTS(nh)%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS, &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                    !Loop over all elements containing the dof
                    DO elem_idx=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                      ne=DOMAIN_NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
                      BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                      DO nn=1,BASIS%NUMBER_OF_NODES
                        mp=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                        DO nnk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                          mk=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                          !Find the local and global column and add the global column to the indices list
                          local_column=FIELD_VARIABLE%COMPONENTS(nh)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(mk,mp,1)
                          global_column=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                          CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_ny)%PTR,global_column,ERR,ERROR,*999)
                        ENDDO !mk
                      ENDDO !nn
                    ENDDO !elem_idx
                    CALL LIST_REMOVE_DUPLICATES(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                    CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                    NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                    ROW_INDICES(local_ny+1)=NUMBER_OF_NON_ZEROS+1
                  ELSE
                    LOCAL_ERROR="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))//" is not a node based dof"
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !local_ny
                !Allocate and setup the column locations
                ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices",ERR,ERROR,*999)
                DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                  CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS,COLUMNS,ERR,ERROR,*999)
                  DO column_idx=1,NUMBER_OF_COLUMNS
                    COLUMN_INDICES(ROW_INDICES(local_ny)+column_idx-1)=COLUMNS(column_idx)
                  ENDDO !column_idx
                  DEALLOCATE(COLUMNS)
                ENDDO !local_ny
                IF(DIAGNOSTICS1) THEN
                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Global matrix structure:",ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Global matrix number : ",matrix_idx,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",DEPENDENT_DOFS_DOMAIN_MAPPING% &
                    & TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",DEPENDENT_DOFS_DOMAIN_MAPPING% &
                    & NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                  IF(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL*DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL/=0) THEN
                    SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL* &
                      DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,DP)*100.0_DP
                    CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY,"F5.2",ERR,ERROR,*999)
                  ENDIF
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DEPENDENT_DOFS_DOMAIN_MAPPING% &
                    & TOTAL_NUMBER_OF_LOCAL+1,8,8,ROW_INDICES,'("  Row indices    :",8(X,I13))','(18X,8(X,I13))', &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8,COLUMN_INDICES, &
                    & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', ERR,ERROR,*999)
                ENDIF
              CASE DEFAULT
                LOCAL_ERROR="The matrix storage type of "// &
                  & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE,"*",ERR,ERROR))// &
                  & " is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The matrix structure type of "// &
                & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%MATRICES(matrix_idx)%STRUCTURE_TYPE,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Column indices is already associated",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Row indieces is already associated",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The matrix index of "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
          & " is invalid. The index must be >= 1 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("PROBLEM_GLOBAL_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(COLUMNS)) DEALLOCATE(COLUMNS)
    IF(ALLOCATED(COLUMN_INDICES_LISTS)) THEN
      DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(COLUMN_INDICES_LISTS(local_ny)%PTR)) &
          & CALL LIST_DESTROY(COLUMN_INDICES_LISTS(local_ny)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_ny
      DEALLOCATE(COLUMN_INDICES_LISTS)
    ENDIF
998 CALL ERRORS("PROBLEM_GLOBAL_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)
    CALL EXITS("PROBLEM_GLOBAL_MATRIX_STRUCTURE_CALCULATE")
    RETURN 1
  END SUBROUTINE PROBLEM_GLOBAL_MATRIX_STRUCTURE_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finish the creation of a solution for a problem.
  SUBROUTINE PROBLEM_SOLUTION_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: SOLUTION
    
    CALL ENTERS("PROBLEM_SOLUTION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SOLUTION=>PROBLEM%SOLUTION
      IF(ASSOCIATED(SOLUTION)) THEN
        IF(SOLUTION%SOLUTION_FINISHED) THEN
          CALL FLAG_ERROR("Problem solution has been finished",ERR,ERROR,*999)
        ELSE
          !Finish the problem specific solution setup.
          CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SOLUTION_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
          !Finish the problem solution creation
          SOLUTION%SOLUTION_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The problem solution is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a solution for the problem. \todo Should this return a pointer to the problem solution???
  SUBROUTINE PROBLEM_SOLUTION_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to create a solution for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLUTION_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOLUTION)) THEN
        CALL FLAG_ERROR("The problem solution is already associated",ERR,ERROR,*999)        
      ELSE
        !Initialise the problem solver
        CALL PROBLEM_SOLUTION_INITIALISE(PROBLEM,ERR,ERROR,*999)
        !Start the problem specific solution setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SOLUTION_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOLUTION_FINALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOLUTION_FINALISE
    !###  Description:
    !###    Finalise the solution for a problem and deallocate all memory.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLUTION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOLUTION)) THEN
        IF(ASSOCIATED(PROBLEM%SOLUTION%INTERPOLATION)) CALL PROBLEM_INTERPOLATION_FINALISE(PROBLEM%SOLUTION,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM%SOLUTION%LINEAR_DATA)) CALL PROBLEM_LINEAR_DATA_FINALISE(PROBLEM%SOLUTION,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM%SOLUTION%NONLINEAR_DATA)) CALL PROBLEM_NONLINEAR_DATA_FINALISE(PROBLEM%SOLUTION,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM%SOLUTION%TIME_DATA)) CALL PROBLEM_TIME_DATA_FINALISE(PROBLEM%SOLUTION,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM%SOLUTION%GLOBAL_MATRICES)) CALL GLOBAL_MATRICES_DESTROY(PROBLEM%SOLUTION%GLOBAL_MATRICES, &
          & ERR,ERROR,*999)
        !IF(ASSOCIATED(PROBLEM%SOLUTION%SOLVER_MATRICES)) CALL PROBLEM_SOLVER_MATRICES_FINALISE(PROBLEM%SOLUTION,ERR,ERROR,*999)
        DEALLOCATE(PROBLEM%SOLUTION)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_FINALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for the problem solution
  SUBROUTINE PROBLEM_SOLUTION_GLOBAL_SPARSITY_TYPE_SET(PROBLEM,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: SPARSITY_TYPE !<The sparsity type to set \see PROBLEM_ROUTINES_SolutionGlobalSparsityTypes,PROBLEM_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_SOLUTION_GLOBAL_SPARSITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOLUTION)) THEN
        IF(PROBLEM%SOLUTION%SOLUTION_FINISHED) THEN
          CALL FLAG_ERROR("Problem solution has already been finished",ERR,ERROR,*999)
        ELSE
          SELECT CASE(SPARSITY_TYPE)
          CASE(PROBLEM_SOLUTION_SPARSE_GLOBAL_MATRICES)
            PROBLEM%SOLUTION%GLOBAL_SPARSITY_TYPE=PROBLEM_SOLUTION_SPARSE_GLOBAL_MATRICES
          CASE(PROBLEM_SOLUTION_FULL_GLOBAL_MATRICES)
            PROBLEM%SOLUTION%GLOBAL_SPARSITY_TYPE=PROBLEM_SOLUTION_FULL_GLOBAL_MATRICES
          CASE DEFAULT
            LOCAL_ERROR="The specified global sparsity type of "//TRIM(NUMBER_TO_VSTRING(SPARSITY_TYPE,"*",ERR,ERROR))// &
              & " is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_GLOBAL_SPARSITY_TYPE_SET")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_GLOBAL_SPARSITY_TYPE_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_GLOBAL_SPARSITY_TYPE_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_GLOBAL_SPARSITY_TYPE_SET
  
  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOLUTION_INITIALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOLUTION_INITIALISE
    !###  Description:
    !###    Initialises the solution for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOLUTION_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOLUTION)) THEN
        CALL FLAG_ERROR("Solution is already associated for this problem",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM%SOLUTION,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution",ERR,ERROR,*999)
        PROBLEM%SOLUTION%PROBLEM=>PROBLEM
        PROBLEM%SOLUTION%OUTPUT_TYPE=PROBLEM_SOLUTION_NO_OUTPUT
        PROBLEM%SOLUTION%GLOBAL_SPARSITY_TYPE=PROBLEM_SOLUTION_SPARSE_GLOBAL_MATRICES
        NULLIFY(PROBLEM%SOLUTION%INTERPOLATION)
        NULLIFY(PROBLEM%SOLUTION%LINEAR_DATA)
        NULLIFY(PROBLEM%SOLUTION%NONLINEAR_DATA)
        NULLIFY(PROBLEM%SOLUTION%TIME_DATA)
        NULLIFY(PROBLEM%SOLUTION%SOLVER)
        NULLIFY(PROBLEM%SOLUTION%GLOBAL_MATRICES)
        !NULLIFY(PROBLEM%SOLUTION%SOLVER_MATRICES)
        PROBLEM%SOLUTION%SOLUTION_FINISHED=.FALSE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for the problem solution
  SUBROUTINE PROBLEM_SOLUTION_OUTPUT_TYPE_SET(PROBLEM,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_SOLUTION_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOLUTION)) THEN
        IF(PROBLEM%SOLUTION%SOLUTION_FINISHED) THEN
          CALL FLAG_ERROR("Problem solution has already been finished",ERR,ERROR,*999)
        ELSE
          SELECT CASE(OUTPUT_TYPE)
          CASE(PROBLEM_SOLUTION_NO_OUTPUT)
            PROBLEM%SOLUTION%OUTPUT_TYPE=PROBLEM_SOLUTION_NO_OUTPUT
          CASE(PROBLEM_SOLUTION_TIMING_OUTPUT)
            PROBLEM%SOLUTION%OUTPUT_TYPE=PROBLEM_SOLUTION_TIMING_OUTPUT
          CASE(PROBLEM_SOLUTION_GLOBAL_MATRIX_OUTPUT)
            PROBLEM%SOLUTION%OUTPUT_TYPE=PROBLEM_SOLUTION_GLOBAL_MATRIX_OUTPUT
          CASE(PROBLEM_SOLUTION_ELEMENT_MATRIX_OUTPUT)
            PROBLEM%SOLUTION%OUTPUT_TYPE=PROBLEM_SOLUTION_ELEMENT_MATRIX_OUTPUT
          CASE DEFAULT
            LOCAL_ERROR="The specified output type of "//TRIM(NUMBER_TO_VSTRING(OUTPUT_TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_OUTPUT_TYPE_SET")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_OUTPUT_TYPE_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_OUTPUT_TYPE_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_OUTPUT_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the solver matrices for the problem solution and deallocates all memory.
  SUBROUTINE PROBLEM_SOLVER_MATRICES_FINALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVER_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      !IF(ASSOCIATED(PROBLEM_SOLUTION%SOLVER_MATRICES)) THEN
      !  IF(ALLOCATED(PROBLEM_SOLUTION%SOLVER_MATRICES%SOLVER_TO_GLOBAL_MAP)) THEN
      !    DO no=1,PROBLEM_SOLUTION%SOLVER_MATRICES%NUMBER_OF_COLUMNS
      !      CALL PROBLEM_SOLVER_TO_GLOBAL_MAP_FINALISE(PROBLEM_SOLUTION%SOLVER_MATRICES%SOLVER_TO_GLOBAL_MAP(no),ERR,ERROR,*999)
      !    ENDDO !no
      !    DEALLOCATE(PROBLEM_SOLUTION%SOLVER_MATRICES%SOLVER_TO_GLOBAL_MAP)
      !  ENDIF
      !  DEALLOCATE(PROBLEM_SOLUTION%SOLVER_MATRICES)
      !ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_MATRICES_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_MATRICES_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_MATRICES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver matrices for the problem solution.
  SUBROUTINE PROBLEM_SOLVER_MATRICES_INITALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVER_MATRICES_INITALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      !IF(ASSOCIATED(PROBLEM_SOLUTION%SOLVER_MATRICES)) THEN
      !  CALL FLAG_ERROR("Solver matrices is already associated for this problem solution",ERR,ERROR,*999)
      !ELSE
      !  ALLOCATE(PROBLEM_SOLUTION%SOLVER_MATRICES,STAT=ERR)
      !  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution solver matrices",ERR,ERROR,*999)
      !  PROBLEM_SOLUTION%SOLVER_MATRICES%PROBLEM_SOLUTION=>PROBLEM_SOLUTION
      !  PROBLEM_SOLUTION%SOLVER_MATRICES%NUMBER_OF_ROWS=0
      !  PROBLEM_SOLUTION%SOLVER_MATRICES%TOTAL_NUMBER_OF_ROWS=0
      !  PROBLEM_SOLUTION%SOLVER_MATRICES%NUMBER_OF_COLUMNS=0
      !  PROBLEM_SOLUTION%SOLVER_MATRICES%UPDATE_MATRIX=.FALSE.
      !ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_MATRICES_INITALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_MATRICES_INITALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_MATRICES_INITALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_MATRICES_INITALISE

  !
  !================================================================================================================================
  !

  !>Finalise a solution to global matrix variable map and deallocate all memory
  SUBROUTINE PROBLEM_SOLVER_TO_GLOBAL_MAP_FINALISE(SOLVER_TO_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TO_GLOBAL_MAP_TYPE):: SOLVER_TO_GLOBAL_MAP !<The solver to global map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    
    CALL ENTERS("PROBLEM_SOLVER_TO_GLOBAL_MAP_FINALISE",ERR,ERROR,*999)

    !IF(ALLOCATED(SOLVER_TO_GLOBAL_MAP%COLUMNS)) DEALLOCATE(SOLVER_TO_GLOBAL_MAP%COLUMNS)
    !IF(ALLOCATED(SOLVER_TO_GLOBAL_MAP%COUPLING_COEFFICIENTS)) DEALLOCATE(SOLVER_TO_GLOBAL_MAP%COUPLING_COEFFICIENTS)
    
    CALL EXITS("PROBLEM_SOLVER_TO_GLOBAL_MAP_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_TO_GLOBAL_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_TO_GLOBAL_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_TO_GLOBAL_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the solver to global matrix variable map.
  SUBROUTINE PROBLEM_SOLVER_TO_GLOBAL_MAP_INITIALISE(SOLVER_TO_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TO_GLOBAL_MAP_TYPE) :: SOLVER_TO_GLOBAL_MAP !<The solver to global map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVER_TO_GLOBAL_MAP_INITIALISE",ERR,ERROR,*999)

    !SOLVER_TO_GLOBAL_MAP%NUMBER_OF_MATRICES=0
       
    CALL EXITS("PROBLEM_SOLVER_TO_GLOBAL_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_TO_GLOBAL_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_TO_GLOBAL_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_TO_GLOBAL_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET(PROBLEM,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET
    !###  Description:
    !###    Sets/changes the field component interpolation for a source field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOURCE)) THEN
        IF(PROBLEM%SOURCE%SOURCE_FINISHED) THEN
          CALL FLAG_ERROR("Problem source has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_COMPONENT_INTERPOLATION_SET(PROBLEM%SOURCE%SOURCE_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
            & COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem source is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET")
    RETURN
999 CALL ERRORS("PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_SOURCE_COMPONENT_INTERPOLATION_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET(PROBLEM,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET
    !###  Description:
    !###    Sets/changes the field component mesh_component for a source field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOURCE)) THEN
        IF(PROBLEM%SOURCE%SOURCE_FINISHED) THEN
          CALL FLAG_ERROR("Problem source has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(PROBLEM%SOURCE%SOURCE_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
            & COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem source is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET")
    RETURN
999 CALL ERRORS("PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_SOURCE_COMPONENT_MESH_COMPONENT_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOURCE_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOURCE_CREATE_FINISH
    !###  Description:
    !###    Finish the creation of a source for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOURCE_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOURCE)) THEN
        !Finish the problem specific source setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SOURCE_SETUP_TYPE,PROBLEM_FINISH_ACTION,ERR,ERROR,*999)
        !Finish the source creation
        PROBLEM%SOURCE%SOURCE_FINISHED=.FALSE.
      ELSE
        CALL FLAG_ERROR("The problem source is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOURCE_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOURCE_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOURCE_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOURCE_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOURCE_CREATE_START(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOURCE_CREATE_START
    !###  Description:
    !###    Start the creation of a source for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOURCE_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOURCE)) THEN
        CALL FLAG_ERROR("The problem source is already associated",ERR,ERROR,*999)        
      ELSE
        !Initialise the problem source
        CALL PROBLEM_SOURCE_INITIALISE(PROBLEM,ERR,ERROR,*999)
        !Start the problem specific source setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SOURCE_SETUP_TYPE,PROBLEM_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOURCE_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOURCE_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOURCE_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOURCE_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOURCE_FINALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOURCE_FINALISE
    !###  Description:
    !###    Finalise the source for a problem and deallocate all memory.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOURCE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOURCE)) THEN
        IF(ASSOCIATED(PROBLEM%SOURCE%SOURCE_FIELD))  &
          & CALL FIELD_DESTROY(PROBLEM%SOURCE%SOURCE_FIELD,PROBLEM%REGION,ERR,ERROR,*999)
        DEALLOCATE(PROBLEM%SOURCE)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOURCE_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOURCE_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOURCE_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOURCE_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOURCE_INITIALISE(PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOURCE_INITIALISE
    !###  Description:
    !###    Initialises the source for a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOURCE_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOURCE)) THEN
        CALL FLAG_ERROR("Source is already associated for this problem",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM%SOURCE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem source",ERR,ERROR,*999)
        PROBLEM%SOURCE%PROBLEM=>PROBLEM
        PROBLEM%SOURCE%SOURCE_FINISHED=.FALSE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOURCE_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOURCE_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOURCE_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOURCE_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SOURCE_SCALING_SET(PROBLEM,SCALING_TYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SOURCE_SCALING_SET
    !###  Description:
    !###    Sets/changes the field scaling for a source field of a problem.

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOURCE_SCALING_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%SOURCE)) THEN
        IF(PROBLEM%SOURCE%SOURCE_FINISHED) THEN
          CALL FLAG_ERROR("Problem source has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_SCALING_TYPE_SET(PROBLEM%SOURCE%SOURCE_FIELD,SCALING_TYPE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem source is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOURCE_SCALING_SET")
    RETURN
999 CALL ERRORS("PROBLEM_SOURCE_SCALING_SET",ERR,ERROR)
    CALL EXITS("PROBLEM_SOURCE_SCALING_SET")
    RETURN 1
  END SUBROUTINE PROBLEM_SOURCE_SCALING_SET

  !
  !================================================================================================================================
  !
  
  !#### Generic-Subroutine: PROBLEM_SPECIFICATION_SET
  !###  Description:
  !###    Sets/changes the problem specification i.e., problem class, type and subtype for a problem.
  !###  Child-subroutines: PROBLEM_SPECIFICATION_SET_NUMBER,PROBLEM_SPECIFICATION_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SPECIFICATION_SET_NUMBER(USER_NUMBER,REGION,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SPECIFICATION_SET_NUMBER
    !###  Description:
    !###    Sets/changes the problem specification i.e., problem class, type  and subtype for a problem identified by a user
    !###    number.
    !###  Parent-subroutine: PROBLEM_SPECIFICATION_SET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(IN) :: PROBLEM_CLASS
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("PROBLEM_SPECIFICATION_SET_NUMBER",ERR,ERROR,*999)

!!TODO: Take in region number here and user FIND_REGION_NUMBER. This would require FIND_REGION_NUMBER to be moved from
!!REGION_ROUTINES otherwise there will be a circular module reference.

    CALL PROBLEM_USER_NUMBER_FIND(USER_NUMBER,REGION,PROBLEM,ERR,ERROR,*999)
    CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
           
    CALL EXITS("PROBLEM_SPECIFICATION_SET_NUMBER")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_SET_NUMBER",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_SET_NUMBER")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_SET_NUMBER
  
  !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_SPECIFICATION_SET_PTR(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_SPECIFICATION_SET_PTR
    !###  Description:
    !###    Sets/changes the problem specification i.e., problem class, type and subtype for a problem identified by a pointer.
    !###  Parent-subroutine: PROBLEM_SPECIFICATION_SET

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(IN) :: PROBLEM_CLASS
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SPECIFICATION_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        CALL FLAG_ERROR("Problem has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(PROBLEM_CLASS)
        CASE(PROBLEM_ELASTICITY_CLASS)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PROBLEM_FLUID_MECHANICS_CLASS)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
          CALL PROBLEM_CLASSICAL_FIELD_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM_CLASS,"*",ERR,ERROR))//" is not valid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SPECIFICATION_SET_PTR")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_SET_PTR",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_SET_PTR")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_SET_PTR
  
  !
  !================================================================================================================================
  !
  
  !> Finalises the time data information for a problem and deallocates all memory
  SUBROUTINE PROBLEM_TIME_DATA_FINALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_TIME_DATA_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%TIME_DATA)) THEN
        DEALLOCATE(PROBLEM_SOLUTION%TIME_DATA)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_TIME_DATA_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_TIME_DATA_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_TIME_DATA_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_TIME_DATA_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the time data information for a problem solution
  SUBROUTINE PROBLEM_TIME_DATA_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<The pointer to the problem solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_TIME_DATA_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%TIME_DATA)) THEN
        CALL FLAG_ERROR("Time data is already associated for this problem solution",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLUTION%TIME_DATA,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution time data",ERR,ERROR,*999)
        PROBLEM_SOLUTION%TIME_DATA%PROBLEM_SOLUTION=>PROBLEM_SOLUTION        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_TIME_DATA_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_TIME_DATA_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_TIME_DATA_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_TIME_DATA_INITIALISE

   !
  !================================================================================================================================
  !

  SUBROUTINE PROBLEM_USER_NUMBER_FIND(USER_NUMBER,REGION,PROBLEM,ERR,ERROR,*)

    !#### Subroutine: PROBLEM_USER_NUMBER_FIND
    !###  Description:
    !###    Finds and returns in PROBLEM a pointer to the problem identified by USER_NUMBER in the given REGION. If no problem
    !##     with that USER_NUMBER exists PROBLEM is left nullified.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_USER_NUMBER_FIND",ERR,ERROR,*999)

    NULLIFY(PROBLEM)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%PROBLEMS)) THEN
        problem_idx=1
        DO WHILE(problem_idx<=REGION%PROBLEMS%NUMBER_OF_PROBLEMS.AND..NOT.ASSOCIATED(PROBLEM))
          IF(REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            PROBLEM=>REGION%PROBLEMS%PROBLEMS(problem_idx)%PTR
          ELSE
            problem_idx=problem_idx+1
          ENDIF
        ENDDO
      ELSE
        LOCAL_ERROR="The problems on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("PROBLEM_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("PROBLEM_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE PROBLEM_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises all problems on a region and deallocates all memory.
  SUBROUTINE PROBLEMS_FINALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finalise the problems for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::USER_NUMBER

    CALL ENTERS("PROBLEMS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%PROBLEMS)) THEN
        DO WHILE(REGION%PROBLEMS%NUMBER_OF_PROBLEMS>0)
          USER_NUMBER=REGION%PROBLEMS%PROBLEMS(1)%PTR%USER_NUMBER
          CALL PROBLEM_DESTROY(USER_NUMBER,REGION,ERR,ERROR,*999)
        ENDDO !problem_idx
        DEALLOCATE(REGION%PROBLEMS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("PROBLEMS_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEMS_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEMS_FINALISE")
    RETURN 1   
  END SUBROUTINE PROBLEMS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises all problems on a region.
  SUBROUTINE PROBLEMS_INITIALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the problems for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEMS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%PROBLEMS)) THEN
        CALL FLAG_ERROR("Region already has associated problems",ERR,ERROR,*998)
      ELSE
!!TODO: Inherit any problems from the parent region???
        ALLOCATE(REGION%PROBLEMS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region problems",ERR,ERROR,*999)
        REGION%PROBLEMS%REGION=>REGION
        REGION%PROBLEMS%NUMBER_OF_PROBLEMS=0
        NULLIFY(REGION%PROBLEMS%PROBLEMS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associted",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("PROBLEMS_INITIALISE")
    RETURN
999 IF(ASSOCIATED(REGION%PROBLEMS)) DEALLOCATE(REGION%PROBLEMS)
998 CALL ERRORS("PROBLEMS_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEMS_INITIALISE")
    RETURN 1   
  END SUBROUTINE PROBLEMS_INITIALISE
  
  !
  !================================================================================================================================
  !
  
END MODULE PROBLEM_ROUTINES
