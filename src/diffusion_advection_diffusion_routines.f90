!> \file
!> \authors Andrew Cookson
!> \brief This module handles all routines pertaining to diffusion coupled to diffusion.
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
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
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

!>TThis module handles all routines pertaining to diffusion coupled to diffusion.


MODULE DIFFUSION_ADVECTION_DIFFUSION_ROUTINES

  USE ADVECTION_DIFFUSION_EQUATION_ROUTINES
  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES  
  USE DIFFUSION_EQUATION_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
!  USE FINITE_ELASTICITY_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
!   USE FITTING_ROUTINES !also in makefiles
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS  
  USE MATRIX_VECTOR
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES


  IMPLICIT NONE

  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_EQUATIONS_SET_SETUP
  PUBLIC DiffusionAdvecDiff_EquationsSetSpecificationSet
  PUBLIC DIFFUSION_ADVEC_DIFF_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP
  PUBLIC DiffusionAdvectionDiffusion_ProblemSpecificationSet
  
  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_FINITE_ELEMENT_CALCULATE

  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE
  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE

  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE DIFFUSION_ADVEC_DIFF_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_ADVEC_DIFF_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a "// &
          & "diffusion and advection-diffusion type equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a diffusion & advection-diffusion equation type of a multi physics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_ADVEC_DIFF_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("DIFFUSION_ADVEC_DIFF_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("DIFFUSION_ADVEC_DIFF_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVEC_DIFF_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion & advection-diffusion coupled equation.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables


    CALL ENTERS("DIFFUSION_ADVECTION_DIFFUSION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
             
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_ADVECTION_DIFFUSION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_EQUATIONS_SET_SETUP")
    RETURN 1

  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a coupled diffusion & advection-diffusion equation finite element equations set.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables



    CALL ENTERS("DIFFUSION_ADVECTION_DIFFUSION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      
  CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("DIFFUSION_ADVECTION_DIFFUSION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE DiffusionAdvecDiff_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("DiffusionAdvecDiff_EquationsSetSpecificationSet",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)

    CALL Exits("DiffusionAdvecDiff_EquationsSetSpecificationSet")
    RETURN
999 CALL Errors("DiffusionAdvecDiff_EquationsSetSpecificationSet",err,error)
    CALL Exits("DiffusionAdvecDiff_EquationsSetSpecificationSet")
    RETURN 1
  END SUBROUTINE DiffusionAdvecDiff_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a coupled diffusion & advection-diffusion problem.
  SUBROUTINE DiffusionAdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    CALL Enters("DiffusionAdvectionDiffusion_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a coupled diffusion & advection-diffusion equation type of a multi physics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("Diffusion advection-diffusion transport problem specification must have 3 entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    CALL Exits("DiffusionAdvectionDiffusion_ProblemSpecificationSet")
    RETURN
999 CALL Errors("DiffusionAdvectionDiffusion_ProblemSpecificationSet",err,error)
    CALL Exits("DiffusionAdvectionDiffusion_ProblemSpecificationSet")
    RETURN 1
  END SUBROUTINE DiffusionAdvectionDiffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the coupled diffusion-diffusion equations problem.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION, SOLVER_ADVECTION_DIFFUSION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION, SOLVER_EQUATIONS_ADVECTION_DIFFUSION
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER_DIFFUSION)
    NULLIFY(SOLVER_ADVECTION_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
          & err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   coupled source diffusion--advection-diffusion
      !--------------------------------------------------------------------
      CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a diffusion & advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled diffusion & advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,ERR,ERROR,*999)
            !
            !Set the first solver to be a linear solver for the advection-diffusion problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !
            !Set the second solver to be a linear solver for the diffusion problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_DIFFUSION,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a coupled diffusion & advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !
            !Get the advection-diffusion solver and create the advection-diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION, & 
              & SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            !Get the diffusion solver and create the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_DIFFUSION,SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_DIFFUSION,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_DIFFUSION, & 
              & SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_DIFFUSION,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !

          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the advection-diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)             
            !
            !Finish the creation of the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_DIFFUSION,SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)             
            !
          
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled diffusion & advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a coupled diffusion & advection-diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   c a s e   d e f a u l t
      !-----------------------------------------------------------------
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a coupled source diffusion & advection-diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)

      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion-diffusion problem pre-solve.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR


    CALL ENTERS("DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)

              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !copy current value of concentration_one to another variable
                CALL ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !Set source term to be updated value of concentration_two
                !CALL ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                
                !Update indpendent data fields
!                 CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
!                 CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
 
              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                !copy current value of concentration_one to another variable
                CALL DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !compute value of constant source term - evaluated from lamdba*(0.5*(c_1^{t+1}+c_1^{t}) - c_2^{t})
                !CALL DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE")
    RETURN
999 CALL ERRORS("DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the diffusion-diffusion problem post solve.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
!              CALL ADVECTION_DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
!              CALL DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE")
    RETURN
999 CALL ERRORS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE",ERR,ERROR)
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE")
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the diffuion-diffusion problem post solve output data.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                !CALL ADVECTION_DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !CALL DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 CALL ERRORS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    CALL EXITS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA")
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA
      
  !   
  !================================================================================================================================
  !


END MODULE DIFFUSION_ADVECTION_DIFFUSION_ROUTINES
