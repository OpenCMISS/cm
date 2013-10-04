!> \file
!> \authors Adam Reeve
!> \brief This module handles all routines pertaining to finite elasticity coupled with fluid pressure for poroelasticity.
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

!>This module handles all routines pertaining to finite elasticity coupled with fluid pressure for poroelasticity.


MODULE FINITE_ELASTICITY_FLUID_PRESSURE_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TYPES

  IMPLICIT NONE

  PUBLIC ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP
  PUBLIC ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP
  PUBLIC ElasticityFluidPressure_ProblemSpecificationSet
  
  PUBLIC ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE

  PUBLIC ELASTICITY_FLUID_PRESSURE_PRE_SOLVE
  PUBLIC ELASTICITY_FLUID_PRESSURE_POST_SOLVE

  PUBLIC ELASTICITY_FLUID_PRESSURE_CONTROL_LOOP_PRE_LOOP

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity fluid pressure equation type of a multi physics equations set class.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a "// &
          & "finite elasticity-fluid pressure class equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
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
          & " is not valid for a finite elasticity fluid pressure equation type of a multi physics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure equation.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string


    CALL ENTERS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    CALL FLAG_ERROR("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP is not implemented.",ERR,ERROR,*999)

    CALL EXITS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a finite elasticity fluid pressure equation finite element equations set.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    CALL FLAG_ERROR("ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE is not implemented.",ERR,ERROR,*999)

    CALL EXITS("ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a finite elasticity fluid pressure equation type.
  SUBROUTINE ElasticityFluidPressure_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    CALL Enters("ElasticityFluidPressure_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a finite elasticity fluid pressure equation type of a multi physics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("Finite elasticity fluid pressure problem specificaion must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    CALL Exits("ElasticityFluidPressure_ProblemSpecificationSet")
    RETURN
999 CALL Errors("ElasticityFluidPressure_ProblemSpecificationSet",err,error)
    CALL Exits("ElasticityFluidPressure_ProblemSpecificationSet")
    RETURN 1
  END SUBROUTINE ElasticityFluidPressure_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure equations problem.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ALLOCATED(problem%specification)) THEN
        IF(.NOT.ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(problem%specification,1)<3) THEN
          CALL FlagError("Problem specification must have three entries for a finite elasticity-Darcy problem.", &
            & err,error,*999)
        END IF
      ELSE
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   Standard finite elasticity fluid pressure
      !--------------------------------------------------------------------
      CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for an finite elasticity ALE fluid pressure  equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a load increment loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity fluid pressure equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation for the solver
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !
            !Set the first solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity fluid pressure equation."
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
            !Get the nonlinear solver and create the solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity fluid pressure equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity ALE fluid pressure equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a standard finite elasticity fluid pressure equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)

      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity fluid pressure problem pre-solve.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_FLUID_PRESSURE_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for an elasticity fluid pressure problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
            IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
              CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a fluid pressure fluid type of a multi physics problem class."
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

    CALL EXITS("ELASTICITY_FLUID_PRESSURE_PRE_SOLVE")
    RETURN
999 CALL ERRORS("ELASTICITY_FLUID_PRESSURE_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure  problem post solve.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_FLUID_PRESSURE_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for an elasticity fluid pressure problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
            CALL FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a finite elasticity fluid pressure type of a multi physics problem class."
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

    CALL EXITS("ELASTICITY_FLUID_PRESSURE_POST_SOLVE")
    RETURN
999 CALL ERRORS("ELASTICITY_FLUID_PRESSURE_POST_SOLVE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_POST_SOLVE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_FLUID_PRESSURE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_FLUID_PRESSURE

    CALL ENTERS("ELASTICITY_FLUID_PRESSURE_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_FLUID_PRESSURE)
    NULLIFY(SOLVER_FLUID_PRESSURE)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        ! Eventually we may want to do different things depending on problem type/subtype
        ! too but for now we can just check the loop type.
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"-- Starting load increment --",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LOAD INCREMENT NUMBER =           ", &
              & CONTROL_LOOP%LOAD_INCREMENT_LOOP%ITERATION_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
          ENDIF
          IF(DIAGNOSTICS1) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"-- Starting load increment --",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"LOAD INCREMENT NUMBER =           ", &
              & CONTROL_LOOP%LOAD_INCREMENT_LOOP%ITERATION_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
          ENDIF

        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ELASTICITY_FLUID_PRESSURE_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("ELASTICITY_FLUID_PRESSURE_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("ELASTICITY_FLUID_PRESSURE_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

END MODULE FINITE_ELASTICITY_FLUID_PRESSURE_ROUTINES
