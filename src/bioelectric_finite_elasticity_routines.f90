!> \file
!> \authors Thomas Heidlauf
!> \brief This module handles all routines pertaining to bioelectrics coupled with finite elasticity.
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
!> Contributor(s): Thomas Klotz
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

!>This module handles all routines pertaining to bioelectrics coupled with finite elasticity.


MODULE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BIOELECTRIC_ROUTINES
  USE BIODOMAIN_EQUATION_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TYPES

#include "macros.h"

  IMPLICIT NONE

  PUBLIC BioelectricFiniteElasticity_EquationsSetSetup
  PUBLIC BioelectricFiniteElasticity_EquationsSetSolutionMethodSet

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP
  PUBLIC BioelectricFiniteElasticity_ProblemSpecificationSet
 
  PUBLIC BioelectricFiniteElasticity_FiniteElementCalculate

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  PUBLIC BioelectricFiniteElasticity_ControlLoopPreLoop
  PUBLIC BioelectricFiniteElasticity_ControlLoopPostLoop
  
  PUBLIC BioelectricFiniteElasticity_UpdateGeometricField

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectrics finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE BioelectricFiniteElasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a "// &
          & "Bioelectric-finite elasticity type equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        !\todo what are problem constants doing here???
      CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a bioelectrics finite elasticity equation type of a multi physics equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity equation.
  SUBROUTINE BioelectricFiniteElasticity_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("BioelectricFiniteElasticity_EquationsSetSetup",ERR,ERROR,*999)

    CALL FlagError("BioelectricFiniteElasticity_EquationsSetSetup is not implemented.",ERR,ERROR,*999)

    EXITS("BioelectricFiniteElasticity_EquationsSetSetup")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_EquationsSetSetup",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a bioelectrics finite elasticity equation finite element equations set.
  SUBROUTINE BioelectricFiniteElasticity_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("BioelectricFiniteElasticity_FiniteElementCalculate",ERR,ERROR,*999)

    CALL FlagError("BioelectricFiniteElasticity_FiniteElementCalculate is not implemented.",ERR,ERROR,*999)

    EXITS("BioelectricFiniteElasticity_FiniteElementCalculate")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_FiniteElementCalculate",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a bioelectric finite elasticity problem type .
  SUBROUTINE BioelectricFiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("BioelectricFiniteElasticity_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE, &
            & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a bioelectric finite elasticity type of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("Bioelectric finite elasticity problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("BioelectricFiniteElasticity_ProblemSpecificationSet")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ProblemSpecificationSet",err,error)
    EXITS("BioelectricFiniteElasticity_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectric finite elasticity problem.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: MONODOMAIN_SUB_LOOP,ELASTICITY_SUB_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS,MONODOMAIN_SOLVERS,ELASTICITY_SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(MONODOMAIN_SUB_LOOP)
    NULLIFY(ELASTICITY_SUB_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(MONODOMAIN_SOLVERS)
    NULLIFY(ELASTICITY_SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(CELLML_EQUATIONS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   Transient Gudunov monodomain, simple finite elasticity  
      !--------------------------------------------------------------------
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
        & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
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
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up the control sub loop for monodomain
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(MONODOMAIN_SUB_LOOP,'MONODOMAIN_TIME_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(MONODOMAIN_SUB_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(MONODOMAIN_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            IF(PROBLEM%specification(3)==PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE) THEN
              !Set up the control sub loop for finite elasicity
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
              CALL CONTROL_LOOP_LABEL_SET(ELASTICITY_SUB_LOOP,'ELASTICITY_WHILE_LOOP',ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(ELASTICITY_SUB_LOOP,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_OUTPUT_TYPE_SET(ELASTICITY_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
            ELSE
              !Set up the control sub loop for finite elasicity
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
              CALL CONTROL_LOOP_LABEL_SET(ELASTICITY_SUB_LOOP,'ELASTICITY_LOAD_INCREMENT_LOOP',ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(ELASTICITY_SUB_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_OUTPUT_TYPE_SET(ELASTICITY_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
            ENDIF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
            !Sub-loops are finished when parent is finished
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the monodomain sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(MONODOMAIN_SOLVERS,2,ERR,ERROR,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"ODE Solver",ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)

            !Get the finite elasticity sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(ELASTICITY_SOLVERS,1,ERR,ERROR,*999)
            !Set the finite elasticity solver to be a nonlinear solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the monodomain solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(MONODOMAIN_SOLVERS,ERR,ERROR,*999)

            !Get the finite elasticity solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(ELASTICITY_SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)

            !Get the monodomain sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Create the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)

            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Get the finite elasticity solver and create the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            
            !Get the monodomain sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Get the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             

            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the creation of the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Get the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectrics finite elasticity equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a transient monodomain quasistatic finite elasticity equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the bioelectrics finite elasticity problem pre-solve.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
              & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              CALL BIODOMAIN_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              CALL BIODOMAIN_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
              CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem post solve.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
            SELECT CASE(SOLVER%SOLVE_TYPE)
            CASE(SOLVER_DAE_TYPE)
              CALL BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_TYPE)
              CALL BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_NONLINEAR_TYPE)
              CALL FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Solver solve type "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem pre-control loop.
  SUBROUTINE BioelectricFiniteElasticity_ControlLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BioelectricFiniteElasticity_ControlLoopPreLoop",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              !do nothing ???
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)
                CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*999)
              CASE(PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
                CALL BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN(CONTROL_LOOP,ERR,ERROR,*999)
                CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The third problem specification of "// &
                  & TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for bioelectric finite elasticity problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL FiniteElasticity_ControlTimeLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
                IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1) THEN
                  CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*999)
                ENDIF
                CALL BioelectricFiniteElasticity_ComputeFibreStretch(CONTROL_LOOP,ERR,ERROR,*999)
                CALL BioelectricFiniteElasticity_ForceLengthVelocityRelation(CONTROL_LOOP,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The third problem specification of "// &
                  & TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for bioelectric finite elasticity problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL FiniteElasticity_ControlTimeLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for bioelectric finite elasticity problem type."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The second problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("BioelectricFiniteElasticity_ControlLoopPreLoop")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ControlLoopPreLoop",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ControlLoopPreLoop")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ControlLoopPreLoop

  !
  !================================================================================================================================
  !

  !>Computes the fibre stretch for bioelectrics finite elasticity problems.
  SUBROUTINE BioelectricFiniteElasticity_ComputeFibreStretch(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(BASIS_TYPE), POINTER :: GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_U1
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_ELEMENTS,FIELD_VAR_TYPE
    INTEGER(INTG) :: equations_set_idx,gauss_idx,dof_idx,element_idx,idx
    REAL(DP) :: DZDNU(3,3),DZDNUT(3,3),AZL(3,3)

    ENTERS("BioelectricFiniteElasticity_ComputeFibreStretch",ERR,ERROR,*999)

    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(FIELD_VAR_U1)
    NULLIFY(DEPENDENT_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
          CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
          CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
          !Loop over the equations sets associated with the solver
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) THEN
                    LOCAL_ERROR="Geometric field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(DEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Dependent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Independent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  EQUATIONS=>EQUATIONS_SET%EQUATIONS
                  IF(.NOT.ASSOCIATED(EQUATIONS)) THEN
                    LOCAL_ERROR="Equations is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  FIBRE_FIELD=>EQUATIONS%INTERPOLATION%FIBRE_FIELD
                  IF(.NOT.ASSOCIATED(FIBRE_FIELD)) THEN
                    LOCAL_ERROR="Fibre field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Equations set is not associated for equations set index "// &
                    & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !equations_set_idx
            ELSE
              CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
          ENDIF

          CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_U1,ERR,ERROR,*999)

          DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER

          IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1) THEN
            !copy the old fibre stretch to the previous values parameter set
            CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
          ENDIF

          NUMBER_OF_ELEMENTS=DEPENDENT_FIELD%GEOMETRIC_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS

          !loop over the elements of the finite elasticity mesh (internal and boundary elements)
          DO element_idx=1,NUMBER_OF_ELEMENTS

            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BASIS       
            DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BASIS

            !Initialise tensors and matrices
            DZDNU=0.0_DP
            DO idx=1,3
              DZDNU(idx,idx)=1.0_DP
            ENDDO

            !Grab interpolation parameters
            FIELD_VAR_TYPE=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%VARIABLE_TYPE
            DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR
            GEOMETRIC_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
            FIBRE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR

            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
              & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
              & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
              & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)

            !Point interpolation pointer
            GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            GEOMETRIC_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
            FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
            DEPENDENT_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_TYPE)%PTR

            !Loop over gauss points
            DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            
              !Interpolate dependent, geometric, fibre and materials fields
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI, &
                & DEPENDENT_INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI, &
                & GEOMETRIC_INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)

              !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
              CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
                & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)
              
              !compute C=F^T F
              CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)
              CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999)

              !store the fibre stretch lambda_f, i.e., sqrt(C_11) or AZL(1,1)
              dof_idx=FIELD_VAR_U1%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & dof_idx,SQRT(AZL(1,1)),ERR,ERROR,*999)

            ENDDO !gauss_idx
          ENDDO !element_idx

          !now the ghost elements -- get the relevant info from the other computational nodes
          CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

        CASE DEFAULT
          LOCAL_ERROR="Control loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
            & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE !the control loop contains subloops
        !do nothing
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_ComputeFibreStretch")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ComputeFibreStretch",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ComputeFibreStretch")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ComputeFibreStretch

  !
  !================================================================================================================================
  !

  !>Computes the bioelectrics finite elasticity force-length and force_velocity relations.
  SUBROUTINE BioelectricFiniteElasticity_ForceLengthVelocityRelation(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_PARENT
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(BASIS_TYPE), POINTER :: GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_U,FIELD_VAR_U1
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_ELEMENTS
    INTEGER(INTG) :: equations_set_idx,gauss_idx,dof_idx,element_idx
    INTEGER(INTG) :: ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS
    REAL(DP) :: LENGTH_HS,LENGTH_HS_0,ACTIVE_STRESS,FIBRE_STRETCH,FIBRE_STRETCH_OLD
    REAL(DP) :: FACTOR_LENGTH,FACTOR_VELO,SARCO_LENGTH,VELOCITY,VELOCITY_MAX,TIME_STEP,kappa,A,S,d,c

!tomo
    REAL(DP) :: VELOCITY_AVERAGE,STRETCH_AVERAGE,OLD_STRETCH_AVERAGE
    INTEGER(INTG) :: counter

    ENTERS("BioelectricFiniteElasticity_ForceLengthVelocityRelation",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(GEOMETRIC_BASIS)
    NULLIFY(DECOMPOSITION)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(FIELD_VAR_U,FIELD_VAR_U1)
    NULLIFY(DEPENDENT_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP%PROBLEM%CONTROL_LOOP,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
          CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
          CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
          !Loop over the equations sets associated with the solver
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(DEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Dependent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Independent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Equations set is not associated for equations set index "// &
                    & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !equations_set_idx
            ELSE
              CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
          ENDIF

          CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VAR_U,ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_U1,ERR,ERROR,*999)

          !get the inital half-sarcomere length
          dof_idx=FIELD_VAR_U1%COMPONENTS(2)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,LENGTH_HS_0,ERR,ERROR,*999)

          !get the maximum contraction velocity
          dof_idx=FIELD_VAR_U1%COMPONENTS(3)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VELOCITY_MAX,ERR,ERROR,*999)

          ITERATION_NUMBER=CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
          MAXIMUM_NUMBER_OF_ITERATIONS=CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS
          !in the first iteration store the unaltered homogenized active stress field 
          IF(ITERATION_NUMBER==1) THEN
            CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
          ELSE
            !restore the solution of the previous time step
            CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
          ENDIF

          TIME_STEP=CONTROL_LOOP_PARENT%TIME_LOOP%TIME_INCREMENT

          DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER

!tomo start
          VELOCITY_AVERAGE=0.0_DP
          STRETCH_AVERAGE=0.0_DP
          OLD_STRETCH_AVERAGE=0.0_DP
          counter=0
!tomo end

          NUMBER_OF_ELEMENTS=DEPENDENT_FIELD%GEOMETRIC_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS

          !loop over the elements of the finite elasticity mesh (internal and boundary elements)
          DO element_idx=1,NUMBER_OF_ELEMENTS

            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BASIS       
            DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

            !Loop over gauss points
            DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS

              !get the unaltered ACTIVE_STRESS at the GP
              dof_idx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                & dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)

              ! FORCE-LENGTH RELATION -------------------------------------------------------------------------------
              
              !get the current fibre stretch at the GP
              dof_idx=FIELD_VAR_U1%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,dof_idx,FIBRE_STRETCH,ERR,ERROR,*999)

              !compute the current half-sarcomere length at the GP: l_hs = lambda_f * l_hs_0
              LENGTH_HS=LENGTH_HS_0*FIBRE_STRETCH
              
              !compute the scale factor (0,1) due to sarcomere F-l relation of Gordon, A. M., A.F. Huxley, and F.J. Julian. 
              !The variation in isometric tension with sarcomere length in vertebrate muscle fibres. 
              !The Journal of Physiology 184.1 (1966): 170-192.
              SARCO_LENGTH=2.0_DP*LENGTH_HS
              IF(SARCO_LENGTH.LE.1.27_DP) THEN
                FACTOR_LENGTH=0.0_DP
              ELSEIF(SARCO_LENGTH.LE.1.7_DP) THEN
                FACTOR_LENGTH=1.6047_DP*SARCO_LENGTH-2.0379_DP
              ELSEIF(SARCO_LENGTH.LE.2.34_DP) THEN
                FACTOR_LENGTH=0.4844_DP*SARCO_LENGTH-0.1334_DP
              ELSEIF(SARCO_LENGTH.LE.2.51_DP) THEN
                FACTOR_LENGTH=1.0_DP
              ELSEIF(SARCO_LENGTH.LE.3.94_DP) THEN
                FACTOR_LENGTH=-0.6993_DP*SARCO_LENGTH+2.7552_DP
              ELSE
                FACTOR_LENGTH=0.0_DP
              ENDIF

              !multiply the ACTIVE_STRESS with the scale factor
              ACTIVE_STRESS=ACTIVE_STRESS*FACTOR_LENGTH

              !update the ACTIVE_STRESS at GP
              dof_idx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)





              ! FORCE-VELOCITY RELATION -------------------------------------------------------------------------------
              
              !get fibre stretch at the GP of the previous time step
              dof_idx=FIELD_VAR_U1%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_VALUES_SET_TYPE,dof_idx,FIBRE_STRETCH_OLD,ERR,ERROR,*999)

              !compute the contraction velocity
              VELOCITY=(FIBRE_STRETCH-FIBRE_STRETCH_OLD)/TIME_STEP

!!!! original
              !NOTE: VELOCITY_MAX is the max shortening velocity, and hence negative
              IF(VELOCITY<VELOCITY_MAX) THEN
                CALL FLAG_WARNING('Exceeded maximum contraction velocity (shortening).',ERR,ERROR,*999)
!                VELOCITY=VELOCITY_MAX
                !damping
                IF(ITERATION_NUMBER<(MAXIMUM_NUMBER_OF_ITERATIONS/2)) THEN
                  VELOCITY=VELOCITY*1.0_DP/DBLE((MAXIMUM_NUMBER_OF_ITERATIONS/2)-ITERATION_NUMBER)
                ENDIF
              ELSEIF(VELOCITY>(ABS(VELOCITY_MAX))) THEN
                CALL FLAG_WARNING('Exceeded maximum contraction velocity (lengthening).',ERR,ERROR,*999)
!!!                VELOCITY=-VELOCITY_MAX
!                !damping
!                IF(ITERATION_NUMBER<(MAXIMUM_NUMBER_OF_ITERATIONS/2)) THEN
!                  VELOCITY=VELOCITY*1.0_DP/DBLE((MAXIMUM_NUMBER_OF_ITERATIONS/2)-ITERATION_NUMBER)
!                ENDIF
              ENDIF

              !compute scale factor (FACTOR_VELO)
              kappa=0.24_DP !make mat param TODO
              A=1.6_DP !make mat param TODO
              S=2.5_DP !make mat param TODO
              IF(VELOCITY<0.0_DP) THEN
                !shortening contraction              
                FACTOR_VELO=(1-VELOCITY/VELOCITY_MAX)/(1+VELOCITY/VELOCITY_MAX/kappa)
              ELSE
                !lengthening contraction
                d=kappa*(1-A)
                c=VELOCITY/VELOCITY_MAX*S*(kappa+1)
                FACTOR_VELO=1+c*(A-1)/(d+c)
              ENDIF
              
              !multiply the ACTIVE_STRESS with the scale factor
              ACTIVE_STRESS=ACTIVE_STRESS*FACTOR_VELO

              !update the ACTIVE_STRESS at GP
              dof_idx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)

            ENDDO !gauss_idx
          ENDDO !element_idx
!!!!end original




!!!!tomo --- try averaging              
!!!              VELOCITY_AVERAGE=VELOCITY_AVERAGE+VELOCITY
!!!              STRETCH_AVERAGE=STRETCH_AVERAGE+FIBRE_STRETCH
!!!              OLD_STRETCH_AVERAGE=OLD_STRETCH_AVERAGE+FIBRE_STRETCH_OLD
!!!              
!!!              counter=counter+1

!!!            ENDDO !gauss_idx
!!!          ENDDO !element_idx

!!!          VELOCITY=VELOCITY_AVERAGE/REAL(counter)
!!!          FIBRE_STRETCH=STRETCH_AVERAGE/REAL(counter)
!!!          FIBRE_STRETCH_OLD=OLD_STRETCH_AVERAGE/REAL(counter)

!!!          VELOCITY_AVERAGE=(FIBRE_STRETCH-FIBRE_STRETCH_OLD)/TIME_STEP

!!!!          VELOCITY=0.0_DP

!!!!          !damping
!!!!          IF(ITERATION_NUMBER<(MAXIMUM_NUMBER_OF_ITERATIONS/2)) THEN
!!!!            VELOCITY=VELOCITY*1.0_DP/DBLE((MAXIMUM_NUMBER_OF_ITERATIONS/2)-ITERATION_NUMBER)
!!!!          ENDIF

!!!          LOCAL_ERROR="######### VELOCITY: "//TRIM(NUMBER_TO_VSTRING(VELOCITY,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)
!!!          LOCAL_ERROR="######### VELOCITY AVERAGE: "//TRIM(NUMBER_TO_VSTRING(VELOCITY_AVERAGE,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)
!!!          LOCAL_ERROR="######### STRETCH: "//TRIM(NUMBER_TO_VSTRING(FIBRE_STRETCH,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)
!!!          LOCAL_ERROR="######### OLD STRETCH: "//TRIM(NUMBER_TO_VSTRING(FIBRE_STRETCH_OLD,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)


!!!          
!!!          
!!!          VELOCITY=VELOCITY_AVERAGE
!!!          
!!!          
!!!          

!!!          !compute scale factor (FACTOR_VELO)
!!!          kappa=0.24_DP !make mat param TODO
!!!          A=1.6_DP !make mat param TODO
!!!          S=2.5_DP !make mat param TODO
!!!          IF(VELOCITY<0.0_DP) THEN
!!!            !shortening contraction              
!!!            FACTOR_VELO=(1-VELOCITY/VELOCITY_MAX)/(1+VELOCITY/VELOCITY_MAX/kappa)
!!!            IF(VELOCITY<VELOCITY_MAX) THEN
!!!              CALL FLAG_WARNING('Exceeded maximum contraction velocity (shortening).',ERR,ERROR,*999)
!!!            ENDIF
!!!          ELSE
!!!            !lengthening contraction
!!!            d=kappa*(1-A)
!!!            c=VELOCITY/VELOCITY_MAX*S*(kappa+1)
!!!            FACTOR_VELO=1+c*(A-1)/(d+c)
!!!            IF(VELOCITY>ABS(VELOCITY_MAX)) THEN
!!!              CALL FLAG_WARNING('Exceeded maximum contraction velocity (lengthening).',ERR,ERROR,*999)
!!!            ENDIF
!!!          ENDIF

!!!          DO element_idx=1,NUMBER_OF_ELEMENTS

!!!            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BASIS       
!!!            DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
!!!            DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

!!!            !Loop over gauss points
!!!            DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS

!!!              !get the ACTIVE_STRESS at the GP
!!!              dof_idx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
!!!              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!!!                & dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)

!!!              !multiply the ACTIVE_STRESS with the scale factor
!!!              ACTIVE_STRESS=ACTIVE_STRESS*FACTOR_VELO

!!!              !update the ACTIVE_STRESS at GP
!!!              dof_idx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,element_idx)
!!!              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!!!                & dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)
!!!              
!!!            ENDDO !gauss_idx
!!!          ENDDO !element_idx
!tomo end

          !now the ghost elements -- get the relevant info from the other computational nodes
          CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

        CASE DEFAULT
          LOCAL_ERROR="Control loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
            & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE !the control loop contains subloops
        !do nothing
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF


    EXITS("BioelectricFiniteElasticity_ForceLengthVelocityRelation")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ForceLengthVelocityRelation",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ForceLengthVelocityRelation")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ForceLengthVelocityRelation

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem post-control loop.
  SUBROUTINE BioelectricFiniteElasticity_ControlLoopPostLoop(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG) :: equations_set_idx
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION   
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,LOCAL_ERROR,METHOD
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ELASTICITY_SUB_LOOP,BIOELECTRIC_SUB_LOOP

    ENTERS("BioelectricFiniteElasticity_ControlLoopPostLoop",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
              !the monodomain time loop - output of the monodomain fields
              CALL BIODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
                & " is not valid for a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            CALL BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL BioelectricFiniteElasticity_ConvergenceCheck(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            !do nothing
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - output the finite elasticity fields 
        IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          !Export the dependent field for this time step
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            PROBLEM=>CONTROL_LOOP%PROBLEM
            IF(ASSOCIATED(PROBLEM)) THEN
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(ELASTICITY_SUB_LOOP)
              !Get the solver. The first solver of the second sub loop will contain the finite elasticity dependent field equation set
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              !Loop over the equations sets associated with the solver
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      NULLIFY(DEPENDENT_REGION)
                      CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                      FILENAME="MainTime_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                        & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                      METHOD="FORTRAN"
                      CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for equations set index "// &
                        & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                        & " in the solver mapping."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                ELSE
                  CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
              ENDIF
              IF((PROBLEM%SPECIFICATION(3)==PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE).OR. &
               & (PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE).OR. &
               & (PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)) THEN
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(BIOELECTRIC_SUB_LOOP)
                !Get the solver. The second solver of the first sub loop will contain the bioelectrics equation set
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,BIOELECTRIC_SUB_LOOP,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(BIOELECTRIC_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                !Loop over the equations sets associated with the solver
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        NULLIFY(DEPENDENT_REGION)
                        CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                        FILENAME="MainTime_M_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                          & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                        METHOD="FORTRAN"
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                        
                        WRITE(*,*) TIME_LOOP%ITERATION_NUMBER
                        
                      ELSE
                        LOCAL_ERROR="Equations set is not associated for equations set index "// &
                          & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                          & " in the solver mapping."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !equations_set_idx
                  ELSE
                    CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF !PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE
            ELSE
              CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_ControlLoopPostLoop")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ControlLoopPostLoop",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ControlLoopPostLoop")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ControlLoopPostLoop


  !
  !================================================================================================================================
  !

  !>Check for the convergence of the bioelectric finite elasticity while loop, i.e., if the force-length and force-velocity relations yielded a different actual configuration.
  SUBROUTINE BioelectricFiniteElasticity_ConvergenceCheck(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG) :: equations_set_idx,NUMBER_OF_NODES,dof_idx,node_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR
    REAL(DP) :: x1,x2,x3,y1,y2,y3,my_sum

    ENTERS("BioelectricFiniteElasticity_ConvergenceCheck",ERR,ERROR,*999)

    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(FIELD_VAR)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            !do nothing
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            !do nothing
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLUTION_UPDATE(SOLVER,ERR,ERROR,*999) !tomo added this
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Loop over the equations sets associated with the solver
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  ELSE
                    LOCAL_ERROR="Equations set is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                      & " in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !equations_set_idx
              ELSE
                CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
            ENDIF

            IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1) THEN
              !
            ELSE
              CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VAR,ERR,ERROR,*999)

              NUMBER_OF_NODES=DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION% &
                & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES%NUMBER_OF_LOCAL

              my_sum=0.0_DP
              DO node_idx=1,NUMBER_OF_NODES

                !get the current node position (x) and the node position of the last iteration (y)
                dof_idx=FIELD_VAR%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)% &
                  & DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,dof_idx,x1,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,dof_idx,y1,ERR,ERROR,*999)
                dof_idx=FIELD_VAR%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)% &
                  & DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,dof_idx,x2,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,dof_idx,y2,ERR,ERROR,*999)
                dof_idx=FIELD_VAR%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)% &
                  & DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,dof_idx,x3,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,dof_idx,y3,ERR,ERROR,*999)

                !sum up
                my_sum=my_sum+SQRT((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
              ENDDO
              
              IF(my_sum<1.0E-06_DP) THEN !if converged then:
                CONTROL_LOOP%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
                CALL BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
                !copy the current solution to the previous solution
                CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
              ELSEIF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS) THEN
                CALL BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
                CALL FLAG_WARNING('----------- Maximum number of iterations in while loop reached. -----------',ERR,ERROR,*999)
                !copy the current solution to the previous solution
                CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
              ENDIF

            ENDIF

            !copy the current solution to the previous iteration solution
            CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)

          CASE DEFAULT
            !do nothing
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE !the main time loop
        !do nothing
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_ConvergenceCheck")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ConvergenceCheck",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ConvergenceCheck")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ConvergenceCheck

  !
  !================================================================================================================================
  !

  !>Update the the bioelectric equation geometric field from the finite elasticity dependent field (deformed geometry)
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,CALC_CLOSEST_GAUSS_POINT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    LOGICAL, INTENT(IN) :: CALC_CLOSEST_GAUSS_POINT !<If true then the closest finite elasticity Gauss point for each bioelectrics node is calculated
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_ELASTICITY,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_ELASTICITY,GEOMETRIC_FIELD_MONODOMAIN,GEOMETRIC_FIELD_ELASTICITY
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_MONODOMAIN,INDEPENDENT_FIELD_MONODOMAIN,DEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_DEP_M,FIELD_VAR_GEO_M,FIELD_VAR_IND_FE,FIELD_VAR_IND_M,FIELD_VAR_IND_M_2
    INTEGER(INTG) :: component_idx,element_idx,ne,start_elem,START_ELEMENT,start_element_idx
    INTEGER(INTG) :: DEPENDENT_FIELD_INTERPOLATION,GEOMETRIC_FIELD_INTERPOLATION
    INTEGER(INTG) :: node_idx,node_idx_2,GAUSS_POINT,gauss_idx,fibre_idx
    INTEGER(INTG) :: nodes_in_Xi_1,nodes_in_Xi_2,nodes_in_Xi_3,n3,n2,n1,dof_idx,dof_idx2,idx,my_element_idx
    REAL(DP) :: VALUE,VALUE_LEFT,DISTANCE,VELOCITY,VELOCITY_MAX,OLD_DIST
    REAL(DP) :: XI(3),PREVIOUS_NODE(3),DIST_INIT,SARCO_LENGTH_INIT,TIME_STEP,DIST
    REAL(DP), POINTER :: GAUSS_POSITIONS(:,:)

    ENTERS("BioelectricFiniteElasticity_UpdateGeometricField",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(PROBLEM)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(INDEPENDENT_FIELD_ELASTICITY)
    NULLIFY(DEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(INDEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(DEPENDENT_FIELD_ELASTICITY)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(GEOMETRIC_FIELD_MONODOMAIN)
    NULLIFY(GEOMETRIC_FIELD_ELASTICITY)
    NULLIFY(ELEMENTS_TOPOLOGY)
    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(FIELD_VAR_DEP_M)
    NULLIFY(FIELD_VAR_GEO_M)
    NULLIFY(FIELD_VAR_IND_FE)
    NULLIFY(FIELD_VAR_IND_M)
    NULLIFY(FIELD_VAR_IND_M_2)
    NULLIFY(GAUSS_POSITIONS)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(3))

            CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)

              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
              !get the monodomain sub loop, solvers, solver, and finally geometric and field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    GEOMETRIC_FIELD_MONODOMAIN=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Geometric field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_MAPPING)
              NULLIFY(EQUATIONS_SET)
              NULLIFY(SOLVER_EQUATIONS)
              !get the finite elasticity sub loop, solvers, solver, and finally the dependent field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    DEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(DEPENDENT_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              DO component_idx=1,GEOMETRIC_FIELD_MONODOMAIN%VARIABLES(1)%NUMBER_OF_COMPONENTS
                !check for identical interpolation of the fields
                GEOMETRIC_FIELD_INTERPOLATION=GEOMETRIC_FIELD_MONODOMAIN%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
                DEPENDENT_FIELD_INTERPOLATION=DEPENDENT_FIELD_ELASTICITY%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
                IF(GEOMETRIC_FIELD_INTERPOLATION==DEPENDENT_FIELD_INTERPOLATION) THEN
                  !copy the dependent field components to the geometric field
                  CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & component_idx,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="The interpolation type of component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR, &
                    & ERROR))//" of field number "//TRIM(NUMBER_TO_VSTRING(GEOMETRIC_FIELD_MONODOMAIN%USER_NUMBER,"*",ERR, &
                    & ERROR))//" does not coincide with the interpolation type of field number " &
                    & //TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_ELASTICITY%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO

            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)

              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
              !get the monodomain sub loop, solvers, solver, and finally geometric field and dependent field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    GEOMETRIC_FIELD_MONODOMAIN=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Geometric field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    ! the Field_V_Variable_Type contains the 3D nodal positions
                    DEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(DEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_MAPPING)
              NULLIFY(EQUATIONS_SET)
              NULLIFY(SOLVER_EQUATIONS)
              !get the finite elasticity sub loop, solvers, solver, and finally the dependent and independent fields
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    GEOMETRIC_FIELD_ELASTICITY=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF


              node_idx=0
              node_idx_2=0
              fibre_idx=0
              CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VAR_DEP_M,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VAR_GEO_M,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VAR_IND_FE,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_IND_M,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE,FIELD_VAR_IND_M_2,ERR,ERROR,*999)

              NODES_MAPPING=>GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION% &
                & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES
              
              ELEMENTS_TOPOLOGY=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%TOPOLOGY%ELEMENTS


              !get the maximum contraction velocity 
              dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(2)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,dof_idx,VELOCITY_MAX,ERR,ERROR,*999)
              !NOTE: VELOCITY_MAX is the max shortening velocity, and hence negative!!!
              !The max lengthening velocity is assumed to be   abs(VELOCITY_MAX)/2.0
              
              !get the time step of the elasticity problem
              TIME_STEP=CONTROL_LOOP_PARENT%TIME_LOOP%TIME_INCREMENT


              !loop over the elements of the finite elasticity mesh (internal and boundary elements)
              !no need to consider ghost elements here since only bioelectrical fields are changed
              DO element_idx=1,ELEMENTS_TOPOLOGY%NUMBER_OF_ELEMENTS
                ne=ELEMENTS_TOPOLOGY%ELEMENTS(element_idx)%LOCAL_NUMBER
                my_element_idx=element_idx

                !the Field_V_Variable_Type of the FE independent field contains the number of nodes in each Xi-direction of the bioelectrics grid
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_1,ERR,ERROR,*999)
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_2,ERR,ERROR,*999)
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_3,ERR,ERROR,*999)
                !beginning of a fibre in this element: 1=yes, 0=no
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,start_elem,ERR,ERROR,*999)

                !if there is no bioelectrics grid in this finite elasticity element, or the fibres don't begin in this element, jump to the next element
                IF((nodes_in_Xi_1==0).OR.(nodes_in_Xi_2==0).OR.(nodes_in_Xi_3==0).OR.(start_elem==0)) CYCLE
                
                START_ELEMENT=ne
                start_element_idx=my_element_idx
                
                !assume Xi(1) to be normal to the seed surface, i.e. the seed points have Xi(1)=0
                XI=[0.0_DP,1.0_DP/(REAL(2*nodes_in_Xi_2)),1.0_DP/(REAL(2*nodes_in_Xi_3))]
                
                !assume that the bioelectrics node numbers are increased in order Xi(1), Xi(2), Xi(3) 
                DO n3=1,nodes_in_Xi_3
                  DO n2=1,nodes_in_Xi_2
                    fibre_idx=fibre_idx+1
                    
                    !loop over the FE elements that contain nodes of the very same fibres
                    DO
                      !get the finite elasticity dependent field interpolation parameters of this element
                      INTERPOLATION_PARAMETERS=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS &
                        & (FIELD_U_VARIABLE_TYPE)%PTR
                      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne,INTERPOLATION_PARAMETERS, &
                        & ERR,ERROR,*999)
                      INTERPOLATED_POINT=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR

                      !get the positions of the Gauss points of the Finite Elasticity element, GAUSS_POSITIONS(components,number_of_Gauss_points)
                      GAUSS_POSITIONS=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION% &
                        & MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP( &
                        & BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_POSITIONS
                      
                      DO n1=1,nodes_in_Xi_1
                        node_idx=node_idx+1
                        
                        !store the fibre number this bioelectrics node belongs to.
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,3,fibre_idx,ERR,ERROR,*999) 

                        !find the interpolated position of the bioelectric grid node from the FE dependent field
                        CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
                        !update the bioelectrics dependent field Field_V_Variable_Type
                        !the Field_V_Variable_Type of the monodomain dependent field contains the nodal positions in 3D
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,INTERPOLATED_POINT%VALUES(1,1),ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,2,INTERPOLATED_POINT%VALUES(2,1),ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,3,INTERPOLATED_POINT%VALUES(3,1),ERR,ERROR,*999)

                        IF((n1==1).AND.(ne==START_ELEMENT)) THEN
                          !a new line of bioelectrics grid nodes begins
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,0.0_DP,ERR,ERROR,*999)
                        ELSE
                          !get the position in 3D of the previous node
                          dof_idx2=FIELD_VAR_DEP_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(1),ERR,ERROR,*999)
                          dof_idx2=FIELD_VAR_DEP_M%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(2),ERR,ERROR,*999)
                          dof_idx2=FIELD_VAR_DEP_M%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(3),ERR,ERROR,*999)

                          !compute the distance between the previous node and the actual node
                          DIST=SQRT( &
                            & (INTERPOLATED_POINT%VALUES(1,1)-PREVIOUS_NODE(1))*(INTERPOLATED_POINT%VALUES(1,1)-PREVIOUS_NODE(1))+ &
                            & (INTERPOLATED_POINT%VALUES(2,1)-PREVIOUS_NODE(2))*(INTERPOLATED_POINT%VALUES(2,1)-PREVIOUS_NODE(2))+ &
                            & (INTERPOLATED_POINT%VALUES(3,1)-PREVIOUS_NODE(3))*(INTERPOLATED_POINT%VALUES(3,1)-PREVIOUS_NODE(3)))


                          !CONTRACTION VELOCITY CALCULATION
                          
                          !get the distance between the 2 nodes in the previous time step
                          dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx,OLD_DIST,ERR,ERROR,*999)
                          
                          !compute the new contraction velocity
                          VELOCITY=(DIST-OLD_DIST)/TIME_STEP
                          IF(.NOT. CALC_CLOSEST_GAUSS_POINT) THEN
                            !NOTE: VELOCITY_MAX is the max shortening velocity, and hence negative!!!
                            IF(VELOCITY<VELOCITY_MAX) THEN
                              CALL FLAG_WARNING('Exceeded maximum contraction velocity (shortening).',ERR,ERROR,*999)
                              VELOCITY=VELOCITY_MAX
                            !The max lengthening velocity is assumed to be VELOCITY_MAX/2.0
                            ELSEIF(VELOCITY>(ABS(VELOCITY_MAX)/2.0_DP)) THEN
                              CALL FLAG_WARNING('Exceeded maximum contraction velocity (lengthening).',ERR,ERROR,*999)
                              VELOCITY=-VELOCITY_MAX/2.0_DP
                            ENDIF
                          ENDIF
                          
                          !store the relative contraction velocity in component 3 of the U2 variable of the monodomain independent field
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,3,VELOCITY/ABS(VELOCITY_MAX),ERR,ERROR,*999)

                          !store the node distance for contraction velocity calculation
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,DIST,ERR,ERROR,*999)



                          !get the position in 1D of the previous node
                          dof_idx2=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,VALUE_LEFT,ERR,ERROR,*999)
                          !update the current 1D node position
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,VALUE_LEFT+DIST,ERR,ERROR,*999)

                          !get the initial sarcomere half length and initial node distance
                          dof_idx=FIELD_VAR_IND_M%COMPONENTS(2)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx,SARCO_LENGTH_INIT,ERR,ERROR,*999)
                          dof_idx=FIELD_VAR_IND_M%COMPONENTS(3)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx,DIST_INIT,ERR,ERROR,*999)
                          !update the current sarcomere half length
                          VALUE=SARCO_LENGTH_INIT*DIST/DIST_INIT
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,VALUE,ERR,ERROR,*999)

                          !update the first node to the same value as the second node (no better info available)
                          IF((n1==2).AND.(ne==START_ELEMENT)) THEN
                            !current sarcomere half length
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,1,1,node_idx-1,1,VALUE,ERR,ERROR,*999)                            
                            !old node distance
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,1,1,node_idx-1,1,DIST,ERR,ERROR,*999)
                            !relative contraction velocity
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,1,1,node_idx-1,3,VELOCITY/ABS(VELOCITY_MAX),ERR,ERROR,*999)
                          ENDIF

                        ENDIF !((n1==1).AND.(ne==START_ELEMENT))
                          
                        IF(CALC_CLOSEST_GAUSS_POINT) THEN
                          !calculate the closest finite elasticity Gauss point of each bioelectrics node
                          DISTANCE=1000000.0_DP
                          GAUSS_POINT=0
                          DO gauss_idx=1,SIZE(GAUSS_POSITIONS,2)
                            !compute the distance between the bioelectrics node and the Gauss point
                            VALUE=SQRT( &
                              & (Xi(1)-GAUSS_POSITIONS(1,gauss_idx))*(Xi(1)-GAUSS_POSITIONS(1,gauss_idx))+ &
                              & (Xi(2)-GAUSS_POSITIONS(2,gauss_idx))*(Xi(2)-GAUSS_POSITIONS(2,gauss_idx))+ &
                              & (Xi(3)-GAUSS_POSITIONS(3,gauss_idx))*(Xi(3)-GAUSS_POSITIONS(3,gauss_idx)))
                            IF(VALUE<DISTANCE) THEN
                              DISTANCE=VALUE
                              GAUSS_POINT=gauss_idx
                            ENDIF
                          ENDDO !gauss_idx
                          IF(GAUSS_POINT==0) CALL FLAG_WARNING("Closest Gauss Point not found",ERR,ERROR,*999)
                          !store the nearest Gauss Point info and the inElement info (local element number!!!)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,4,GAUSS_POINT,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,5,ne,ERR,ERROR,*999)
                        ENDIF !CALC_CLOSEST_GAUSS_POINT
                        
                        IF(start_elem==1) THEN
                          !fibres start in this element
                          XI(1)=XI(1)+1.0_DP/(REAL(nodes_in_Xi_1-1))
                        ELSEIF(start_elem==0) THEN
                          !fibres don't start in this element
                          XI(1)=XI(1)+1.0_DP/(REAL(nodes_in_Xi_1))
                        ELSE
                          LOCAL_ERROR="The start element index is incorrect. The index is "// &
                            & TRIM(NUMBER_TO_VSTRING(start_elem,"*",ERR,ERROR))// &
                            & " and should be zero or one." 
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                        
                      ENDDO !n1
                      



!tomo new
                      !smooth of the velocity field

                      !arithmetic mean of all rel_velo values within one FE element
!                      VELOCITY=0.0_DP
!                      DO n1=1,nodes_in_Xi_1
!                        node_idx_2=node_idx_2+1
!                        dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx_2)% &
!                          & DERIVATIVES(1)%VERSIONS(1)
!                        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                          & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                        VELOCITY=VELOCITY+VALUE
!                      ENDDO
!                      VELOCITY=VELOCITY/nodes_in_Xi_1

!                      node_idx_2=node_idx_2-nodes_in_Xi_1
!                      DO n1=1,nodes_in_Xi_1
!                        node_idx_2=node_idx_2+1
!                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                          & FIELD_VALUES_SET_TYPE,1,1,node_idx_2,3,VELOCITY,ERR,ERROR,*999)
!                      ENDDO

!--------------------------------------------------------------------------

!                      !moving average
!                      offset=3
!                      
!                      !do the first three nodes of a fibre manually - arithmetic mean
!                      VELOCITY=0.0_DP

!                      node_idx_2=node_idx_2+1
!                      
!                      dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx_2)% &
!                        & DERIVATIVES(1)%VERSIONS(1)
!                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                      VELOCITY=VELOCITY+VALUE

!                      node_idx_2=node_idx_2+1
!                      dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx_2)% &
!                        & DERIVATIVES(1)%VERSIONS(1)
!                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                      VELOCITY=VELOCITY+VALUE

!                      node_idx_2=node_idx_2+1
!                      dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx_2)% &
!                        & DERIVATIVES(1)%VERSIONS(1)
!                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                      VELOCITY=VELOCITY+VALUE
!                      
!                      VELOCITY=VELOCITY/offset

!                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,1,1,node_idx_2-2,3,VELOCITY,ERR,ERROR,*999)
!                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,1,1,node_idx_2-1,3,VELOCITY,ERR,ERROR,*999)
!                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,1,1,node_idx_2,3,VELOCITY,ERR,ERROR,*999)

!                      !do the major part as moving average
!                      DO n1=1+offset,nodes_in_Xi_1-offset
!                        node_idx_2=node_idx_2+1
!                        VELOCITY=0.0_DP
!                        DO n4=node_idx_2-offset,node_idx_2+offset
!                          dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(n4)% &
!                            & DERIVATIVES(1)%VERSIONS(1)
!                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                          VELOCITY=VELOCITY+VALUE
!                        ENDDO !n4
!                        VELOCITY=VELOCITY/(2*offset+1)
!                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                          & FIELD_VALUES_SET_TYPE,1,1,node_idx_2,3,VELOCITY,ERR,ERROR,*999)
!                      ENDDO !n1
!                      
!                      !do the last three nodes of a fibre manually - arithmetic mean
!                      VELOCITY=0.0_DP

!                      node_idx_2=node_idx_2+1
!                      dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx_2)% &
!                        & DERIVATIVES(1)%VERSIONS(1)
!                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                      VELOCITY=VELOCITY+VALUE

!                      node_idx_2=node_idx_2+1
!                      dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx_2)% &
!                        & DERIVATIVES(1)%VERSIONS(1)
!                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                      VELOCITY=VELOCITY+VALUE

!                      node_idx_2=node_idx_2+1
!                      dof_idx=FIELD_VAR_IND_M_2%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx_2)% &
!                        & DERIVATIVES(1)%VERSIONS(1)
!                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
!                      VELOCITY=VELOCITY+VALUE
!                      
!                      VELOCITY=VELOCITY/offset

!                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,1,1,node_idx_2-2,3,VELOCITY,ERR,ERROR,*999)
!                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,1,1,node_idx_2-1,3,VELOCITY,ERR,ERROR,*999)
!                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
!                        & FIELD_VALUES_SET_TYPE,1,1,node_idx_2,3,VELOCITY,ERR,ERROR,*999)
!!tomo new end                       


                      !if there is not an adjacent element in positive XI_1 direction, go to the next FE element
                      IF(ELEMENTS_TOPOLOGY%ELEMENTS(my_element_idx)%ADJACENT_ELEMENTS(1)%NUMBER_OF_ADJACENT_ELEMENTS==0) EXIT
                      
                      !consider the adjacent element in positive XI_1 direction
                      ne=ELEMENTS_TOPOLOGY%ELEMENTS(my_element_idx)%ADJACENT_ELEMENTS(1)%ADJACENT_ELEMENTS(1)

                      !if a fibre starts in the next element, go to the next FE elem
                      dof_idx=FIELD_VAR_IND_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,dof_idx,start_elem,ERR,ERROR,*999)
                      !beginning of a fibre in this element: 1=yes, 0=no
                      IF (start_elem==1) THEN
                        ne=ELEMENTS_TOPOLOGY%ELEMENTS(element_idx)%LOCAL_NUMBER
                        EXIT
                      ENDIF
                      
                      !find the element_idx that corresponds to ne
                      my_element_idx=0
                      DO idx=1,ELEMENTS_TOPOLOGY%NUMBER_OF_ELEMENTS
                        IF(ne==ELEMENTS_TOPOLOGY%ELEMENTS(idx)%ADJACENT_ELEMENTS(0)%ADJACENT_ELEMENTS(1)) THEN
                          my_element_idx=idx
                          EXIT
                        ENDIF
                      ENDDO
                      IF(my_element_idx==0) CALL FlagError("my_element_idx not found.",ERR,ERROR,*999)                      

                      dof_idx=FIELD_VAR_IND_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,dof_idx,nodes_in_Xi_1,ERR,ERROR,*999)
                      
                      start_elem=0 !fibres don't start in this element
                      
                      XI(1)=1.0_DP/(REAL(nodes_in_Xi_1))

                    ENDDO !
                    !for the beginning of the next fibre, go back to the element in which the last fibre started
                    ne=START_ELEMENT

                    dof_idx=FIELD_VAR_IND_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,nodes_in_Xi_1,ERR,ERROR,*999)

                    my_element_idx=start_element_idx
                    start_elem=1 !fibres start in this element
                    XI(1)=0.0_DP
                    XI(2)=XI(2)+1.0_DP/(REAL(nodes_in_Xi_2))
                  ENDDO !n2
                  XI(1)=0.0_DP
                  XI(2)=1.0_DP/(REAL(2*nodes_in_Xi_2))
                  XI(3)=Xi(3)+1.0_DP/(REAL(nodes_in_Xi_3))
                ENDDO !n3

              ENDDO !element_idx

            CASE(PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)

              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
              !get the monodomain sub loop, solvers, solver, and finally geometric field and dependent field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    GEOMETRIC_FIELD_MONODOMAIN=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Geometric field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    ! the Field_V_Variable_Type contains the 3D nodal positions
                    DEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(DEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_MAPPING)
              NULLIFY(EQUATIONS_SET)
              NULLIFY(SOLVER_EQUATIONS)
              !get the finite elasticity sub loop, solvers, solver, and finally the dependent and independent fields
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    GEOMETRIC_FIELD_ELASTICITY=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF


              node_idx=0
              node_idx_2=0
              fibre_idx=0
              CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VAR_DEP_M,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VAR_GEO_M,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VAR_IND_FE,ERR,ERROR,*999)

              NODES_MAPPING=>GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION% &
                & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES
              
              ELEMENTS_TOPOLOGY=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%TOPOLOGY%ELEMENTS

              !loop over the elements of the finite elasticity mesh (internal and boundary elements)
              !no need to consider ghost elements here since only bioelectrical fields are changed
              DO element_idx=1,ELEMENTS_TOPOLOGY%NUMBER_OF_ELEMENTS
                ne=ELEMENTS_TOPOLOGY%ELEMENTS(element_idx)%LOCAL_NUMBER
                my_element_idx=element_idx

                !the Field_V_Variable_Type of the FE independent field contains the number of nodes in each Xi-direction of the bioelectrics grid
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_1,ERR,ERROR,*999)
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_2,ERR,ERROR,*999)
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_3,ERR,ERROR,*999)
                !beginning of a fibre in this element: 1=yes, 0=no
                dof_idx=FIELD_VAR_IND_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,start_elem,ERR,ERROR,*999)

                !if there is no bioelectrics grid in this finite elasticity element, or the fibres don't begin in this element, jump to the next element
                IF((nodes_in_Xi_1==0).OR.(nodes_in_Xi_2==0).OR.(nodes_in_Xi_3==0).OR.(start_elem==0)) CYCLE
                
                START_ELEMENT=ne
                start_element_idx=my_element_idx
                
                !assume Xi(1) to be normal to the seed surface, i.e. the seed points have Xi(1)=0
                XI=[0.0_DP,1.0_DP/(REAL(2*nodes_in_Xi_2)),1.0_DP/(REAL(2*nodes_in_Xi_3))]
                
                !assume that the bioelectrics node numbers are increased in order Xi(1), Xi(2), Xi(3) 
                DO n3=1,nodes_in_Xi_3
                  DO n2=1,nodes_in_Xi_2
                    fibre_idx=fibre_idx+1
                    
                    !loop over the FE elements that contain nodes of the very same fibres
                    DO
                      !get the finite elasticity dependent field interpolation parameters of this element
                      INTERPOLATION_PARAMETERS=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS &
                        & (FIELD_U_VARIABLE_TYPE)%PTR
                      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne,INTERPOLATION_PARAMETERS, &
                        & ERR,ERROR,*999)
                      INTERPOLATED_POINT=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR

                      !get the positions of the Gauss points of the Finite Elasticity element, GAUSS_POSITIONS(components,number_of_Gauss_points)
                      GAUSS_POSITIONS=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION% &
                        & MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP( &
                        & BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_POSITIONS
                      
                      DO n1=1,nodes_in_Xi_1
                        node_idx=node_idx+1
                        
                        !store the fibre number this bioelectrics node belongs to.
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,3,fibre_idx,ERR,ERROR,*999) 

                        !find the interpolated position of the bioelectric grid node from the FE dependent field
                        CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
                        !update the bioelectrics dependent field Field_V_Variable_Type
                        !the Field_V_Variable_Type of the monodomain dependent field contains the nodal positions in 3D
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,INTERPOLATED_POINT%VALUES(1,1),ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,2,INTERPOLATED_POINT%VALUES(2,1),ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,3,INTERPOLATED_POINT%VALUES(3,1),ERR,ERROR,*999)

                        IF((n1==1).AND.(ne==START_ELEMENT)) THEN
                          !a new line of bioelectrics grid nodes begins
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,0.0_DP,ERR,ERROR,*999)
                        ELSE
                          !get the position in 3D of the previous node
                          dof_idx2=FIELD_VAR_DEP_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(1),ERR,ERROR,*999)
                          dof_idx2=FIELD_VAR_DEP_M%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(2),ERR,ERROR,*999)
                          dof_idx2=FIELD_VAR_DEP_M%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(3),ERR,ERROR,*999)

                          !compute the distance between the previous node and the actual node
                          DIST=SQRT( &
                            & (INTERPOLATED_POINT%VALUES(1,1)-PREVIOUS_NODE(1))*(INTERPOLATED_POINT%VALUES(1,1)-PREVIOUS_NODE(1))+ &
                            & (INTERPOLATED_POINT%VALUES(2,1)-PREVIOUS_NODE(2))*(INTERPOLATED_POINT%VALUES(2,1)-PREVIOUS_NODE(2))+ &
                            & (INTERPOLATED_POINT%VALUES(3,1)-PREVIOUS_NODE(3))*(INTERPOLATED_POINT%VALUES(3,1)-PREVIOUS_NODE(3)))

                          !get the position in 1D of the previous node
                          dof_idx2=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,dof_idx2,VALUE_LEFT,ERR,ERROR,*999)
                          !update the current 1D node position
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,VALUE_LEFT+DIST,ERR,ERROR,*999)

                        ENDIF !((n1==1).AND.(ne==START_ELEMENT))
                          
                        IF(CALC_CLOSEST_GAUSS_POINT) THEN
                          !calculate the closest finite elasticity Gauss point of each bioelectrics node
                          DISTANCE=1000000.0_DP
                          GAUSS_POINT=0
                          DO gauss_idx=1,SIZE(GAUSS_POSITIONS,2)
                            !compute the distance between the bioelectrics node and the Gauss point
                            VALUE=SQRT( &
                              & (Xi(1)-GAUSS_POSITIONS(1,gauss_idx))*(Xi(1)-GAUSS_POSITIONS(1,gauss_idx))+ &
                              & (Xi(2)-GAUSS_POSITIONS(2,gauss_idx))*(Xi(2)-GAUSS_POSITIONS(2,gauss_idx))+ &
                              & (Xi(3)-GAUSS_POSITIONS(3,gauss_idx))*(Xi(3)-GAUSS_POSITIONS(3,gauss_idx)))
                            IF(VALUE<DISTANCE) THEN
                              DISTANCE=VALUE
                              GAUSS_POINT=gauss_idx
                            ENDIF
                          ENDDO !gauss_idx
                          IF(GAUSS_POINT==0) CALL FLAG_WARNING("Closest Gauss Point not found",ERR,ERROR,*999)
                          !store the nearest Gauss Point info and the inElement info (local element number!!!)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,4,GAUSS_POINT,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,5,ne,ERR,ERROR,*999)
                        ENDIF !CALC_CLOSEST_GAUSS_POINT
                        
                        IF(start_elem==1) THEN
                          !fibres start in this element
                          XI(1)=XI(1)+1.0_DP/(REAL(nodes_in_Xi_1-1))
                        ELSEIF(start_elem==0) THEN
                          !fibres don't start in this element
                          XI(1)=XI(1)+1.0_DP/(REAL(nodes_in_Xi_1))
                        ELSE
                          LOCAL_ERROR="The start element index is incorrect. The index is "// &
                            & TRIM(NUMBER_TO_VSTRING(start_elem,"*",ERR,ERROR))// &
                            & " and should be zero or one." 
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                        
                      ENDDO !n1


                      !if there is not an adjacent element in positive XI_1 direction, go to the next FE element
                      IF(ELEMENTS_TOPOLOGY%ELEMENTS(my_element_idx)%ADJACENT_ELEMENTS(1)%NUMBER_OF_ADJACENT_ELEMENTS==0) EXIT
                      
                      !consider the adjacent element in positive XI_1 direction
                      ne=ELEMENTS_TOPOLOGY%ELEMENTS(my_element_idx)%ADJACENT_ELEMENTS(1)%ADJACENT_ELEMENTS(1)

                      !if a fibre starts in the next element, go to the next FE elem
                      dof_idx=FIELD_VAR_IND_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,dof_idx,start_elem,ERR,ERROR,*999)
                      !beginning of a fibre in this element: 1=yes, 0=no
                      IF (start_elem==1) THEN
                        ne=ELEMENTS_TOPOLOGY%ELEMENTS(element_idx)%LOCAL_NUMBER
                        EXIT
                      ENDIF
                      
                      !find the element_idx that corresponds to ne
                      my_element_idx=0
                      DO idx=1,ELEMENTS_TOPOLOGY%NUMBER_OF_ELEMENTS
                        IF(ne==ELEMENTS_TOPOLOGY%ELEMENTS(idx)%ADJACENT_ELEMENTS(0)%ADJACENT_ELEMENTS(1)) THEN
                          my_element_idx=idx
                          EXIT
                        ENDIF
                      ENDDO
                      IF(my_element_idx==0) CALL FlagError("my_element_idx not found.",ERR,ERROR,*999)                      

                      dof_idx=FIELD_VAR_IND_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,dof_idx,nodes_in_Xi_1,ERR,ERROR,*999)
                      
                      start_elem=0 !fibres don't start in this element
                      
                      XI(1)=1.0_DP/(REAL(nodes_in_Xi_1))

                    ENDDO !
                    !for the beginning of the next fibre, go back to the element in which the last fibre started
                    ne=START_ELEMENT

                    dof_idx=FIELD_VAR_IND_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,nodes_in_Xi_1,ERR,ERROR,*999)

                    my_element_idx=start_element_idx
                    start_elem=1 !fibres start in this element
                    XI(1)=0.0_DP
                    XI(2)=XI(2)+1.0_DP/(REAL(nodes_in_Xi_2))
                  ENDDO !n2
                  XI(1)=0.0_DP
                  XI(2)=1.0_DP/(REAL(2*nodes_in_Xi_2))
                  XI(3)=Xi(3)+1.0_DP/(REAL(nodes_in_Xi_3))
                ENDDO !n3

              ENDDO !element_idx

            CASE DEFAULT
              LOCAL_ERROR="The third problem specification of "// &
                & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(2),"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity of a multi physics problem."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The second problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(2),"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("BioelectricFiniteElasticity_UpdateGeometricField")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_UpdateGeometricField",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_UpdateGeometricField")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_UpdateGeometricField

  !
  !================================================================================================================================
  !

  !>Interpolates the finite elasticity independent field from the biolectrics independent field.
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE BioelectricFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_ELASTICITY,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_MONODOMAIN,INDEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING,NODES_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE_U,FIELD_VARIABLE_V,FIELD_VARIABLE_FE
    INTEGER(INTG) :: node_idx,element_idx,gauss_idx,ne
    INTEGER(INTG) :: nearestGP,inElement,dof_idx
    INTEGER(INTG) :: NUMBER_OF_GAUSS_POINTS
    REAL(DP) :: ACTIVE_STRESS
    REAL(DP) :: TITIN_STRESS_UNBOUND,TITIN_STRESS_BOUND,TITIN_STRESS_CROSS_FIBRE_UNBOUND,TITIN_STRESS_CROSS_FIBRE_BOUND,ACTIVATION
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_GAUSS_POINTS=64
    INTEGER(INTG) :: NUMBER_OF_NODES(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: ACTIVE_STRESS_VALUES(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: TITIN_STRESS_VALUES_UNBOUND(MAX_NUMBER_OF_GAUSS_POINTS),TITIN_STRESS_VALUES_BOUND(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: ACTIVATION_VALUES(MAX_NUMBER_OF_GAUSS_POINTS)

    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(FIELD_VARIABLE_U)
    NULLIFY(FIELD_VARIABLE_V)
    NULLIFY(FIELD_VARIABLE_FE)
    
    ENTERS("BioelectricFiniteElasticity_IndependentFieldInterpolate",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        IF(.NOT.ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(problem%specification,1)<3) THEN
          CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
            & err,error,*999)
        END IF
        SELECT CASE(PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
          IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
            !--- MONODOMAIN ---
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) CALL FlagError("Independent field is not associated.", &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF

            !--- FINITE ELASTICITY ---
            NULLIFY(SOLVERS)
            NULLIFY(SOLVER)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) CALL FlagError("Independent field is not associated.",ERR, &
                    & ERROR,*999)
                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF

            !--- NOW INTERPOLATE ---
            ELEMENTS_MAPPING=>INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
            NODES_MAPPING=>INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES

            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_U,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE_V,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_FE,ERR,ERROR,*999)

            !loop over the finite elasticity elements
            !first process the internal and boundary elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%BOUNDARY_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              
              NUMBER_OF_GAUSS_POINTS=INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY% &
                & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP &
                & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS

              IF(NUMBER_OF_GAUSS_POINTS>MAX_NUMBER_OF_GAUSS_POINTS) CALL FlagError( & 
                & "NUMBER_OF_GAUSS_POINTS is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)
              NUMBER_OF_NODES=0
              ACTIVE_STRESS_VALUES=0.0_DP
              TITIN_STRESS_VALUES_UNBOUND=0.0_DP
              TITIN_STRESS_VALUES_BOUND=0.0_DP
              TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND=0.0_DP
              TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND=0.0_DP
              ACTIVATION_VALUES=0.0_DP
              
              !loop over the bioelectrics nodes
              DO node_idx=1,NODES_MAPPING%NUMBER_OF_LOCAL
                dof_idx=FIELD_VARIABLE_V%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,inElement,ERR,ERROR,*999) !component 5 of variable V contains inElem info (LOCAL NUMBERING!!!)

                !check if the bioelectrics node is located within the finite elasticity element
                IF(inElement==ne) THEN
                  !component 4 of variable V contains Nearest Gauss Point info
                  dof_idx=FIELD_VARIABLE_V%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                    & VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & dof_idx,nearestGP,ERR,ERROR,*999)
                  IF(nearestGP>MAX_NUMBER_OF_GAUSS_POINTS) CALL FlagError( &
                    & "Nearest Gauss Point is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)
                  !component 1 of variable U contains the active stress
                  dof_idx=FIELD_VARIABLE_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                    & VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)
                  
                  !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                  NUMBER_OF_NODES(nearestGP)=NUMBER_OF_NODES(nearestGP)+1
                  !add up the active stress value
                  ACTIVE_STRESS_VALUES(nearestGP)=ACTIVE_STRESS_VALUES(nearestGP)+ACTIVE_STRESS

                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                    !component 2 of variable U contains the titin stress unbound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_UNBOUND,ERR,ERROR,*999)
                    !component 3 of variable U contains the titin stress bound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_BOUND,ERR,ERROR,*999)
                    !component 4 of variable U contains the titin XF-stress (cross-fibre directions) unbound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_UNBOUND,ERR,ERROR,*999)
                    !component 5 of variable U contains the titin XF-stress (cross-fibre directions) bound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_BOUND,ERR,ERROR,*999)
                    !component 6 of variable U contains the titin activation
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(6)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,ACTIVATION,ERR,ERROR,*999)

                    TITIN_STRESS_VALUES_UNBOUND(nearestGP)=TITIN_STRESS_VALUES_UNBOUND(nearestGP)+TITIN_STRESS_UNBOUND
                    TITIN_STRESS_VALUES_BOUND(nearestGP)=TITIN_STRESS_VALUES_BOUND(nearestGP)+TITIN_STRESS_BOUND
                    TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(nearestGP)=TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(nearestGP) + &
                      & TITIN_STRESS_CROSS_FIBRE_UNBOUND
                    TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(nearestGP)=TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(nearestGP) + &
                      & TITIN_STRESS_CROSS_FIBRE_BOUND
                    ACTIVATION_VALUES(nearestGP)=ACTIVATION_VALUES(nearestGP)+ACTIVATION
                  ENDIF

                ENDIF
              ENDDO

              !loop over the finite elasticity Gauss points
              DO gauss_idx=1,NUMBER_OF_GAUSS_POINTS
                !make sure we don't divide by zero
                IF(NUMBER_OF_NODES(gauss_idx)<=0) THEN
                  ACTIVE_STRESS=0.0_DP
                  TITIN_STRESS_UNBOUND=0.0_DP
                  TITIN_STRESS_BOUND=0.0_DP
                  TITIN_STRESS_CROSS_FIBRE_UNBOUND=0.0_DP
                  TITIN_STRESS_CROSS_FIBRE_BOUND=0.0_DP
                  ACTIVATION=0.0_DP
                ELSE
                  ACTIVE_STRESS=ACTIVE_STRESS_VALUES(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_UNBOUND=TITIN_STRESS_VALUES_UNBOUND(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_BOUND=TITIN_STRESS_VALUES_BOUND(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_CROSS_FIBRE_UNBOUND=TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_CROSS_FIBRE_BOUND=TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(gauss_idx)/ &
                    & NUMBER_OF_NODES(gauss_idx)
                  ACTIVATION=ACTIVATION_VALUES(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                ENDIF

                dof_idx=FIELD_VARIABLE_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)

                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_UNBOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_BOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_UNBOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(5)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_BOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(6)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,ACTIVATION,ERR,ERROR,*999)
                ENDIF

              ENDDO !gauss_idx
            ENDDO !element_idx

            !now the ghost elements -- get the relevant info from the other computational nodes
            CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Independent field interpolation is not implemented for problem subtype " &
            & //TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_IndependentFieldInterpolate")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_IndependentFieldInterpolate",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_IndependentFieldInterpolate")
    RETURN 1

  END SUBROUTINE BioelectricFiniteElasticity_IndependentFieldInterpolate

  !
  !================================================================================================================================
  !

  !>Computes force enhancement based on the titin model of C Rode et al. (2009). Force depression is not yet implemented.
  !>Titin-induced force enhancement and force depression: A 'sticky-spring' mechanism in muscle contractions? C Rode, T Siebert, R Blickhan - Journal of theoretical biology, 2009.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_MONODOMAIN
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_IND_M
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    INTEGER(INTG) :: node_idx,dof_idx
    REAL(DP), PARAMETER :: LENGTH_ACTIN=1.04_DP,LENGTH_MBAND=0.0625_DP,LENGTH_MYOSIN=0.7375_DP
    REAL(DP), PARAMETER :: LENGTH_ZERO=0.635_DP,LENGTH_ZDISC=0.05_DP
    REAL(DP) :: F0,FORCE,TITIN_BOUND,TITIN_XF_BOUND,TITIN_UNBOUND,TITIN_XF_UNBOUND
    REAL(DP) :: ELONGATION,ELONGATION_NEW
    REAL(DP) :: ELONGATION_DIST_IG,ELONGATION_PEVK
    REAL(DP) :: LENGTH_TITIN,LENGTH_INIT_TITIN,LENGTH_DIST_IG_F0,LENGTH_DIST_IG
    REAL(DP) :: STIFFNESS_PEVK
    REAL(DP) :: SARCO_LENGTH_AT_ACTIVATION,SARCO_LENGTH
    REAL(DP) :: ACTIN_MYOSIN_DISTANCE,d10
    REAL(DP) :: SLOPE,STIFFNESS_DIST,FORCE_DISTAL_IG,DELTA_F,DIFF_QUOT
    REAL(DP), PARAMETER, DIMENSION(5) :: COEFF_MATRIX=[5.0239_DP,-0.6717_DP,-2.5841_DP,-5.0128_DP,-5.0239_DP]
    REAL(DP), DIMENSION(250) :: LENGTHS_DIST_IG,FORCES_DIST_IG
    REAL(DP), PARAMETER :: DX=0.001_DP
    REAL(DP), PARAMETER :: FORCE_INCREMENT=1.e-5_DP,TOL=1.e-5_DP
    INTEGER(INTG) :: INDEX_REF,INDEX_PSEUDO,INDEX_I
    INTEGER(INTG), PARAMETER :: DIM_DATA=250
    INTEGER(INTG) :: SWITCH_MODEL
    
    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(PROBLEM)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(INDEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(FIELD_VAR_IND_M)
    
    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN",ERR,ERROR,*999)

    !the realtion between length and force of the distal_Ig region is very nonlinear. Linear interpolation of Rode's data is used here.
    FORCES_DIST_IG= &
      & [0.0_DP,0.001_DP,0.002_DP,0.003_DP,0.004_DP,0.005_DP,0.006_DP,0.007_DP,0.008_DP,0.009_DP,0.01_DP,0.011_DP,0.012_DP, &
      & 0.013_DP,0.014_DP,0.015_DP,0.016_DP,0.017_DP,0.018_DP,0.019_DP,0.02_DP,0.021_DP,0.022_DP,0.023_DP,0.024_DP, &
      & 0.025_DP,0.026_DP,0.027_DP,0.028_DP,0.029_DP,0.03_DP,0.031_DP,0.032_DP,0.033_DP,0.034_DP,0.035_DP,0.036_DP, &
      & 0.037_DP,0.038_DP,0.039_DP,0.04_DP,0.041_DP,0.042_DP,0.043_DP,0.044_DP,0.045_DP,0.046_DP,0.047_DP,0.048_DP, &
      & 0.049_DP,0.05_DP,0.051_DP,0.052_DP,0.053_DP,0.054_DP,0.055_DP,0.056_DP,0.057_DP,0.058_DP,0.059_DP,0.06_DP, &
      & 0.061_DP,0.062_DP,0.063_DP,0.064_DP,0.065_DP,0.066_DP,0.067_DP,0.068_DP,0.069_DP,0.070_DP,0.071_DP,0.072_DP, &
      & 0.073_DP,0.074_DP,0.075_DP,0.076_DP,0.077_DP,0.078_DP,0.079_DP,0.080_DP,0.081_DP,0.082_DP,0.083_DP,0.084_DP, &
      & 0.085_DP,0.086_DP,0.087_DP,0.088_DP,0.089_DP,0.090_DP,0.091_DP,0.092_DP,0.093_DP,0.094_DP,0.095_DP,0.096_DP, &
      & 0.097_DP,0.098_DP,0.099_DP,0.1_DP,0.101_DP,0.102_DP,0.103_DP,0.104_DP,0.105_DP,0.106_DP,0.107_DP,0.108_DP, &
      & 0.109_DP,0.11_DP,0.111_DP,0.112_DP,0.113_DP,0.114_DP,0.115_DP,0.116_DP,0.117_DP,0.118_DP,0.119_DP,0.12_DP, &
      & 0.121_DP,0.122_DP,0.123_DP,0.124_DP,0.125_DP,0.126_DP,0.127_DP,0.128_DP,0.129_DP,0.13_DP,0.131_DP,0.132_DP, &
      & 0.133_DP,0.134_DP,0.135_DP,0.136_DP,0.137_DP,0.138_DP,0.139_DP,0.14_DP,0.141_DP,0.142_DP,0.143_DP,0.144_DP, &
      & 0.145_DP,0.146_DP,0.147_DP,0.148_DP,0.149_DP,0.15_DP,0.151_DP,0.152_DP,0.153_DP,0.154_DP,0.155_DP,0.156_DP, &
      & 0.157_DP,0.158_DP,0.159_DP,0.16_DP,0.161_DP,0.162_DP,0.163_DP,0.164_DP,0.165_DP,0.166_DP,0.167_DP, 0.168_DP, &
      & 0.169_DP,0.17_DP,0.171_DP,0.172_DP,0.173_DP,0.174_DP,0.175_DP,0.176_DP,0.177_DP,0.178_DP,0.179_DP,0.18_DP, &
      & 0.181_DP,0.182_DP,0.183_DP,0.184_DP,0.185_DP,0.186_DP,0.187_DP,0.188_DP,0.189_DP,0.19_DP,0.191_DP,0.192_DP, &
      & 0.193_DP,0.194_DP,0.195_DP,0.196_DP,0.197_DP,0.198_DP,0.199_DP,0.2_DP,0.201_DP,0.202_DP,0.203_DP,0.204_DP, &
      & 0.205_DP,0.206_DP,0.207_DP,0.208_DP,0.209_DP,0.21_DP,0.211_DP,0.212_DP,0.213_DP,0.214_DP,0.215_DP,0.216_DP, &
      & 0.217_DP,0.218_DP,0.219_DP,0.22_DP,0.221_DP,0.222_DP,0.223_DP,0.224_DP,0.225_DP,0.226_DP,0.227_DP,0.228_DP, &
      & 0.229_DP,0.23_DP,0.231_DP,0.232_DP,0.233_DP,0.234_DP,0.235_DP,0.236_DP,0.237_DP,0.238_DP,0.239_DP,0.24_DP, &
      & 0.241_DP,0.242_DP,0.243_DP,0.244_DP,0.245_DP,0.246_DP,0.247_DP,0.248_DP,0.249_DP]
    LENGTHS_DIST_IG= &
      & [0.0_DP,0.03461753545561_DP,0.049729169766010_DP,0.058506390323323_DP,0.064606296848594_DP,0.06922519775133_DP, &
      & 0.0729080120998386_DP,0.0759458896446241_DP,0.0785230355395668_DP,0.0807314335143191_DP,0.0826674161660979_DP, &
      & 0.0843819721302_DP,0.0859161360822_DP,0.087300738288_DP,0.0885510536196_DP,0.08970061165_DP,0.090751366_DP, &
      & 0.0917285714001_DP,0.0926271710799_DP,0.093467026018_DP,0.094254010845_DP,0.094992014919_DP,0.09568255451_DP, &
      & 0.0963346932312_DP,0.0969518155718_DP,0.097537000419_DP,0.098093047899_DP,0.098622503780_DP,0.09912768169_DP, &
      & 0.0996103583026_DP,0.1000734170008_DP,0.100517614158_DP,0.100942907967_DP,0.101351270601_DP,0.101745244913_DP, &
      & 0.1021260518375_DP,0.1024947934706_DP,0.102851281365_DP,0.103194317381_DP,0.103528086731_DP,0.103853341524_DP, &
      & 0.1041686065999_DP,0.1044739635997_DP,0.104772842609_DP,0.105064758806_DP,0.105347078318_DP,0.105624524436_DP, &
      & 0.1058959740402_DP,0.1061596374590_DP,0.106419674124_DP,0.106673222977_DP,0.106921779223_DP,0.107167251003_DP, &
      & 0.1074055929211_DP,0.1076418938850_DP,0.107872798106_DP,0.108100580266_DP,0.108324935479_DP,0.108545154704_DP, &
      & 0.1087634145225_DP,0.1089769308376_DP,0.109189523632_DP,0.109397109133_DP,0.109604335645_DP,0.109806785604_DP, &
      & 0.1100090293896_DP,0.1102069598668_DP,0.110404790711_DP,0.110598542577_DP,0.110792294444_DP,0.110982362315_DP, &
      & 0.1111722575399_DP,0.1113591719379_DP,0.111545514634_DP,0.111729654458_DP,0.111912731875_DP,0.112094428489_DP, &
      & 0.1122745119339_DP,0.1124540532647_DP,0.112631398979_DP,0.112808744694_DP,0.112983883284_DP,0.113158733268_DP, &
      & 0.1133324054841_DP,0.1135049882670_DP,0.113677360515_DP,0.113847891889_DP,0.114018423263_DP,0.114187784974_DP, &
      & 0.1143564686814_DP,0.1145249703398_DP,0.114691998719_DP,0.114859027099_DP,0.115025270473_DP,0.115190825075_DP, &
      & 0.1153563796784_DP,0.1155207617063_DP,0.115685013868_DP,0.115849024045_DP,0.116012135436_DP,0.116175246828_DP, &
      & 0.1163378963393_DP,0.1165000194746_DP,0.116662142610_DP,0.116823704015_DP,0.116984982740_DP,0.117146261465_DP, &
      & 0.1173069696799_DP,0.1174675396300_DP,0.117628109580_DP,0.117788165045_DP,0.117948154078_DP,0.118108143111_DP, &
      & 0.1182677147299_DP,0.1184272433357_DP,0.118586771941_DP,0.118745999791_DP,0.118905181478_DP,0.119064363164_DP, &
      & 0.1192233610196_DP,0.1193823026790_DP,0.119541244338_DP,0.119700101992_DP,0.119858904246_DP,0.120017706500_DP, &
      & 0.1201764919257_DP,0.1203352494533_DP,0.120494006981_DP,0.120652768322_DP,0.120811570168_DP,0.120970372014_DP, &
      & 0.1211291738603_DP,0.1212880693042_DP,0.121446999172_DP,0.121605929041_DP,0.121764923098_DP,0.121924059629_DP, &
      & 0.1220831961613_DP,0.1222423326927_DP,0.122401700265_DP,0.122561117298_DP,0.122720534331_DP,0.122880045624_DP, &
      & 0.1230398124447_DP,0.1231995792652_DP,0.123359346085_DP,0.123519381317_DP,0.123679562894_DP,0.123839744470_DP, &
      & 0.1239999260467_DP,0.1241605622478_DP,0.124321219453_DP,0.124481876659_DP,0.124642637543_DP,0.124803827368_DP, &
      & 0.1249650171945_DP,0.1251262070202_DP,0.125287609380_DP,0.125449385133_DP,0.125611160886_DP,0.125772936638_DP, &
      & 0.1259350043157_DP,0.1260974158090_DP,0.126259827302_DP,0.126422238795_DP,0.126584979901_DP,0.126748073636_DP, &
      & 0.1269111673704_DP,0.1270742611047_DP,0.127237669513_DP,0.127401488846_DP,0.127565308178_DP,0.127729127511_DP, &
      & 0.1278931842348_DP,0.1280577695422_DP,0.128222354849_DP,0.128386940157_DP,0.128551614623_DP,0.128717003455_DP, &
      & 0.1288823922862_DP,0.1290477811173_DP,0.129213169948_DP,0.129379259581_DP,0.129545486803_DP,0.129711714025_DP, &
      & 0.1298779412472_DP,0.1300445897140_DP,0.130211687650_DP,0.130378785586_DP,0.130545883522_DP,0.130713029866_DP, &
      & 0.1308810284274_DP,0.1310490269884_DP,0.131217025549_DP,0.131385024110_DP,0.131553528414_DP,0.131722455223_DP, &
      & 0.1318913820322_DP,0.1320603088411_DP,0.132229235650_DP,0.132399074312_DP,0.132568954822_DP,0.132738835332_DP, &
      & 0.1329087158420_DP,0.1330788764544_DP,0.133249734060_DP,0.133420591667_DP,0.133591449273_DP,0.133762306879_DP, &
      & 0.1339336992411_DP,0.1341055553883_DP,0.134277411535_DP,0.134449267682_DP,0.134621123829_DP,0.134793694483_DP, &
      & 0.1349665687658_DP,0.1351394430487_DP,0.135312317331_DP,0.135485191614_DP,0.135658878561_DP,0.135832788820_DP, &
      & 0.1360066990800_DP,0.1361806093395_DP,0.136354519599_DP,0.136529253243_DP,0.136704215658_DP,0.136879178072_DP, &
      & 0.1370541404878_DP,0.1372291029027_DP,0.137404806924_DP,0.137580836098_DP,0.137756865271_DP,0.137932894445_DP, &
      & 0.1381089236188_DP,0.1382855157880_DP,0.138462624830_DP,0.138639733872_DP,0.138816842915_DP,0.138993951957_DP, &
      & 0.1391713448843_DP,0.1393495454909_DP,0.139527746097_DP,0.139705946704_DP,0.139884147310_DP,0.140062347917_DP, &
      & 0.1402415516727_DP,0.1404208541990_DP,0.140600156725270_DP,0.140779459251523_DP,0.140958761777775_DP]

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        SELECT CASE(PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
          IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
            !The first control_loop is the one for monodomain
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
            !The second solver is associated with the diffusion part of the monodomain equation
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) THEN
                    CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                  ENDIF

                  CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_IND_M,ERR,ERROR,*999)

                  NODES_MAPPING=>INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION% &
                    & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES

                  ! Initialization
                  INDEX_REF=1

                  DO node_idx=1,NODES_MAPPING%NUMBER_OF_LOCAL
                  
                    !the fourth component of the U1 variable contains the half sarcomere length at activation
                    dof_idx=FIELD_VAR_IND_M%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,SARCO_LENGTH_AT_ACTIVATION,ERR,ERROR,*999)

                    !the first component of the U1 variable contains the actual half sarcomere length
                    dof_idx=FIELD_VAR_IND_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,SARCO_LENGTH,ERR,ERROR,*999)
                    
                    ELONGATION=SARCO_LENGTH-SARCO_LENGTH_AT_ACTIVATION
                    
                    IF(ELONGATION.LT.0) THEN
                      LENGTH_TITIN=SARCO_LENGTH-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                      !function to approximate the relation between the initial titin length and the initial passive force F0	
                      TITIN_UNBOUND=COEFF_MATRIX(1)*EXP(LENGTH_TITIN)+COEFF_MATRIX(2)*LENGTH_TITIN**3+COEFF_MATRIX(3)* &
                        & LENGTH_TITIN**2+COEFF_MATRIX(4)*LENGTH_TITIN+COEFF_MATRIX(5) 
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,2,TITIN_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,3,TITIN_UNBOUND,ERR,ERROR,*999)

                      ! calculate x-fibre titin stress (trigonometry):                      
                      ! fitting function to calculate the filament-lattice parameter d10 in nanometer (from Elliot 1963 -mammalian muscle)
                      d10=-13.39_DP*SARCO_LENGTH+58.37_DP
                      ! calculate the distance between actin and myosin filament in micro meter (geometrical relationship)
                      ACTIN_MYOSIN_DISTANCE=0.001_DP*(2.0_DP/3.0_DP*d10)  
                      ! calculate x-fibre stress with tangens-function (unbound titin)
                      TITIN_XF_UNBOUND=0.5_DP*TITIN_UNBOUND*ACTIN_MYOSIN_DISTANCE/LENGTH_TITIN
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,4,TITIN_XF_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,5,TITIN_XF_UNBOUND,ERR,ERROR,*999)

                    ELSE  ! Force enhancement  

                      ! Calculate Titin-Force for unbound titin filaments
                      LENGTH_TITIN=SARCO_LENGTH-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                      !function to approximate the relation between the initial titin length and the initial passive force F0	
                      TITIN_UNBOUND=COEFF_MATRIX(1)*EXP(LENGTH_TITIN)+COEFF_MATRIX(2)*LENGTH_TITIN**3+COEFF_MATRIX(3)* &
                        & LENGTH_TITIN**2+COEFF_MATRIX(4)*LENGTH_TITIN+COEFF_MATRIX(5)

                      ! Calculate Titin-Force for bound titin filaments
                      ! Switch between different force-enhancent models/implementations 
                      ! 1=linear approximation 2=square approximation 3=simple iterative solution of Rode's model
                      ! 4=Newton's-method to solve Rode's model
                      SWITCH_MODEL=4
                      !linear approximation
                      IF(SWITCH_MODEL.EQ.1) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
                        STIFFNESS_PEVK=1000.0_DP*(0.1880_DP*SARCO_LENGTH_AT_ACTIVATION**4-0.8694_DP*SARCO_LENGTH_AT_ACTIVATION**3+ &
                          & 1.5084_DP*SARCO_LENGTH_AT_ACTIVATION**2-1.1577_DP*SARCO_LENGTH_AT_ACTIVATION+0.3345_DP)
                        IF(ELONGATION.LT.0.02) THEN ! 0.02 -> offset value 
                          FORCE=0.0_DP
                        ELSE 
                          FORCE=(ELONGATION-0.02_DP)*STIFFNESS_PEVK
                        ENDIF
                        TITIN_BOUND=F0+FORCE

                      !square approximation
                      ELSEIF(SWITCH_MODEL.EQ.2) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        FORCE=2.0_DP*ELONGATION**2
                        TITIN_BOUND=F0+FORCE

                      !Rode's model -> Ekin's implementation
                      ELSEIF(SWITCH_MODEL.EQ.3) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
                        STIFFNESS_PEVK=1000.0_DP*(0.1880_DP*SARCO_LENGTH_AT_ACTIVATION**4-0.8694_DP*SARCO_LENGTH_AT_ACTIVATION**3+ &
                          & 1.5084_DP*SARCO_LENGTH_AT_ACTIVATION**2-1.1577_DP*SARCO_LENGTH_AT_ACTIVATION+0.3345_DP)
                      
                        IF(F0.LE.0) THEN
                          LENGTH_DIST_IG_F0=-35.63_DP+39.58889_DP*(F0+0.9_DP)
                        ELSEIF(F0.GE.0.24_DP) THEN
                          LENGTH_DIST_IG_F0=0.1411_DP+0.196576763485477_DP*(F0-0.2495_DP)
                        ELSE                    
                          INDEX_PSEUDO=CEILING(F0/DX)
                          INDEX_I=INDEX_REF+INDEX_PSEUDO-1
                          LENGTH_DIST_IG_F0=LENGTHS_DIST_IG(INDEX_I)-(LENGTHS_DIST_IG(INDEX_I+1)-LENGTHS_DIST_IG(INDEX_I))* &
                            & (FORCES_DIST_IG(INDEX_I)-F0)/(FORCES_DIST_IG(INDEX_I+1)-FORCES_DIST_IG(INDEX_I))
                        ENDIF          
	               
                        ELONGATION_NEW=-1.0_DP
                        FORCE=0.0_DP                  
                        DO WHILE (ELONGATION_NEW.LT.ELONGATION)
                          FORCE=FORCE+FORCE_INCREMENT
                          TITIN_BOUND=FORCE+F0
                          ELONGATION_PEVK=FORCE/STIFFNESS_PEVK  
                          IF (TITIN_BOUND.LE.0.0_DP) THEN
                            LENGTH_DIST_IG=-35.63_DP+39.58889_DP*(TITIN_BOUND+0.9_DP)
                          ELSE IF (TITIN_BOUND.GE.0.24_DP) THEN
                            LENGTH_DIST_IG=0.1411_DP+0.196576763485477_DP*(TITIN_BOUND-0.2495_DP)
                          ELSE                  
                            INDEX_PSEUDO=CEILING(TITIN_BOUND/DX)
                            INDEX_I=INDEX_REF+INDEX_PSEUDO-1
                            LENGTH_DIST_IG=LENGTHS_DIST_IG(INDEX_I)-(LENGTHS_DIST_IG(INDEX_I+1)-LENGTHS_DIST_IG( &
                              & INDEX_I))*(FORCES_DIST_IG(INDEX_I)-TITIN_BOUND)/(FORCES_DIST_IG(INDEX_I+1)-FORCES_DIST_IG(INDEX_I))
                          END IF
                          ELONGATION_DIST_IG=LENGTH_DIST_IG-LENGTH_DIST_IG_F0                  
                          ELONGATION_NEW=ELONGATION_PEVK+ELONGATION_DIST_IG
                        ENDDO

                      !Rode's titin model -> solve with Newton's method
                      ELSEIF(SWITCH_MODEL.EQ.4) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
                        STIFFNESS_PEVK=1000.0_DP*(0.1880_DP*SARCO_LENGTH_AT_ACTIVATION**4-0.8694_DP*SARCO_LENGTH_AT_ACTIVATION**3+ &
                          & 1.5084_DP*SARCO_LENGTH_AT_ACTIVATION**2-1.1577_DP*SARCO_LENGTH_AT_ACTIVATION+0.3345_DP)

                        !calculate LENGTH_DIST_IG_F0 with linear inter- or extrapolation
                        INDEX_I=2
                        SLOPE=0.0_DP
                        LENGTH_DIST_IG_F0=0.0_DP
                        IF(F0.GE.FORCES_DIST_IG(DIM_DATA)) THEN
                          INDEX_I=DIM_DATA
                        ELSE
                          DO WHILE (F0.GT.FORCES_DIST_IG(INDEX_I))
                            INDEX_I=INDEX_I+1
                          ENDDO
                        ENDIF
                        SLOPE=(FORCES_DIST_IG(INDEX_I)-FORCES_DIST_IG(INDEX_I-1))/ &
                          & (LENGTHS_DIST_IG(INDEX_I)-LENGTHS_DIST_IG(INDEX_I-1))
                        LENGTH_DIST_IG_F0=LENGTHS_DIST_IG(INDEX_I-1)+SLOPE*(F0-FORCES_DIST_IG(INDEX_I-1))

                        !initialize Newton-method to calculate the titin force
                        INDEX_I=2                
                        STIFFNESS_DIST=1.0_DP   
                        FORCE_DISTAL_IG=F0      
                        FORCE=0.0_DP                  ! delta P
                        TITIN_BOUND=0.0_DP            ! total titin force (= P_PEVK)
                        DELTA_F=10.0_DP               ! start value for iteration
                        DIFF_QUOT=1.0_DP              ! numerical derivative              
                        ELONGATION_PEVK=ELONGATION/2.0_DP                               
                        LENGTH_DIST_IG=LENGTH_DIST_IG_F0+ELONGATION-ELONGATION_PEVK     

                        DO WHILE(ABS(DELTA_F).GT.TOL) !Newton-method (solve FORCE_DISTAL_IG-FORCE_0=0)

                          IF(LENGTH_DIST_IG.GE.LENGTHS_DIST_IG(DIM_DATA)) THEN  !Extrapolation if LENGTHS_DIST_IG(end)>LENGTH_DIST_IG
                            INDEX_I=DIM_DATA
                          ELSE 
                            DO WHILE(LENGTH_DIST_IG.GT.LENGTHS_DIST_IG(INDEX_I))
                              INDEX_I=INDEX_I+1
                            ENDDO
                          ENDIF
                          STIFFNESS_DIST=(FORCES_DIST_IG(INDEX_I)-FORCES_DIST_IG(INDEX_I-1))/ &
                            & (LENGTHS_DIST_IG(INDEX_I)-LENGTHS_DIST_IG(INDEX_I-1))
                          FORCE_DISTAL_IG=FORCES_DIST_IG(INDEX_I-1)+STIFFNESS_DIST* &
                            & (LENGTH_DIST_IG-LENGTHS_DIST_IG(INDEX_I-1))                

                          FORCE=STIFFNESS_PEVK*ELONGATION_PEVK
                          TITIN_BOUND=FORCE+F0

                          DELTA_F=TITIN_BOUND-FORCE_DISTAL_IG !new iteration for FORCE_DISTAL_IG-TITIN_BOUND
                          DIFF_QUOT=STIFFNESS_PEVK+STIFFNESS_DIST                 
                          ELONGATION_PEVK=ELONGATION_PEVK-DELTA_F/DIFF_QUOT                 
                          LENGTH_DIST_IG=LENGTH_DIST_IG_F0+ELONGATION-ELONGATION_PEVK       
                        ENDDO
                      ENDIF ! switch model

                      ! calculate x-fibre titin stress                       
                      ! fitting function to calculate the filament-lattice parameter d10 in nanometer (from Elliot 1963 -mammalian muscle)
                      d10=-13.39_DP*SARCO_LENGTH+58.37_DP
                      ! calculate the distance between actin and myosin filament in micro meter (geometrical relationship)
                      ACTIN_MYOSIN_DISTANCE=0.001_DP*(2.0_DP/3.0_DP*d10) 
                      ! calculate x-fibre stress with tangent (unbound titin)
                      TITIN_XF_UNBOUND=0.5_DP*TITIN_UNBOUND*ACTIN_MYOSIN_DISTANCE/LENGTH_TITIN 
                      ! calculate x-fibre stress with tangent (for bound titin filaments)
                      TITIN_XF_BOUND=0.5_DP*TITIN_BOUND*ACTIN_MYOSIN_DISTANCE/LENGTH_DIST_IG  
                      
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,2,TITIN_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,3,TITIN_BOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,4,TITIN_XF_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,5,TITIN_XF_BOUND,ERR,ERROR,*999)
                    ENDIF ! Check if elongation is positive or not
                  ENDDO ! Over the nodes

                  !now the ghost elements -- get the relevant info from the other computational nodes
                  CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD_MONODOMAIN, & 
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD_MONODOMAIN, & 
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        CASE DEFAULT
          CALL FlagError("Problem subtype not implemented for titin",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Problem not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN",ERR,ERROR)
    RETURN 1

  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN

  !
  !================================================================================================================================
  !

END MODULE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES

