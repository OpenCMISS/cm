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
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TYPES

  IMPLICIT NONE

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET
  
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectrics finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE)
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
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectrics finite elasticity equation type of a multi physics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity equation.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string


    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    CALL FLAG_ERROR("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP is not implemented.",ERR,ERROR,*999)

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a bioelectrics finite elasticity equation finite element equations set.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    CALL FLAG_ERROR("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE is not implemented.",ERR,ERROR,*999)

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a bioelectrics finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)

    CALL FLAG_ERROR("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET is not implemented.",ERR,ERROR,*999)

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a bioelectric finite elasticity problem type .
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric finite elasticity problem type of a multi physics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure equations problem.
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

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(MONODOMAIN_SUB_LOOP)
    NULLIFY(ELASTICITY_SUB_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(MONODOMAIN_SOLVERS)
    NULLIFY(ELASTICITY_SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)

      !--------------------------------------------------------------------
      !   Transient Gudunov monodomain, simple finite elasticity  
      !--------------------------------------------------------------------
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)
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
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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

            !Set up the control sub loop for finite elasicity
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(ELASTICITY_SUB_LOOP,'ELASTICITY_SIMPLE_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(ELASTICITY_SUB_LOOP,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(ELASTICITY_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
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
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectrics finite elasticity equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a transient monodomain quasistatic finite elasticity equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP")
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

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              CALL BIODOMAIN_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
              CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
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

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE")
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

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)
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
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
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

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration???
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    CALL FLAG_ERROR("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP is not implemented.",ERR,ERROR,*999)
    
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

END MODULE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES
