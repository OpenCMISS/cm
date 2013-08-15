!> \file
!> \authors Andreas Hessenthaler
!> \brief This module handles all routines pertaining to finite elasticity coupled with navier stokes for fsi problems
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

!>This module handles all routines pertaining to finite elasticity coupled with navier stokes for fsi problems


MODULE FSI_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
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
  USE NAVIER_STOKES_EQUATIONS_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TYPES

  IMPLICIT NONE

  PUBLIC FSI_EQUATIONS_SET_SETUP
  PUBLIC FSI_EQUATIONS_SET_SUBTYPE_SET
  PUBLIC FSI_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC FSI_PROBLEM_SETUP
  PUBLIC FSI_PROBLEM_SUBTYPE_SET
  
  PUBLIC FSI_FINITE_ELEMENT_CALCULATE

  PUBLIC FSI_PRE_SOLVE
  PUBLIC FSI_POST_SOLVE

  PUBLIC FSI_CONTROL_LOOP_PRE_LOOP
  PUBLIC FSI_CONTROL_LOOP_POST_LOOP

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity navier stokes equation type of a multi physics equations set class.
  SUBROUTINE FSI_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,Err,Error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FSI_EQUATIONS_SET_SOLUTION_METHOD_SET",Err,Error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",Err,Error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",Err,Error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",Err,Error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",Err,Error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",Err,Error,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",Err,Error))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",Err,Error))// &
          & " is not valid for a finite elasticity navier stokes equation type of a multi physics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",Err,Error,*999)
    ENDIF
       
    CALL EXITS("FSI_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("FSI_EQUATIONS_SET_SOLUTION_METHOD_SET",Err,Error)
    CALL EXITS("FSI_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE FSI_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity navier stokes equation.
  SUBROUTINE FSI_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,Err,Error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string


    CALL ENTERS("FSI_EQUATIONS_SET_SETUP",Err,Error,*999)
    
    CALL FLAG_ERROR("FSI_EQUATIONS_SET_SETUP is not implemented.",Err,Error,*999)

    CALL EXITS("FSI_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("FSI_EQUATIONS_SET_SETUP",Err,Error)
    CALL EXITS("FSI_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE FSI_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a finite elasticity navier stokes equation finite element equations set.
  SUBROUTINE FSI_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,Err,Error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string

    CALL ENTERS("FSI_FINITE_ELEMENT_CALCULATE",Err,Error,*999)

    CALL FLAG_ERROR("FSI_FINITE_ELEMENT_CALCULATE is not implemented.",Err,Error,*999)

    CALL EXITS("FSI_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("FSI_FINITE_ELEMENT_CALCULATE",Err,Error)
    CALL EXITS("FSI_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE FSI_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a finite elasticity navier stokes  equation type of a multi physics equations set class.
  SUBROUTINE FSI_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,Err,Error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string
    !Local Variables

    CALL ENTERS("FSI_EQUATIONS_SET_SUBTYPE_SET",Err,Error,*999)

    CALL FLAG_ERROR("FSI_EQUATIONS_SET_SUBTYPE_SET is not implemented.",Err,Error,*999)

    CALL EXITS("FSI_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("FSI_EQUATIONS_SET_SUBTYPE_SET",Err,Error)
    CALL EXITS("FSI_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE FSI_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a finite elasticity navier stokes equation type .
  SUBROUTINE FSI_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,Err,Error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FSI_PROBLEM_SUBTYPE_SET",Err,Error,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE
        PROBLEM%SUBTYPE=PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",Err,Error))// &
          & " is not valid for a finite elasticity navier stokes type of a multi physics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",Err,Error,*999)
    ENDIF
       
    CALL EXITS("FSI_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("FSI_PROBLEM_SUBTYPE_SET",Err,Error)
    CALL EXITS("FSI_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE FSI_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity navier stokes equations problem.
  SUBROUTINE FSI_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,Err,Error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: MovingMeshSubLoop,ElasticitySubLoop
    TYPE(SOLVER_TYPE), POINTER :: SOLVER,MOVING_MESH_SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,MovingMeshSolverEquations,MOVING_MESH_SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: MovingMeshSolvers,SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FSI_PROBLEM_SETUP",Err,Error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(MovingMeshSubLoop)
    NULLIFY(ElasticitySubLoop)
    NULLIFY(SOLVER)
    NULLIFY(MOVING_MESH_SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(MOVING_MESH_SOLVER_EQUATIONS)
    NULLIFY(MovingMeshSolverEquations)
    NULLIFY(SOLVERS)
    NULLIFY(MovingMeshSolvers)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
        !Standard FiniteElasticity NavierStokes ALE
        CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",Err,Error))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",Err,Error))// &
                    & " is invalid for an finite elasticity ALE navier stokes equation."
                  CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a time control loop as parent loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,Err,Error,*999)
                  CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,Err,Error,*999)
                  CALL CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,Err,Error,*999)       
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,Err,Error,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,Err,Error,*999)            
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",Err,Error))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",Err,Error))// &
                    & " is invalid for a finite elasticity navier stokes equation."
                  CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,Err,Error,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                
                  !Start the solvers creation
                  CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_NUMBER_SET(SOLVERS,2,ERR,ERROR,*999)
                  !Set the first solver to be a linear solver for the Laplace mesh movement problem
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,MOVING_MESH_SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_TYPE_SET(MOVING_MESH_SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
                  !Set solver defaults
                  CALL SOLVER_LIBRARY_TYPE_SET(MOVING_MESH_SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                  !Set the solver to be a first order dynamic solver 
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                  !Set solver defaults
                  CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
                  CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
                
                
                
                
                
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                
                
                  !Get the solvers
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  !Finish the solvers creation
                  CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                  
                  
                  
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",Err,Error))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",Err,Error))// &
                      & " is invalid for a finite elasticity navier stokes equation."
                  CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,Err,Error,*999)
                  
                  
                  
                  
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,MOVING_MESH_SOLVER,ERR,ERROR,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(MOVING_MESH_SOLVER,MOVING_MESH_SOLVER_EQUATIONS,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(MOVING_MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(MOVING_MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC, &
                    & ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(MOVING_MESH_SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                  & ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                  
                  
                  
                  
                  
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,Err,Error,*999)
                  
                  
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,MOVING_MESH_SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(MOVING_MESH_SOLVER,MOVING_MESH_SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(MOVING_MESH_SOLVER_EQUATIONS,ERR,ERROR,*999)             

                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)    
                  
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",Err,Error))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",Err,Error))// &
                    & " is invalid for a finite elasticity navier stokes equation."
                  CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",Err,Error))// &
                & " is invalid for a finite elasticity ALE navier stokes equation."
              CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",Err,Error))// &
            & " does not equal a standard finite elasticity navier stokes equation subtype."
          CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",Err,Error,*999)
    ENDIF
       
    CALL EXITS("FSI_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("FSI_PROBLEM_SETUP",Err,Error)
    CALL EXITS("FSI_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE FSI_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity navier stokes problem pre-solve.
  SUBROUTINE FSI_PRE_SOLVE(ControlLoop,Solver,Err,Error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ControlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: Solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FSI_PRE_SOLVE",Err,Error,*999)

    IF(ASSOCIATED(ControlLoop)) THEN
      IF(ASSOCIATED(Solver)) THEN
        IF(ASSOCIATED(ControlLoop%PROBLEM)) THEN
          SELECT CASE(ControlLoop%PROBLEM%SUBTYPE)
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
              IF(ControlLoop%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Running pre-solve steps.",Err,Error,*999)
                !Pre solve for ALE NavierStokes equations set
                CALL NAVIER_STOKES_PRE_SOLVE(ControlLoop,Solver,Err,Error,*999)
                !Pre solve for FiniteElasticity equations set
                !Nothing to be done???
              ELSE
                CALL FLAG_ERROR("Incorrect loop type. Must be time loop.",Err,Error,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(ControlLoop%PROBLEM%SUBTYPE,"*",Err,Error))// &
                & " is not valid for a finite elasticity navier stokes type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",Err,Error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",Err,Error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",Err,Error,*999)
    ENDIF

    CALL EXITS("FSI_PRE_SOLVE")
    RETURN
999 CALL ERRORS("FSI_PRE_SOLVE",Err,Error)
    CALL EXITS("FSI_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE FSI_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the finite elasticity navier stokes  problem post solve.
  SUBROUTINE FSI_POST_SOLVE(ControlLoop,Solver,Err,Error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ControlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: Solver,Solver2!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR,Method,FileName

    CALL ENTERS("FSI_POST_SOLVE",Err,Error,*999)
    
    NULLIFY(Solver2)

    IF(ASSOCIATED(ControlLoop)) THEN
      IF(ASSOCIATED(Solver)) THEN
        IF(ASSOCIATED(ControlLoop%PROBLEM)) THEN 
          SELECT CASE(ControlLoop%PROBLEM%SUBTYPE)
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
              !Post solve for the linear solver
              IF(Solver%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement post solve... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(Solver%SOLVERS,1,Solver2,ERR,ERROR,*999)
           !     CALL SOLVERS_SOLVER_GET(Solver%SOLVERS,2,Solver2,ERR,ERROR,*999)
                IF(ASSOCIATED(Solver2%DYNAMIC_SOLVER)) THEN
                  Solver2%DYNAMIC_SOLVER%ALE=.TRUE.
                ELSE  
                  CALL FLAG_ERROR("Dynamic solver is not associated for ALE problem.",ERR,ERROR,*999)
                END IF
                
                
                
                
                
                
                
                
              !Post solve for the dynamic solver
              ELSE IF(Solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Navier-Stokes post solve... ",ERR,ERROR,*999)
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
     !           IF(.NOT.ASSOCIATED(ControlLoop%TIME_LOOP)) CALL FLAG_ERROR("Time loop is not associated.",Err,Error,*999)
     !           !Export solid fields
     !           FileName="SolidStep00"//TRIM(NUMBER_TO_VSTRING( &
     !             & INT(ControlLoop%TIME_LOOP%CURRENT_TIME/ControlLoop%TIME_LOOP%TIME_INCREMENT),"*",Err,Error))// &
     !             & "output"
     !           Method="FORTRAN"
     !           CALL FIELD_IO_NODES_EXPORT(Solver%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%REGION%FIELDS, &
     !             & FileName,Method,ERR,ERROR,*999)
     !           CALL FIELD_IO_ELEMENTS_EXPORT(Solver%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%REGION%FIELDS, &
     !             & FileName,Method,ERR,ERROR,*999)
     !           !Export fluid fields
     !           FileName="FluidStep00"//TRIM(NUMBER_TO_VSTRING( &
     !             & INT(ControlLoop%TIME_LOOP%CURRENT_TIME/ControlLoop%TIME_LOOP%TIME_INCREMENT),"*",Err,Error))// &
     !             & "output"
     !           Method="FORTRAN"
     !           CALL FIELD_IO_NODES_EXPORT(Solver%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(2)%PTR%REGION%FIELDS, &
     !             & FileName,Method,ERR,ERROR,*999)
     !           CALL FIELD_IO_ELEMENTS_EXPORT(Solver%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(2)%PTR%REGION%FIELDS, &
     !             & FileName,Method,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(ControlLoop%PROBLEM%SUBTYPE,"*",Err,Error))// &
                  & " for a FiniteElasticity-NavierStokes type of a multi physics problem class has unknown solver solve type."
                CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(ControlLoop%PROBLEM%SUBTYPE,"*",Err,Error))// &
                & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",Err,Error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",Err,Error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",Err,Error,*999)
    ENDIF

    CALL EXITS("FSI_POST_SOLVE")
    RETURN
999 CALL ERRORS("FSI_POST_SOLVE",Err,Error)
    CALL EXITS("FSI_POST_SOLVE")
    RETURN 1
  END SUBROUTINE FSI_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE FSI_CONTROL_LOOP_PRE_LOOP(ControlLoop,Err,Error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ControlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string

    !Local Variables

    CALL ENTERS("FSI_CONTROL_LOOP_PRE_LOOP",Err,Error,*999)
    
    CALL FLAG_ERROR("FSI_CONTROL_LOOP_PRE_LOOP not implemented.",Err,Error,*999)

    CALL EXITS("FSI_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("FSI_CONTROL_LOOP_PRE_LOOP",Err,Error)
    CALL EXITS("FSI_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE FSI_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration. Updates interface and fluid geometric fields and exports fields.
  SUBROUTINE FSI_CONTROL_LOOP_POST_LOOP(ControlLoop,Err,Error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ControlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: DynamicSolver,LinearSolver!<A pointer to the solvers
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TimeLoop
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: DynamicSolverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: DynamicSolverMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: SolidEquationsSet,FluidEquationsSet,EquationsSet
    TYPE(FIELD_TYPE), POINTER :: SolidGeometricField,InterfaceGeometricField,SolidDependentField
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: InterfaceCondition
    TYPE(INTERFACE_TYPE), POINTER :: FSInterface
    TYPE(NODES_TYPE), POINTER :: InterfaceNodes
    TYPE(VARYING_STRING) :: Method,SolidFileName,FluidFileName,InterfaceFileName
    REAL(DP) :: StartTime,CurrentTime,TimeIncrement,TimeStepNumber,Value
    INTEGER(INTG) :: EquationsSetIndex,InterfaceNodeNumber,InterfaceNodeComponent,NumberOfComponents
    LOGICAL :: FluidEquationsSetFound,SolidEquationsSetFound=.FALSE.

    CALL ENTERS("FSI_CONTROL_LOOP_POST_LOOP",Err,Error,*999)
    
    NULLIFY(DynamicSolver)
    NULLIFY(LinearSolver)
    NULLIFY(TimeLoop)
    NULLIFY(DynamicSolverEquations)
    NULLIFY(SolidEquationsSet)
    NULLIFY(FluidEquationsSet)
    NULLIFY(SolidGeometricField)
    NULLIFY(InterfaceGeometricField)
    NULLIFY(SolidDependentField)
    NULLIFY(InterfaceCondition)
    NULLIFY(FSInterface)
    NULLIFY(InterfaceNodes)
    
    !Check pointers
    IF(.NOT.ASSOCIATED(ControlLoop)) CALL FLAG_ERROR("Main control loop not associated.",Err,Error,*999)
    IF(.NOT.ASSOCIATED(ControlLoop%SOLVERS)) CALL FLAG_ERROR("Solvers are not associated.",Err,Error,*999)
    !Get solvers for FSI
    CALL SOLVERS_SOLVER_GET(ControlLoop%SOLVERS,1,DynamicSolver,Err,Error,*999)
    CALL SOLVERS_SOLVER_GET(ControlLoop%SOLVERS,2,LinearSolver,Err,Error,*999)
    TimeLoop=>ControlLoop%TIME_LOOP
    IF(.NOT.ASSOCIATED(TimeLoop)) CALL FLAG_ERROR("Time loop not associated.",Err,Error,*999)
    !Get times
    StartTime=TimeLoop%START_TIME
    CurrentTime=TimeLoop%CURRENT_TIME
    TimeIncrement=TimeLoop%TIME_INCREMENT
    TimeStepNumber=(CurrentTime-StartTime)/TimeIncrement!GLOBAL_ITERATION_NUMBER???
    !===============================================================================================================================
    !First update mesh and calculate boundary velocity values
    CALL NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH(ControlLoop,DynamicSolver,Err,Error,*999)
    !===============================================================================================================================
    !Update interface geometric field and export results
    DynamicSolverEquations=>DynamicSolver%SOLVER_EQUATIONS
    IF(ASSOCIATED(DynamicSolverEquations)) THEN
      DynamicSolverMapping=>DynamicSolverEquations%SOLVER_MAPPING
      IF(ASSOCIATED(DynamicSolverMapping)) THEN
        EquationsSetIndex=1
        FluidEquationsSetFound=.FALSE.
        SolidEquationsSetFound=.FALSE.
        DO WHILE((EquationsSetIndex<=DynamicSolverMapping%NUMBER_OF_EQUATIONS_SETS &
          & .AND..NOT.SolidEquationsSetFound) &
          & .OR.(EquationsSetIndex<=DynamicSolverMapping%NUMBER_OF_EQUATIONS_SETS &
          & .AND..NOT.FluidEquationsSetFound))
          EquationsSet=>DynamicSolverMapping%EQUATIONS_SETS(EquationsSetIndex)%PTR
          IF(EquationsSet%CLASS==EQUATIONS_SET_ELASTICITY_CLASS &
            & .AND.EquationsSet%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE &
            & .AND.EquationsSet%SUBTYPE==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE) THEN
            SolidEquationsSet=>EquationsSet
            SolidEquationsSetFound=.TRUE.
          ELSEIF(EquationsSet%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
            & .AND.EquationsSet%TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
            & .AND.EquationsSet%SUBTYPE==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE) THEN
            FluidEquationsSet=>EquationsSet
            FluidEquationsSetFound=.TRUE.
          ELSE
            CALL FLAG_ERROR("Invalid equations sets associated with dynamic solver for FSI.", Err,Error,*999)
          ENDIF
          EquationsSetIndex=EquationsSetIndex+1
        ENDDO
        IF(.NOT.SolidEquationsSetFound) CALL FLAG_ERROR("Could not find solid equations set for FSI.",Err,Error,*999)
        IF(.NOT.FluidEquationsSetFound) CALL FLAG_ERROR("Could not find fluid equations set for FSI.",Err,Error,*999)
        SolidGeometricField=>SolidEquationsSet%GEOMETRY%GEOMETRIC_FIELD
        IF(ASSOCIATED(SolidGeometricField)) THEN
          CALL FIELD_NUMBER_OF_COMPONENTS_GET(SolidGeometricField,FIELD_U_VARIABLE_TYPE,NumberOfComponents,Err,Error,*999)
          IF(DynamicSolverMapping%NUMBER_OF_INTERFACE_CONDITIONS>1) THEN
            CALL FLAG_ERROR("Invalid number of interface conditions. Must be 1 for FSI.",Err,Error,*999)
          ENDIF
          SolidDependentField=>SolidEquationsSet%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(SolidDependentField)) THEN
            InterfaceCondition=>DynamicSolverMapping%INTERFACE_CONDITIONS(1)%PTR
            IF(ASSOCIATED(InterfaceCondition)) THEN
              FSInterface=>InterfaceCondition%INTERFACE
              IF(ASSOCIATED(FSInterface)) THEN
                InterfaceNodes=>FSInterface%NODES
                IF(ASSOCIATED(InterfaceNodes)) THEN
                  InterfaceGeometricField=>InterfaceCondition%GEOMETRY%GEOMETRIC_FIELD
                  IF(ASSOCIATED(InterfaceGeometricField)) THEN
                    !===============================================================================================================
                    !Update interface geometric field
                    DO InterfaceNodeNumber=1,InterfaceNodes%NUMBER_OF_NODES
                      DO InterfaceNodeComponent=1,NumberOfComponents
                        !Default to version 1, derivative 1
                        CALL FIELD_PARAMETER_SET_GET_NODE(SolidDependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                          & 1,1,InterfaceNodes%COUPLED_NODES(1,InterfaceNodeNumber),InterfaceNodeComponent,Value, &
                          & Err,Error,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_NODE(InterfaceGeometricField,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,InterfaceNodeNumber,InterfaceNodeComponent,Value,Err,Error,*999)
                      ENDDO
                    ENDDO
                    CALL FIELD_PARAMETER_SET_UPDATE_START(InterfaceGeometricField, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,Err,Error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(InterfaceGeometricField, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,Err,Error,*999)
                    !===============================================================================================================
                    !Export fields
                    IF(TimeStepNumber>=0.AND.TimeStepNumber<10) THEN
                      SolidFileName="Solid000"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      FluidFileName="Fluid000"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      InterfaceFileName="Interface000"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                    ELSE IF(TimeStepNumber>=10.AND.TimeStepNumber<100) THEN
                      SolidFileName="Solid00"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      FluidFileName="Fluid00"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      InterfaceFileName="Interface00"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                    ELSE IF(TimeStepNumber>=100.AND.TimeStepNumber<1000) THEN
                      SolidFileName="Solid0"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      FluidFileName="Fluid0"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      InterfaceFileName="Interface0"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                    ELSE IF(TimeStepNumber>=1000) THEN
                      SolidFileName="Solid"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      FluidFileName="Fluid"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                      InterfaceFileName="Interface"//TRIM(NUMBER_TO_VSTRING(INT(TimeStepNumber),"*",Err,Error))
                    ENDIF
                    Method="FORTRAN"
                    !Export solid fields
                    IF(.NOT.ASSOCIATED(SolidEquationsSet%REGION)) CALL FLAG_ERROR("Solid region not associated.", &
                      & Err,Error,*999)
                    IF(.NOT.ASSOCIATED(SolidEquationsSet%REGION%FIELDS)) CALL FLAG_ERROR("Solid fields not associated.", &
                      & Err,Error,*999)
                    CALL FIELD_IO_NODES_EXPORT(SolidEquationsSet%REGION%FIELDS, &
                      & SolidFileName,Method,Err,Error,*999)
                    CALL FIELD_IO_ELEMENTS_EXPORT(SolidEquationsSet%REGION%FIELDS, &
                      & SolidFileName,Method,Err,Error,*999)
                    IF(.NOT.ASSOCIATED(FluidEquationsSet%REGION)) CALL FLAG_ERROR("Fluid region not associated.", &
                      & Err,Error,*999)
                    IF(.NOT.ASSOCIATED(FluidEquationsSet%REGION%FIELDS)) CALL FLAG_ERROR("Fluid fields not associated.", &
                      & Err,Error,*999)
                    !Export fluid fields
                    CALL FIELD_IO_NODES_EXPORT(FluidEquationsSet%REGION%FIELDS, &
                      & FluidFileName,Method,Err,Error,*999)
                    CALL FIELD_IO_ELEMENTS_EXPORT(FluidEquationsSet%REGION%FIELDS, &
                      & FluidFileName,Method,Err,Error,*999)
                    IF(.NOT.ASSOCIATED(FSInterface%FIELDS)) CALL FLAG_ERROR("Interface fields not associated.",Err,Error,*999)
                    !Export interface fields
                    CALL FIELD_IO_NODES_EXPORT(FSInterface%FIELDS,InterfaceFileName,Method,Err,Error,*999)
                    CALL FIELD_IO_ELEMENTS_EXPORT(FSInterface%FIELDS,InterfaceFileName,Method,Err,Error,*999)
                  ELSE
                    CALL FLAG_ERROR("Interface geometric field not associated.",Err,Error,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Interface nodes not associated.",Err,Error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface not associated.",Err,Error,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface condition not associated.",Err,Error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solid dependent field not associated.",Err,Error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solid geometric field not associated.",Err,Error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dynamic solver mapping not associated.",Err,Error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Dynamic solver equations not associated.",Err,Error,*999)
    ENDIF

    CALL EXITS("FSI_CONTROL_LOOP_POST_LOOP")
    RETURN
999 CALL ERRORS("FSI_CONTROL_LOOP_POST_LOOP",Err,Error)
    CALL EXITS("FSI_CONTROL_LOOP_POST_LOOP")
    RETURN 1
  END SUBROUTINE FSI_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

END MODULE FSI_ROUTINES
