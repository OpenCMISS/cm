!> \file
!> $Id$
!> \authors Christian Michler, Jack Lee
!> \brief This module handles all routines pertaining to finite elasticity coupled with Darcy.
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

!>This module handles all routines pertaining to finite elasticity coupled with Darcy.


MODULE FINITE_ELASTICITY_DARCY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES  
  USE DARCY_EQUATIONS_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE FINITE_ELASTICITY_ROUTINES
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

  PUBLIC ELASTICITY_DARCY_EQUATIONS_SET_SETUP
  PUBLIC ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET
  PUBLIC ELASTICITY_DARCY_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC ELASTICITY_DARCY_PROBLEM_SETUP
  PUBLIC ELASTICITY_DARCY_PROBLEM_SUBTYPE_SET
  
  PUBLIC ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE

  PUBLIC ELASTICITY_DARCY_PRE_SOLVE
  PUBLIC ELASTICITY_DARCY_POST_SOLVE

  PUBLIC ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP
  PUBLIC ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity Darcy equation type of a multi physics equations set class.
  SUBROUTINE ELASTICITY_DARCY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_DARCY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_ELASTICITY_DARCY_SUBTYPE)
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
          & " is not valid for a finite elasticity Darcy  equation type of a multi physics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_DARCY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy equation.
  SUBROUTINE ELASTICITY_DARCY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string


    CALL ENTERS("ELASTICITY_DARCY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    CALL FLAG_ERROR("ELASTICITY_DARCY_EQUATIONS_SET_SETUP still needs to be implemented.",ERR,ERROR,*999)

    !=================================================================
    ! This routine still needs to be implemented.
    ! It will be used to setup the equations set of a monolithic
    ! finite-elasticity Darcy system.
    ! For the partitioned solution this routine is not called,
    ! since EQUATIONS_SET_SETUP of respective equations_set is called.
    !=================================================================

             
    CALL EXITS("ELASTICITY_DARCY_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a finite elasticity Darcy equation finite element equations set.
  SUBROUTINE ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string


    CALL ENTERS("ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    CALL FLAG_ERROR("ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE still needs to be implemented.",ERR,ERROR,*999)

    !=================================================================
    ! This routine still needs to be implemented.
    ! It will be used to calculate the finite-element matrices and vectors
    ! of a monolithic finite-elasticity Darcy system.
    ! For the partitioned solution this routine is not called,
    ! since FINITE_ELEMENT_CALCULATE of respective equations_set is called.
    !=================================================================


  CALL EXITS("ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a finite elasticity Darcy  equation type of a multi physics equations set class.
  SUBROUTINE ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    

    CALL ENTERS("ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)

    CALL FLAG_ERROR("ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET still needs to be implemented.",ERR,ERROR,*999)

    !=================================================================
    ! This routine still needs to be implemented.
    ! It will be used to set the equations_set_subtype
    ! of a monolithic finite-elasticity Darcy system.
    ! For the partitioned solution this routine is not called,
    ! since EQUATIONS_SET_SUBTYPE_SET of respective equations_set is called.
    !=================================================================

       
    CALL EXITS("ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a finite elasticity Darcy equation type .
  SUBROUTINE ELASTICITY_DARCY_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_DARCY_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_DARCY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE     
      CASE(PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_DARCY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE     
      CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_DARCY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE     
      CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_DARCY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE     
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity Darcy equation type of a multi physics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_DARCY_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy equations problem.
  SUBROUTINE ELASTICITY_DARCY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT,SOLID_SUB_LOOP,FLUID_SUB_LOOP,SUBITERATION_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER, SOLVER_MAT_PROPERTIES, SOLVER_SOLID
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS, SOLVER_EQUATIONS_MAT_PROPERTIES, SOLVER_EQUATIONS_SOLID
    TYPE(SOLVERS_TYPE), POINTER :: SOLID_SOLVERS,FLUID_SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_DARCY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SUBITERATION_LOOP)
    NULLIFY(SOLID_SUB_LOOP)
    NULLIFY(FLUID_SUB_LOOP)
    NULLIFY(SOLID_SOLVERS)
    NULLIFY(FLUID_SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_MAT_PROPERTIES)
    NULLIFY(SOLVER_SOLID)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_EQUATIONS_MAT_PROPERTIES)
    NULLIFY(SOLVER_EQUATIONS_SOLID)
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)

      !--------------------------------------------------------------------
      !   s t a n d a r d   f i n i t e   e l a s t i c i t y   D a r c y
      !--------------------------------------------------------------------
      CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
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
              & " is invalid for an finite elasticity ALE Darcy  equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
            !Solid, load incremented control loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(SOLID_SUB_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
            !Fluid control loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(FLUID_SUB_LOOP,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
            !Sub-loops are finished when parent is finished
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation for the solid solver
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL SOLVERS_CREATE_START(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLID_SOLVERS,1,ERR,ERROR,*999)
            !
            !Set the first solver to be a nonlinear solver for the finite elasticity
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_SOLID,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_SOLID,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
            !
            !Start the solvers creation for the fluid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL SOLVERS_CREATE_START(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(FLUID_SOLVERS,2,ERR,ERROR,*999)
            !
            !Set the second solver to be a linear solver for the material update
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
            !
            !Set the third solver to be a linear solver for ALE Darcy
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLID_SOLVERS,ERR,ERROR,*999)
            !Get the fluid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(FLUID_SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !
            !Get the finite elasticity solver and create the finite elasticity solver equations
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_SOLID,SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !
            !Get the material-properties solver and create the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_QUASISTATIC, &
              & ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the finite elasticity solver equations
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_SOLID,SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)             
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_MAT_PROPERTIES,ERR,ERROR,*999)             
            !
            !Finish the creation of the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !--------------------------------------------------------------------
      !   q u a s i s t a t i c   f i n i t e   e l a s t i c i t y   t r a n s i e n t   D a r c y
      !--------------------------------------------------------------------
      CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
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
              & " is invalid for an finite elasticity ALE Darcy  equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up a subiteration loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(SUBITERATION_LOOP,'SUBITERATION_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(SUBITERATION_LOOP,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(SUBITERATION_LOOP,100,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(SUBITERATION_LOOP,2,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(SUBITERATION_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up load incremented control loop for Solid
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(SOLID_SUB_LOOP,'FINITE_ELASTICITY_LOAD_INCREMENT_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(SOLID_SUB_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
            !For problems that require it, the user can get the solid subloop using:
            !CALL CMISSProblemControlLoopGet(Problem,[1,1,CMISSControlLoopNode],ControlLoopSolid,Err)
            !And then set the number of load increments to 3 for example with:
            !CALL CMISSControlLoopMaximumIterationsSet(ControlLoopSolid,3,Err)
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(SOLID_SUB_LOOP,1,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(SOLID_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up control loop for Fluid
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(FLUID_SUB_LOOP,'DARCY_SIMPLE_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(FLUID_SUB_LOOP,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(FLUID_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
            !Sub-loops are finished when parent is finished
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation for the solid solver
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL SOLVERS_CREATE_START(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLID_SOLVERS,1,ERR,ERROR,*999)
            !
            !Set the solid solver to be a nonlinear solver for the finite elasticity
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_SOLID,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_SOLID,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
            !
            !Start the solvers creation for the fluid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL SOLVERS_CREATE_START(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(FLUID_SOLVERS,1,ERR,ERROR,*999)
            !
            !Set the solver to be a first-order dynamic solver for Darcy
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLID_SOLVERS,ERR,ERROR,*999)

            !Get the fluid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(FLUID_SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !
            !Get the finite elasticity solver and create the finite elasticity solver equations
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_SOLID,SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !
            !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the finite elasticity solver equations
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_SOLID,SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)             
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !--------------------------------------------------------------------
      !   q u a s i s t a t i c   e l a s t i c i t y   t r a n s i e n t   D a r c y   M A T E R I A L   S O L V E
      !--------------------------------------------------------------------
      CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
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
              & " is invalid for an finite elasticity ALE Darcy  equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up a subiteration loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(SUBITERATION_LOOP,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(SUBITERATION_LOOP,9,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(SUBITERATION_LOOP,2,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(SUBITERATION_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up load incremented control loop for Solid
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(SOLID_SUB_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
            !For problems that require it, the user can get the solid subloop using:
            !CALL CMISSProblemControlLoopGet(Problem,[1,1,CMISSControlLoopNode],ControlLoopSolid,Err)
            !And then set the number of load increments to 3 for example with:
            !CALL CMISSControlLoopMaximumIterationsSet(ControlLoopSolid,3,Err)
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(SOLID_SUB_LOOP,1,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(SOLID_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up control loop for Fluid
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(FLUID_SUB_LOOP,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(FLUID_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
            !Sub-loops are finished when parent is finished
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation for the solid solver
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL SOLVERS_CREATE_START(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLID_SOLVERS,1,ERR,ERROR,*999)
            !
            !Set the solid solver to be a nonlinear solver for the finite elasticity
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_SOLID,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_SOLID,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
            !
            !Start the solvers creation for the fluid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL SOLVERS_CREATE_START(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(FLUID_SOLVERS,2,ERR,ERROR,*999)
            !
            !Set the solver to be a linear solver for the material update
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_MAT_PROPERTIES,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
            !
            !Set the other solver to be a first-order dynamic solver for Darcy
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLID_SOLVERS,ERR,ERROR,*999)

            !Get the fluid solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(FLUID_SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !
            !Get the finite elasticity solver and create the finite elasticity solver equations
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_SOLID,SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_SOLID,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !
            !Get the material-properties solver and create the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_EQUATIONS_QUASISTATIC, &
              & ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_MAT_PROPERTIES,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUBITERATION_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,1,SOLID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(SOLID_SUB_LOOP,SOLID_SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the finite elasticity solver equations
            CALL SOLVERS_SOLVER_GET(SOLID_SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_SOLID,SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_SOLID,ERR,ERROR,*999)             
            !
            CALL CONTROL_LOOP_SUB_LOOP_GET(SUBITERATION_LOOP,2,FLUID_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(FLUID_SUB_LOOP,FLUID_SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the material-properties solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,1,SOLVER_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_MAT_PROPERTIES,SOLVER_EQUATIONS_MAT_PROPERTIES,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_MAT_PROPERTIES,ERR,ERROR,*999)             
            !
            !Finish the creation of the Darcy-ALE solver equations
            CALL SOLVERS_SOLVER_GET(FLUID_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity ALE Darcy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   c a s e   d e f a u l t
      !-----------------------------------------------------------------
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard finite elasticity Darcy equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)

      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_DARCY_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity Darcy problem pre-solve.
  SUBROUTINE ELASTICITY_DARCY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_DARCY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
              IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE.AND.SOLVER%GLOBAL_NUMBER==1) THEN
                CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ELSE IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_SIMPLE_TYPE) THEN
                IF(SOLVER%GLOBAL_NUMBER==1) THEN
!                   IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
!                     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now working on material parameters",ERR,ERROR,*999)
!                   ENDIF
                ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                  IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now working on Darcy",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                CALL DARCY_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF


            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE.AND.SOLVER%GLOBAL_NUMBER==1) THEN
                CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ELSE IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_SIMPLE_TYPE) THEN
                IF(SOLVER%GLOBAL_NUMBER==1) THEN
                  IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now working on material parameters",ERR,ERROR,*999)
                  ENDIF
                ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                  IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now working on Darcy",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                CALL DARCY_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF



            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a Darcy fluid type of a multi physics problem class."
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

    CALL EXITS("ELASTICITY_DARCY_PRE_SOLVE")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy  problem post solve.
  SUBROUTINE ELASTICITY_DARCY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_DARCY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
!               CALL ELASTICITY_DARCY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              CALL FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              CALL DARCY_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a finite elasticity Darcy type of a multi physics problem class."
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

    CALL EXITS("ELASTICITY_DARCY_POST_SOLVE")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_POST_SOLVE",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_POST_SOLVE")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DARCY
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_DARCY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_DARCY)
    NULLIFY(SOLVER_DARCY)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        ! Eventually we may want to do different things depending on problem type/subtype
        ! too but for now we can just check the loop type.
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
            IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"==================================================",ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"=============== Starting time step ===============",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"CURRENT_TIME          = ",CURRENT_TIME,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"TIME_INCREMENT        = ",TIME_INCREMENT,ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"==================================================",ERR,ERROR,*999)
            ENDIF
            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"==================================================",ERR,ERROR,*999)
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"=============== Starting time step ===============",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"CURRENT_TIME          = ",CURRENT_TIME,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"TIME_INCREMENT        = ",TIME_INCREMENT,ERR,ERROR,*999)
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"==================================================",ERR,ERROR,*999)
            ENDIF
            CALL DARCY_CONTROL_TIME_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            CALL FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)

          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            !Subiteration loop
            IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"++++ Starting subiteration ++++",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"SUBITERATION_NUMBER       =   ", &
                & CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",ERR,ERROR,*999)
            ENDIF
            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",ERR,ERROR,*999)
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"++++ Starting subiteration ++++",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"SUBITERATION_NUMBER       =   ", &
                & CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",ERR,ERROR,*999)
            ENDIF
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,(/2,CONTROL_LOOP_NODE/),CONTROL_LOOP_DARCY,ERR,ERROR,*999)

            SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
              CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
                CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,1,SOLVER_DARCY,ERR,ERROR,*999)
              CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,2,SOLVER_DARCY,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                  & " is not valid for ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT

            CALL DARCY_EQUATION_PRE_SOLVE_STORE_PREVIOUS_ITERATE(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)

          CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
            IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"-- Starting fluid solve iteration --",ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
            ENDIF
            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"-- Starting fluid solve iteration --",ERR,ERROR,*999)
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
            ENDIF

          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"-- Starting solid solve iteration --",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LOAD INCREMENT NUMBER =           ", &
                & CONTROL_LOOP%LOAD_INCREMENT_LOOP%ITERATION_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
            ENDIF
            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"-- Starting solid solve iteration --",ERR,ERROR,*999)
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

    CALL EXITS("ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DARCY
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_DARCY

    NULLIFY(SOLVER_DARCY)
    NULLIFY(CONTROL_LOOP_DARCY)

    CALL ENTERS("ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"End of time step",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
            !subiteration
            IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"End of subiteration",ERR,ERROR,*999)
            ENDIF
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,(/2,CONTROL_LOOP_NODE/),CONTROL_LOOP_DARCY,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,1,SOLVER_DARCY,ERR,ERROR,*999)
            !CALL DARCY_EQUATION_POST_SOLVE_DETERMINE_DARCY_VELOCITY(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
            !CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
            CALL DARCY_EQUATION_MONITOR_CONVERGENCE(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
            !CALL DARCY_EQUATION_ACCELERATE_CONVERGENCE(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
          CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
            !subiteration
            IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"End of subiteration",ERR,ERROR,*999)
            ENDIF
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,(/2,CONTROL_LOOP_NODE/),CONTROL_LOOP_DARCY,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_DARCY%SOLVERS,2,SOLVER_DARCY,ERR,ERROR,*999)
            !CALL DARCY_EQUATION_POST_SOLVE_DETERMINE_DARCY_VELOCITY(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
            !CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
            CALL DARCY_EQUATION_MONITOR_CONVERGENCE(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
            !CALL DARCY_EQUATION_ACCELERATE_CONVERGENCE(CONTROL_LOOP,SOLVER_DARCY,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a Darcy fluid type of a multi physics problem class with a while control loop."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"End of fluid solve iteration",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"End of solid solve iteration",ERR,ERROR,*999)
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

    CALL EXITS("ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy problem post solve output data.
  SUBROUTINE ELASTICITY_DARCY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_DARCY_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)

              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                CALL FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ELSE IF(SOLVER%GLOBAL_NUMBER==2.OR.SOLVER%GLOBAL_NUMBER==3) THEN
                CALL DARCY_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF

            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a Darcy fluid type of a multi physics problem class."
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
      
    CALL EXITS("ELASTICITY_DARCY_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 CALL ERRORS("ELASTICITY_DARCY_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    CALL EXITS("ELASTICITY_DARCY_POST_SOLVE_OUTPUT_DATA")
    RETURN 1
  END SUBROUTINE ELASTICITY_DARCY_POST_SOLVE_OUTPUT_DATA
      
  !   
  !================================================================================================================================
  !


END MODULE FINITE_ELASTICITY_DARCY_ROUTINES
