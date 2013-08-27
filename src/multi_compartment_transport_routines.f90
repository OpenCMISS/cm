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
!> The Original Code is openCMISS
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

!>TThis module handles all routines pertaining to (advection-)diffusion coupled to (advection-)diffusion.


MODULE MULTI_COMPARTMENT_TRANSPORT_ROUTINES

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

  PUBLIC MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SETUP
  PUBLIC MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SUBTYPE_SET
  PUBLIC MULTI_COMPARTMENT_TRANSPORT_SET_SOLUTION_METHOD_SET

  PUBLIC MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP
  PUBLIC MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SUBTYPE_SET
  
  PUBLIC MULTI_COMPARTMENT_TRANSPORT_FINITE_ELEMENT_CALCULATE

  PUBLIC MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE
  PUBLIC MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE

  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)

              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!     IF(ASSOCIATED(EQUATIONS_SET)) THEN
!       SELECT CASE(EQUATIONS_SET%SUBTYPE)
!       CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
!         SELECT CASE(SOLUTION_METHOD)
!         CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
!           EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
!         CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
!           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!         CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
!           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!         CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
!           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!         CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
!           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!         CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
!           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!         CASE DEFAULT
!           LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
!           CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
!         END SELECT
!       CASE DEFAULT
!         LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
!           & " is not valid for a diffusion & advection-diffusion equation type of a multi physics equations set class."
!         CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
!       END SELECT
!     ELSE
!       CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
!     ENDIF
       
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets up the multi-compartment coupled advection-diffusion & diffusion transport equation.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables


    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
             
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SETUP")
    RETURN 1

  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a multi-compartment coupled advection-diffusion & diffusion transport equation.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables



    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      
  CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a multi-compartment coupled advection-diffusion & diffusion transport equation type of a multi physics equations set class.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    
    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    
    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
       
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a coupled diffusion & advection-diffusion equation type .
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE     
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a multi-compartment coupled transport equation type of a multi physics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the coupled diffusion-diffusion equations problem.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER_DIFFUSION)
    NULLIFY(SOLVER_ADVECTION_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)
     IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)

      !--------------------------------------------------------------------
      !   monolithic coupled source diffusion-diffusion problem
      !--------------------------------------------------------------------
      CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
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
              & " is invalid for a multi-compartment transport equation."
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
              & " is invalid for a multi-compartment transport equation."
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
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !Set the solver to be a linear solver for the diffusion problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_DIFFUSION,ERR,ERROR,*999)
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
                & " is invalid for a multi-compartment transport equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Get the diffusion solver and create the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_DIFFUSION,ERR,ERROR,*999)
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
            !Finish the creation of the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_DIFFUSION,SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)             
            !
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a multi-compartment transport equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a multi-compartment transport equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   c a s e   d e f a u l t
      !-----------------------------------------------------------------
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard multi-component transport equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)

      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the multi-compartment coupled transport problem pre-solve.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR


    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
             SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
             EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
             IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
               IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN

                CALL MULTI_COMP_TRANSPORT_PRE_SOLVE_UPDATE_ANALYTIC_VALUES(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
!               IF(SOLVER%GLOBAL_NUMBER==1) THEN
!                 !copy current value of concentration_one to another variable
!                 !CALL ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
!                 !Set source term to be updated value of concentration_two
!                 !CALL ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
!               ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
!                 !compute value of constant source term - evaluated from lamdba*(0.5*(c_1^{t+1}+c_1^{t}) - c_2^{t})
!                 !CALL DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
!               ENDIF
                ENDIF
               ENDIF
              ENDIF
             ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a multi-compartment transport type of a multi physics problem class."
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

    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !
  !updates the boundary conditions and source term to the required analytic values
  SUBROUTINE MULTI_COMP_TRANSPORT_PRE_SOLVE_UPDATE_ANALYTIC_VALUES(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD
!    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
!    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
!    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
!    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:),GEOMETRIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP) :: NORMAL(3),TANGENTS(3,3),VALUE,X(3),VALUE_SOURCE !<The value to add
!     REAL(DP) :: k_xx, k_yy, k_zz
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,eqnset_idx
    INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX
    REAL(DP) :: A1,A2,A3,A4,D1,D2,D3,D4,LAMBDA_12,LAMBDA_13,LAMBDA_23
!    INTEGER(INTG) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
!    INTEGER(INTG) :: DERIVATIVE_NUMBER !<The node derivative number
!    INTEGER(INTG) :: COMPONENT_NUMBER !<The field variable component number
!    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of (geometry) nodes
!    INTEGER(INTG) :: LOCAL_NODE_NUMBER
!    INTEGER(INTG) :: EQUATIONS_SET_IDX
!    INTEGER(INTG) :: equations_row_number

    CALL ENTERS("MULTI_COMP_TRANSPORT_PRE_SOLVE_UPDATE_ANALYTIC_VALUES",ERR,ERROR,*999)


    A1 = 0.4_DP
    A2 = 0.3_DP
    A3 = 0.2_DP
    A4 = 0.1_DP
    D1=1.0_DP
    D2=1.0_DP
    D3=1.0_DP
    D4=1.0_DP
    LAMBDA_12=0.1_DP
    LAMBDA_13=0.1_DP
    LAMBDA_23=0.1_DP

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
       !write(*,*)'CURRENT_TIME = ',CURRENT_TIME
       !write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            !do nothing?! 
            CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                !loop over all the equation sets and set the appropriate field variable type BCs and
                !the source field associated with each equation set
                DO eqnset_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(eqnset_idx)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                     IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
                            ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
                              & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                            NULLIFY(GEOMETRIC_VARIABLE)
                            NULLIFY(GEOMETRIC_PARAMETERS)
                            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,& 
                              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                             NULLIFY(ANALYTIC_VARIABLE)
                             NULLIFY(ANALYTIC_PARAMETERS)
                             IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                               CALL FIELD_VARIABLE_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,ANALYTIC_VARIABLE,ERR,ERROR,*999)
                               CALL FIELD_PARAMETER_SET_DATA_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                 & ANALYTIC_PARAMETERS,ERR,ERROR,*999)           
                             ENDIF
                             NULLIFY(MATERIALS_FIELD)
                             NULLIFY(MATERIALS_VARIABLE)
                             NULLIFY(MATERIALS_PARAMETERS)
                             IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
                               MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                               CALL FIELD_VARIABLE_GET(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,MATERIALS_VARIABLE,ERR,ERROR,*999)
                               CALL FIELD_PARAMETER_SET_DATA_GET(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                 & MATERIALS_PARAMETERS,ERR,ERROR,*999)           
                             ENDIF
                             EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)=CURRENT_TIME
!                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                              variable_type=DEPENDENT_FIELD%VARIABLES(2*eqnset_idx-1)%VARIABLE_TYPE
                              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
                                    & FIELD_NODE_BASED_INTERPOLATION) THEN
                                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                    IF(ASSOCIATED(DOMAIN)) THEN
                                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                            & FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                            !Loop over the local nodes excluding the ghosts.
                                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                              !!TODO \todo We should interpolate the geometric field here and the node position.
                                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                                !Default to version 1 of each node derivative
                                                local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP% &
                                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                                X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                              ENDDO !dim_idx
                                              !Loop over the derivatives
                                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                                ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                                GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                                  & GLOBAL_DERIVATIVE_INDEX
                                                CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET%SUBTYPE, &
                                                  & ANALYTIC_FUNCTION_TYPE,X,TANGENTS,NORMAL,CURRENT_TIME,variable_type, &
                                                  & GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS, &
                                                  & VALUE,ERR,ERROR,*999)
                                                !Default to version 1 of each node derivative
                                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                                BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                  & CONDITION_TYPES(local_ny)
                                                IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                                                 CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, & 
                                                   & variable_type,FIELD_VALUES_SET_TYPE,local_ny, & 
                                                   & VALUE,ERR,ERROR,*999)
                                                ENDIF

!                                              IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
!                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
!                                                  !If we are a boundary node then set the analytic value on the boundary
!                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
!                                                    & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
!                                                ENDIF
!                                              ENDIF
                                              ENDDO !deriv_idx
                                            ENDDO !node_idx
                                          ELSE
                                            CALL FLAG_ERROR("Boundary conditions variable is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ELSE
                                          CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !component_idx
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                              ELSE
                                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
                              ENDIF

!                              ENDDO !variable_idx
                             CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
                              & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                          ENDIF            
                        ELSE
                          CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        !CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
                  END IF                
!                 ELSE
!                   CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
!                 END IF  
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

            !>Set the source field to a specified analytical function
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
                IF(ASSOCIATED(SOURCE_FIELD)) THEN
                  GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                  IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                    NULLIFY(GEOMETRIC_VARIABLE)
                    CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
                    CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                      variable_type=FIELD_U_VARIABLE_TYPE
                      FIELD_VARIABLE=>SOURCE_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                          IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                            DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                            IF(ASSOCIATED(DOMAIN)) THEN
                              IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                  !Loop over the local nodes excluding the ghosts.
                                  DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                    !!TODO \todo We should interpolate the geometric field here and the node position.
                                    DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                      !Default to version 1 of each node derivative
                                      local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                      X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                    ENDDO !dim_idx
                                    !Loop over the derivatives
                                    DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                      SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM)
                                        SELECT CASE(eqnset_idx)
                                        CASE(1)
                                          VALUE_SOURCE=EXP(-1*CURRENT_TIME)*(-1*A1*(X(1)*X(1)+X(2)*X(2))-4*D1*A1+LAMBDA_12*(A1-A2)*&
                                          & (X(1)*X(1)+X(2)*X(2)))
                                        CASE(2)
                                          VALUE_SOURCE=EXP(-1*CURRENT_TIME)*(-1*A2*(X(1)*X(1)+X(2)*X(2))-4*D2*A2+LAMBDA_12*(A2-A1)*&
                                          & (X(1)*X(1)+X(2)*X(2)))
                                        END SELECT
                                      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM)
                                        SELECT CASE(eqnset_idx)
                                        CASE(1)
                                          VALUE_SOURCE=EXP(-1*CURRENT_TIME)*(-1*A1*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))-&
                                          & 6*D1*A1+LAMBDA_13*(A1-A3)*&
                                          & (X(1)*X(1)+X(2)*X(2)+X(3)*X(3))+LAMBDA_12*(A1-A2)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)))
                                        CASE(2)
                                          VALUE_SOURCE=EXP(-1*CURRENT_TIME)*(-1*A2*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))-&
                                          & 6*D2*A2+LAMBDA_12*(A2-A1)*&
                                          & (X(1)*X(1)+X(2)*X(2)+X(3)*X(3))+LAMBDA_23*(A2-A3)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)))
                                        CASE(3)
                                          VALUE_SOURCE=EXP(-1*CURRENT_TIME)*(-1*A3*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))-&
                                          & 6*D3*A3+LAMBDA_13*(A3-A1)*&
                                          & (X(1)*X(1)+X(2)*X(2)+X(3)*X(3))+LAMBDA_23*(A3-A2)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)))
                                        END SELECT
                                      CASE DEFAULT
                                        LOCAL_ERROR="The analytic function type of "// &
                                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))//&
                                          & " is invalid."
                                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                      END SELECT
                                      !Default to version 1 of each node derivative
                                      local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                        & FIELD_VALUES_SET_TYPE,local_ny,VALUE_SOURCE,ERR,ERROR,*999)
                                    ENDDO !deriv_idx
                                  ENDDO !node_idx
                                ELSE
                                  CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !component_idx
                        CALL FIELD_PARAMETER_SET_UPDATE_START(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                          & ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_FINISH(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                          & ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                  ENDIF            
                ELSE
                  CALL FLAG_ERROR("Equations set source field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
            ENDIF
            ENDDO !eqnset_idx
                ELSE
                  CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
                END IF  
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a multi-physics coupled diffusion equation type of a multi-physics problem class."
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
    CALL EXITS("MULTI_COMP_TRANSPORT_PRE_SOLVE_UPDATE_ANALYTIC_VALUES")
    RETURN
999 CALL ERRORS("MULTI_COMP_TRANSPORT_PRE_SOLVE_UPDATE_ANALYTIC_VALUES",ERR,ERROR)
    CALL EXITS("MULTI_COMP_TRANSPORT_PRE_SOLVE_UPDATE_ANALYTIC_VALUES")
    RETURN 1
  END SUBROUTINE MULTI_COMP_TRANSPORT_PRE_SOLVE_UPDATE_ANALYTIC_VALUES
  !   
  !================================================================================================================================
  !
  !>Sets up the multi-compartment transport problem post solve.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP), POINTER :: OUTPUT_DATA1(:),OUTPUT_DATA2(:),OUTPUT_DATA3(:),OUTPUT_DATA4(:),OUTPUT_DATA5(:)
    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE",ERR,ERROR,*999)
     NULLIFY(OUTPUT_DATA1)
     NULLIFY(OUTPUT_DATA2)
     NULLIFY(OUTPUT_DATA3)
     NULLIFY(OUTPUT_DATA4)
     NULLIFY(OUTPUT_DATA5)
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
!              CALL ADVECTION_DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                  
                   SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                   IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                    EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                    IF(ASSOCIATED(EQUATIONS)) THEN
                     EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                     IF(ASSOCIATED(EQUATIONS_SET)) THEN

!                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                         & FIELD_VALUES_SET_TYPE,OUTPUT_DATA1,ERR,ERROR,*999)
! 
!                        WRITE (*,*) OUTPUT_DATA1
!                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, & 
!                         & FIELD_VALUES_SET_TYPE,OUTPUT_DATA2,ERR,ERROR,*999)
! 
!                        WRITE (*,*) OUTPUT_DATA2
!                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, & 
!                         & FIELD_VALUES_SET_TYPE,OUTPUT_DATA3,ERR,ERROR,*999)
! 
!                        WRITE (*,*) OUTPUT_DATA3
!                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, & 
!                         & FIELD_VALUES_SET_TYPE,OUTPUT_DATA4,ERR,ERROR,*999)
! 
!                        WRITE (*,*) OUTPUT_DATA4
! 
!                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, & 
!                         & FIELD_VALUES_SET_TYPE,OUTPUT_DATA5,ERR,ERROR,*999)
! 
!                        WRITE (*,*) OUTPUT_DATA5

                     ENDIF
                    endif
                   ENDIF
              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
!              CALL DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a multi-compartment type of a multi physics problem class."
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

    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the diffuion-diffusion problem post solve output data.
  SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
                !CALL ADVECTION_DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !CALL DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a multi-compartment transport type of a multi physics problem class."
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
      
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 CALL ERRORS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    CALL EXITS("MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE_OUTPUT_DATA")
    RETURN 1
  END SUBROUTINE MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE_OUTPUT_DATA
      
  !   
  !================================================================================================================================
  !


END MODULE MULTI_COMPARTMENT_TRANSPORT_ROUTINES
