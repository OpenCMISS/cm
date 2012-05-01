!> \file
!> \author Chris Bradley
!> \brief This module handles all Poisson equations routines.
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

!>This module handles all Poisson equations routines.
MODULE POISSON_EQUATIONS_ROUTINES

  USE BASE_ROUTINES 
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES
! temporary input for vector-source
  USE FLUID_MECHANICS_IO_ROUTINES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC POISSON_EQUATION_EQUATIONS_SET_SETUP,POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & POISSON_EQUATION_ANALYTIC_CALCULATE, &
    & POISSON_EQUATION_EQUATIONS_SET_SUBTYPE_SET,POISSON_EQUATION_FINITE_ELEMENT_CALCULATE,POISSON_PRE_SOLVE, &
    & POISSON_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE,POISSON_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE, &
    & POISSON_EQUATION_PROBLEM_SUBTYPE_SET,POISSON_EQUATION_PROBLEM_SETUP,POISSON_POST_SOLVE
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE POISSON_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type,global_ny
    INTEGER(INTG) :: COUNT_DOF
    REAL(DP) :: VALUE,X(3)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR    
    INTEGER(INTG), ALLOCATABLE :: CONDITION(:)
    INTEGER(INTG), ALLOCATABLE :: DOF_NUMBER(:)
    REAL(DP), ALLOCATABLE :: VALUE_BC(:)
    
    CALL ENTERS("POISSON_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)

            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN

              !For testing Integrated Neumann Boundary Conditions
              IF((EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN).OR.&
                & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN_WITHOUT_SOURCE)) THEN
                ALLOCATE(CONDITION(12)) !For 2x2 quad
                ALLOCATE(DOF_NUMBER(12)) !For 2x2 quad
                ALLOCATE(VALUE_BC(12)) !For 2x2 quad
              ELSEIF((EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_1).OR.&
                & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_2)) THEN
                ALLOCATE(CONDITION(3)) !For 2x2 quad
                ALLOCATE(DOF_NUMBER(3)) !For 2x2 quad
                ALLOCATE(VALUE_BC(3)) !For 2x2 quad
              ELSEIF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_3) THEN
                ALLOCATE(CONDITION(6)) !For 2x2 quad
                ALLOCATE(DOF_NUMBER(6)) !For 2x2 quad
                ALLOCATE(VALUE_BC(6)) !For 2x2 quad
              ELSEIF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_6) THEN
                ALLOCATE(CONDITION(17)) !34 Originally
                ALLOCATE(DOF_NUMBER(17))
                ALLOCATE(VALUE_BC(17))
              ELSEIF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_7) THEN
                !ALLOCATE(CONDITION(578)) !Quad 8x8x8
                !ALLOCATE(DOF_NUMBER(578)) !Quad 8x8x8
                !ALLOCATE(VALUE_BC(578)) !Quad 8x8x8
                ALLOCATE(CONDITION(1250)) !Cubic 8x8x8
                ALLOCATE(DOF_NUMBER(1250)) !Cubic 8x8x8
                ALLOCATE(VALUE_BC(1250)) !Cubic 8x8x8
              ELSEIF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN_CUBIC) THEN
                ALLOCATE(CONDITION(46)) !Cubic 8x8
                ALLOCATE(DOF_NUMBER(46)) !Cubic 8x8
                ALLOCATE(VALUE_BC(46)) !Cubic 8x8
              ENDIF
              COUNT_DOF=0

              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES !U and deludeln
                variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS !u,v,w
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
                                local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                  & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                              ENDDO !dim_idx
                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1)
                                  !u=ln(4/(x+y+1^2))
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=LOG(4.0_DP/((X(1)+X(2)+1.0_DP)**2))
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !This is believed to be incorrect, should be: VALUE=-2.0_DP/(X(1)+X(2)+1.0_DP)
                                      VALUE=-2.0_DP*(X(1)+X(2))/(X(1)+X(2)+1.0_DP)
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  END SELECT
                                CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_2)
                                  CALL FLAG_ERROR("The analytic function type is not implemented yet.",ERR,ERROR,*999)
                                CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_3)
                                  CALL FLAG_ERROR("The analytic function type is not implemented yet.",ERR,ERROR,*999)
                                CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1)
                                  !u=ln(6/(x+y+z+1^2))
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=LOG(6.0_DP/((X(1)+X(2)+x(3)+1.0_DP)**2))
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=-3.0_DP/(X(1)+X(2)+X(3)+1.0_DP)
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE(TEST_CASE_NEUMANN)
                                  !u=x^2+3y^2+z^2 Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=X(1)**2+3*X(2)**2+X(3)**2  3D version - also needs a source of 10
                                      VALUE=X(1)**2+3*X(2)**2     !2D version - needs a source of 8
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=2*X(1)+6*X(2)
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 2x2 quad
                                      SELECT CASE(node_idx)
                                      CASE(2,3,4)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(6,11,16)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=2*X(1)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(10,15,20)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-2*X(1)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(22,23,24)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE DEFAULT
                                        VALUE=0
                                      END SELECT
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_NEUMANN_WITHOUT_SOURCE)
                                  !u=-3x^2+3y^2 Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=-3*X(1)**2+3*X(2)**2     !2D version - no source needed
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=-6*X(1)+6*X(2)
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 2x2 quad
                                      SELECT CASE(node_idx)
                                      CASE(2,3,4)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(6,11,16)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-6*X(1)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(10,15,20)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=6*X(1)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(22,23,24)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE DEFAULT
                                        VALUE=0
                                      END SELECT
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_DIRICHLET)
                                  !u=x^2+3y^2 Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !!VALUE=X(1)**2+3*X(2)**2+X(3)**2  3D version - also needs source of 10
                                      VALUE=X(1)**2+3*X(2)**2
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !Do nothing
                                      !VALUE=2*X(1)+6*X(2)+2*X(3)
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_1)
                                  !u=x^2+3y^2+z^2 Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=X(1)**2+3*X(2)**2+X(3)**2  3D version - also needs a source of 10
                                      VALUE=X(1)**2+3*X(2)**2     !2D version - needs a source of 8
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=2*X(1)+6*X(2)
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 2x2 quad
                                      SELECT CASE(node_idx)
                                      CASE(2,3,4)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE DEFAULT
                                        VALUE=0
                                      END SELECT
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_2)
                                  !u=x^2+3y^2+z^2 Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=X(1)**2+3*X(2)**2+X(3)**2  3D version - also needs a source of 10
                                      VALUE=X(1)**2+3*X(2)**2     !2D version - needs a source of 8
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=2*X(1)+6*X(2)
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 2x2 quad
                                      SELECT CASE(node_idx)
                                      CASE(22,23,24)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE DEFAULT
                                        VALUE=0
                                      END SELECT
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_3)
                                  !u=x^2+3y^2+z^2 Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=X(1)**2+3*X(2)**2+X(3)**2  3D version - also needs a source of 10
                                      VALUE=X(1)**2+3*X(2)**2     !2D version - needs a source of 8
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=2*X(1)+6*X(2)
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 2x2 quad
                                      SELECT CASE(node_idx)
                                      CASE(2,3,4)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(22,23,24)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE DEFAULT
                                        VALUE=0
                                      END SELECT
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_4)
                                  !u=x^2+3y^2+z^2 Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=X(1)**2+3*X(2)**2+X(3)**2  3D version - also needs a source of 10
                                      VALUE=X(1)**2+3*X(2)**2     !2D version - needs a source of 8
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=2*X(1)+6*X(2)
                                      !This is for testing Integrated Neumann Boundary Conditions
                                      !Do nothing
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_6)
                                  !u=y Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=X(2) !2D version - no source
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !VALUE=1
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 8x8 quad
                                      SELECT CASE(node_idx)
  !                                     CASE(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
  !                                       COUNT_DOF = COUNT_DOF+1
  !                                       VALUE=1
  !                                       VALUE_BC(COUNT_DOF)=VALUE
  !                                       DOF_NUMBER(COUNT_DOF)=global_ny
  !                                       CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                    !  CASE(239,256,255,272)
                                    !    COUNT_DOF = COUNT_DOF+1
                                    !    VALUE=5
                                    !    VALUE_BC(COUNT_DOF)=VALUE
                                    !    DOF_NUMBER(COUNT_DOF)=global_ny
                                    !    CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_INTEGRATED
                                      CASE(273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-1
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE DEFAULT
                                        VALUE=0
                                      END SELECT
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_7)
                                  !u=y Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=X(2)    !3D version - no source
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 8x8x8 quad
                                      IF((X(1)>=0.AND.X(1)<=1.0_dp).AND.(X(3)>=0.AND.X(3)<=1.0_dp)) THEN
                                      !IF((X(1)>=0.AND.X(1)<=0.5_dp).AND.(X(3)>=0.AND.X(3)<=0.5_dp)) THEN
                                        IF(X(2)==1.0_dp) THEN
                                          COUNT_DOF = COUNT_DOF+1
                                          VALUE=-1
                                          VALUE_BC(COUNT_DOF)=VALUE
                                          DOF_NUMBER(COUNT_DOF)=global_ny
                                          CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                        ELSEIF(X(2)==0.0_dp) THEN
                                          COUNT_DOF = COUNT_DOF+1
                                          VALUE=1
                                          VALUE_BC(COUNT_DOF)=VALUE
                                          DOF_NUMBER(COUNT_DOF)=global_ny
                                          CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                        ENDIF
                                      ENDIF
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(TEST_CASE_NEUMANN_CUBIC)
                                  !u=y Test case example used
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=X(1)**2+3*X(2)**2     !2D version - needs a source of 8
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !This is for testing Integrated Neumann Boundary Conditions

                                      !---------------This is for 8x8 quad
                                      SELECT CASE(node_idx)
                                      CASE(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE(602,603,604,605,606,607,608,609,610,611,612, &
                                         & 613,614,615,616,617,618,619,620,621,622,623,624)
                                        COUNT_DOF = COUNT_DOF+1
                                        VALUE=-6*X(2)
                                        VALUE_BC(COUNT_DOF)=VALUE
                                        DOF_NUMBER(COUNT_DOF)=global_ny
                                        CONDITION(COUNT_DOF)=BOUNDARY_CONDITION_NEUMANN_POINT
                                      CASE DEFAULT
                                        VALUE=0
                                      END SELECT
                                      !---------------

                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT

                                CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_2)
                                  CALL FLAG_ERROR("The analytic function type is not implemented yet.",ERR,ERROR,*999)
                                CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_3)
                                  CALL FLAG_ERROR("The analytic function type is not implemented yet.",ERR,ERROR,*999)
                                CASE(EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1,EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2)
  !\todo: This test case has been set up for the Pressure Poisson equation - needs to be formulated in a more general way
                                  !u=ln(4/(x+y+1^2))
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=SIN(2.0_DP*PI*X(1)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)*SIN(2.0_DP*PI*X(3)/10.0_DP)
                                      !  VALUE=2*X(1)*X(1)+2*X(2)
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
  !                                     do nothing, tbd
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE DEFAULT
                                  LOCAL_ERROR="The analytic function type of "// &
                                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                    & " is invalid."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                END SELECT
                                    !Default to version 1 of each node derivative
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)

                                global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)

                                SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                CASE(TEST_CASE_NEUMANN,TEST_CASE_NEUMANN_WITHOUT_SOURCE)
                                  SELECT CASE(node_idx)
                                  CASE(1,5,21,25) !For 2x2 Quad
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ENDIF
                                  CASE DEFAULT
                                  END SELECT
                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_1,TEST_CASE_MIXED_NEUMANN_DIRICHLET_4)
                                  SELECT CASE(node_idx)
                                  CASE(1,6,11,16,5,10,15,20,21,22,23,24,25) !For 2x2 Quad
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ENDIF
                                  CASE DEFAULT
                                  END SELECT
                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_2)
                                  SELECT CASE(node_idx)
                                  CASE(1,2,3,4,5,6,11,16,10,15,20,21,25) !For 2x2 Quad
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ENDIF
                                  CASE DEFAULT
                                  END SELECT
                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_3)
                                  SELECT CASE(node_idx)
                                  CASE(1,6,11,16,21,5,10,15,20,25) !For 2x2 Quad
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ENDIF
                                  CASE DEFAULT
                                  END SELECT
                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_6)
                                  SELECT CASE(node_idx)
                                  CASE(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                      VALUE=0
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ENDIF
                                  CASE DEFAULT
                                  END SELECT
                                CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_7)

  !                                 IF((X(1)>=0.AND.X(1)<=1.0_dp).AND.(X(3)>=0.AND.X(3)<=1.0_dp)) THEN
  !                                   IF(X(2)==0.0_dp) THEN
  !                                     IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
  !                                       CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type,local_ny, &
  !                                         & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
  !                                     ENDIF
  !                                   ENDIF
  !                                 ENDIF

                                  SELECT CASE(node_idx)
                                  CASE(301)
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ENDIF
                                  CASE DEFAULT
                                  END SELECT
                                CASE(TEST_CASE_NEUMANN_CUBIC)
                                  SELECT CASE(node_idx)
                                  CASE(1,26,51,76,101,126,151,176,201,226,251,276,301,326,351,376,401, &
                                       & 426,451,476,501,526,551,576,601, &
                                       & 25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400, &
                                       & 425,450,475,500,525,550,575,600,625) !For 8x8 Cubic
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ENDIF
                                  CASE DEFAULT
                                  END SELECT
                                CASE DEFAULT
                                  IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                      & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                  ENDIF
                                END SELECT

                                IF(variable_type==FIELD_DELUDELN_VARIABLE_TYPE.and.node_idx/=1) THEN
                                  IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                    !If we are a boundary node then set the analytic value on the boundary
                                    !CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type,local_ny, &
                                    !  & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    !Do nothing at present
                                  ENDIF
                                ENDIF

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
                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
                ENDIF

                !If FIELD_DELUDELN_VARIABLE_TYPE and TEST_CASE then Neumann
                IF((EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN_WITHOUT_SOURCE.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_1.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_2.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_3.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_6.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_7.OR. &
                  & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN_CUBIC))))))) &
                  & .AND.variable_type==FIELD_DELUDELN_VARIABLE_TYPE) THEN
                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,DOF_NUMBER,CONDITION,VALUE_BC,ERR,ERROR,*999)
                ENDIF

              ENDDO !variable_idx

              IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN.OR. &
                & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN_WITHOUT_SOURCE.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_1.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_2.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_3.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_6.OR. &
                  & (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_MIXED_NEUMANN_DIRICHLET_7.OR. &
                  & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==TEST_CASE_NEUMANN_CUBIC))))))) THEN
                DEALLOCATE(CONDITION)
                DEALLOCATE(DOF_NUMBER)
                DEALLOCATE(VALUE_BC)
              ENDIF

              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)

            ELSE
              CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("POISSON_EQUATION_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_ANALYTIC_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Sets up the Poisson equation type of a classical field equations set class.
  SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Poisson equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_EQUATIONS_SET_LINEAR_SOURCE_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_EQUATIONS_SET_LINEAR_SOURCE_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
        & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_EQUATIONS_SET_PRESSURE_POISSON_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_EQUATIONS_SET_NONLINEAR_SOURCE_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_EQUATIONS_SET_NONLINEAR_SOURCE_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Poisson equation type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Poisson equation type of an classical field equations set class.
  SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE)        
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
      CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, & 
        & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE, &
        & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE)                
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
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)        
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
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)        
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
          & " is not valid for a Poisson equation type of an classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a Poisson equation type of a classical field equations set class.
  SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_FIELD_USER_NUMBER
    TYPE(FIELD_TYPE), POINTER :: DUMMY_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DUMMY_FIELD_USER_NUMBER=0
      NULLIFY(DUMMY_FIELD)
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE
      CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE
      CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE
      CASE(EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE
      CASE(EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE
      CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_POISSON_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Poisson equation type of a classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_SUBTYPE_SET
  !
  !================================================================================================================================
  !

  !>Sets up the standard Poisson equation for Pressure Poisson Equation (PPE).
  SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_PRESSURE_POISSON_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS,GEOMETRIC_MESH_COMPONENT,SOURCE_FIELD_NUMBER_OF_COMPONENTS,I, &
      & SOURCE_FIELD_NUMBER_OF_VARIABLES
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: analytic_field,DEPENDENT_FIELD,geometric_field
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_EQUATIONS_SET_PRESSURE_POISSON_SETUP",ERR,ERROR,*999)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE.OR. & 
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE.OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a velocity source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              !These 3 parameter sets will contain the original velocity field with x,y,z components
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_INPUT_VEL1_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_INPUT_VEL2_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_INPUT_VEL3_SET_TYPE,ERR,ERROR,*999)
              !This parameter set will contain the inside/outside labelling information for the PPE approach
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_INPUT_LABEL_SET_TYPE,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a velocity source Poisson equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                  !Constant source. Materials field components are 1 for each dimension and 1 for the constant source
                  !i.e., k and c in div(k.grad(u(x)))=c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  !Linear source. Materials field components are 1 for each dimension and 2 for the linear source
                  !i.e., k and a and c in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2
                ENDIF
                !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                !Default the source materials components to the first component geometric interpolation with constant interpolation
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !components_idx
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+2, &
                    & ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
            
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                  !Constant source. Materials field components are 1 for each dimension and 1 for the constant source
                  !i.e., k and c in div(k.grad(u(x)))=c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  !Linear source. Materials field components are 1 for each dimension and 2 for the linear source
                  !i.e., k and a and c in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                !Now set the source values to 1.0
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a velocity source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
              !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Create the auto created source field
                !start field creation with name 'SOURCE_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                !label the field
                CALL FIELD_LABEL_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,"Source Field",ERR,ERROR, & 
                  & *999)
                !define new created field to be source
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                      & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,EQUATIONS_SET% & 
                  & GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                !set number of variables to 1 (1 for U)
                SOURCE_FIELD_NUMBER_OF_VARIABLES=1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                  & SOURCE_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                  & (/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension
                SOURCE_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,SOURCE_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                DO I=1,SOURCE_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                END DO
                
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,SOURCE_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                  !Other solutions not defined yet
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension and one for pressure
                SOURCE_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & SOURCE_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                    &"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
              !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                !These 2 parameter sets will contain the fitted hermite/lagrange velocity field
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)
!                 CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_INPUT_DATA3_SET_TYPE,ERR,ERROR,*999)
!                 CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_BOUNDARY_SET_TYPE,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Navier-Poisson fluid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for a Navier-Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
            !define an independent field for ALE information
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE)
          !do nothing!
          CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
          !do nothing!
          CASE(EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
             !Set start action
              CASE(EQUATIONS_SET_SETUP_START_ACTION)
                IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                  !Create the auto created independent field
                  !start field creation with name 'INDEPENDENT_FIELD'
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                    & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                  !start creation of a new field
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                  !label the field
                  CALL FIELD_LABEL_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",ERR,ERROR,*999)
                  !define new created field to be independent
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  !look for decomposition rule already defined
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & ERR,ERROR,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% & 
                    & GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                  !set number of variables to 1 (1 for U)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & 1,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & (/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  !calculate number of components with one component for each dimension
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  !Default to the geometric interpolation setup
                  DO I=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  END DO
                    SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                      !Specify fem solution method
                      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                        DO I=1,NUMBER_OF_DIMENSIONS
                          CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                          & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                        END DO
                        CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                          & ERR,ERROR,*999)
                        CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                          & ERR,ERROR,*999)
                        !Other solutions not defined yet
                      CASE DEFAULT
                        LOCAL_ERROR="The solution method of " &
                          & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT 
                  ELSE
                    !Check the user specified field
                    CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                      & ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                    !calculate number of components with one component for each dimension and one for pressure
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                    SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                        CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                          & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                        CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                          & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                          &"*",ERR,ERROR))//" is invalid."
                         CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ENDIF    
              !Specify finish action
              CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_MESH_VELOCITY_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_BOUNDARY_SET_TYPE,ERR,ERROR,*999)
                ENDIF
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard PPE fluid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for a PPE equation."
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
!                   CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1)
!                     !Set analtyic function type
!                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1
!                   CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2)
!                     !Set analtyic function type
!                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2
!                   CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3)
!                     !Set analtyic function type
!                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3
!                   CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1)
!                     !Set analtyic function type
!                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1
!                   CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2)
!                     !Set analtyic function type
!                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2
!                   CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_3)
!                     !Set analtyic function type
!                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_3
                  CASE(EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1)
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1
                  CASE(EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2)
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a velocity source Poisson equation."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a velocity source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Create the equations
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a velocity source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a velocity source Poisson equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not a velocity source Poisson equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_PRESSURE_POISSON_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_EQUATIONS_SET_PRESSURE_POISSON_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_PRESSURE_POISSON_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_PRESSURE_POISSON_SETUP

  !
  !================================================================================================================================
  !
  !>Update mesh velocity and move mesh for ALE PPE problem

  SUBROUTINE POISSON_PRE_SOLVE_UPDATE_PPE_MESH(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_PPE, SOLVER_LAPLACE !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_ALE_PPE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_ALE_PPE  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_ALE_PPE !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_ALE_PPE !<A pointer to the equations set
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS_ALE_PPE,GEOMETRIC_MESH_COMPONENT
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION,component_idx,deriv_idx,local_ny,node_idx,variable_idx,variable_type

    CALL ENTERS("POISSON_PRE_SOLVE_UPDATE_PPE_MESH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP%PARENT_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
      NULLIFY(SOLVER_LAPLACE)
      NULLIFY(SOLVER_ALE_PPE)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PARENT_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PARENT_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
              !Update mesh within the dynamic solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                !Get the independent field for the ALE PPE problem
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ALE_PPE,ERR,ERROR,*999)
                SOLVER_EQUATIONS_ALE_PPE=>SOLVER_ALE_PPE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_ALE_PPE)) THEN
                  SOLVER_MAPPING_ALE_PPE=>SOLVER_EQUATIONS_ALE_PPE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_ALE_PPE)) THEN
                    EQUATIONS_SET_ALE_PPE=>SOLVER_MAPPING_ALE_PPE%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_ALE_PPE)) THEN
                      INDEPENDENT_FIELD_ALE_PPE=>EQUATIONS_SET_ALE_PPE%INDEPENDENT%INDEPENDENT_FIELD
                    ELSE
                      CALL FLAG_ERROR("ALE PPE equations set is not associated.",ERR,ERROR,*999)
                    END IF
                    !Get the data
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_ALE_PPE%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS_ALE_PPE,ERR,ERROR,*999)
!\todo: Introduce flags set by the user (42/1 only for testings purpose)
                    !Copy input to PPE' independent field
                    INPUT_TYPE=42
                    INPUT_OPTION=1
                    NULLIFY(MESH_DISPLACEMENT_VALUES)
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_ALE_PPE%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, & 
                      & NUMBER_OF_DIMENSIONS_ALE_PPE,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%PARENT_LOOP%TIME_LOOP% &
                      & ITERATION_NUMBER,1.0_DP)
                    CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_PPE%INDEPENDENT%INDEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_PPE%INDEPENDENT%INDEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("ALE PPE solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("ALE PPE solver equations are not associated.",ERR,ERROR,*999)
                END IF
                 !Use calculated values to update mesh
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET_ALE_PPE%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
!                 CALL FIELD_PARAMETER_SET_DATA_GET(INDEPENDENT_FIELD_ALE_PPE,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                EQUATIONS=>SOLVER_MAPPING_ALE_PPE%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                  IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                    DO variable_idx=1,EQUATIONS_SET_ALE_PPE%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                      variable_type=EQUATIONS_SET_ALE_PPE%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                      FIELD_VARIABLE=>EQUATIONS_SET_ALE_PPE%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                          DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                          IF(ASSOCIATED(DOMAIN)) THEN
                            IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                !Loop over the local nodes excluding the ghosts.
                                DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                  DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(EQUATIONS_SET_ALE_PPE%GEOMETRY%GEOMETRIC_FIELD, & 
                                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                      & MESH_DISPLACEMENT_VALUES(local_ny),ERR,ERROR,*999)
                                  ENDDO !deriv_idx
                                ENDDO !node_idx
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDDO !component_idx
                      ENDIF
                    ENDDO !variable_idx
                  ELSE
                    CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_PPE%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_PPE%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                !Now use displacement values to calculate velocity values
                TIME_INCREMENT=CONTROL_LOOP%PARENT_LOOP%TIME_LOOP%TIME_INCREMENT
                ALPHA=1.0_DP/TIME_INCREMENT
                CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD_ALE_PPE,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,ERR,ERROR,*999)
              ELSE  
                CALL FLAG_ERROR("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a PPE equation fluid type of a fluid mechanics problem class."
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
    CALL EXITS("POISSON_PRE_SOLVE_UPDATE_PPE_MESH")
    RETURN
999 CALL ERRORS("POISSON_PRE_SOLVE_UPDATE_PPE_MESH",ERR,ERROR)
    CALL EXITS("POISSON_PRE_SOLVE_UPDATE_PPE_MESH")
    RETURN 1
  END SUBROUTINE POISSON_PRE_SOLVE_UPDATE_PPE_MESH

  !
  !================================================================================================================================
  !
  !>Update source for fitted PPE problem

  SUBROUTINE POISSON_PRE_SOLVE_UPDATE_PPE_SOURCE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_FITTED, SOLVER_PPE !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_FITTED,SOURCE_FIELD_PPE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_FITTED, SOLVER_EQUATIONS_PPE  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_FITTED, SOLVER_MAPPING_PPE !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_FITTED, EQUATIONS_SET_PPE !<A pointer to the equations set
!     TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
!     TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
!     TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
!     TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
!     TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: I,NUMBER_OF_DIMENSIONS_PPE,NUMBER_OF_DIMENSIONS_FITTED,GEOMETRIC_MESH_COMPONENT

    CALL ENTERS("POISSON_PRE_SOLVE_UPDATE_PPE_SOURCE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
      NULLIFY(SOLVER_PPE)
      NULLIFY(SOLVER_FITTED)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !do nothing for the fitting process
              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                !Get the dependent field for the three component Laplace problem
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_FITTED,ERR,ERROR,*999)
                SOLVER_EQUATIONS_FITTED=>SOLVER_FITTED%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_FITTED)) THEN
                  SOLVER_MAPPING_FITTED=>SOLVER_EQUATIONS_FITTED%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_FITTED)) THEN
                    EQUATIONS_SET_FITTED=>SOLVER_MAPPING_FITTED%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_FITTED)) THEN
                      DEPENDENT_FIELD_FITTED=>EQUATIONS_SET_FITTED%DEPENDENT%DEPENDENT_FIELD
                    ELSE
                      CALL FLAG_ERROR("Fitted equations set is not associated.",ERR,ERROR,*999)
                    END IF
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_FITTED%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS_FITTED,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Fitted solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Fitted solver equations are not associated.",ERR,ERROR,*999)
                END IF
                !Get the source field for the PPE problem
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_PPE,ERR,ERROR,*999)
                SOLVER_EQUATIONS_PPE=>SOLVER_PPE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_PPE)) THEN
                  SOLVER_MAPPING_PPE=>SOLVER_EQUATIONS_PPE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_PPE)) THEN
                    EQUATIONS_SET_PPE=>SOLVER_MAPPING_PPE%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_PPE)) THEN
                      SOURCE_FIELD_PPE=>EQUATIONS_SET_PPE%SOURCE%SOURCE_FIELD
                    ELSE
                      CALL FLAG_ERROR("PPE equations set is not associated.",ERR,ERROR,*999)
                    END IF
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_PPE%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS_PPE,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("PPE solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("PPE solver equations are not associated.",ERR,ERROR,*999)
                END IF
                !Copy result from FITTING to PPE's source field
                IF(NUMBER_OF_DIMENSIONS_PPE==NUMBER_OF_DIMENSIONS_FITTED) THEN
                  IF(CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER/=1) THEN
                    DO I=1,NUMBER_OF_DIMENSIONS_PPE
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Save old source field... ",ERR,ERROR,*999)
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(SOURCE_FIELD_PPE, & 
                          & FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,I,SOURCE_FIELD_PPE, & 
                          & FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,I,ERR,ERROR,*999)
                    END DO
                  ENDIF
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Update new source field... ",ERR,ERROR,*999)
                  DO I=1,NUMBER_OF_DIMENSIONS_PPE
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_FITTED, & 
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,SOURCE_FIELD_PPE, & 
                        & FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,I,ERR,ERROR,*999)
                  END DO
                  IF(CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER==1) THEN
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Old source field is new source field!... ",ERR,ERROR,*999)
                    DO I=1,NUMBER_OF_DIMENSIONS_PPE
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(SOURCE_FIELD_PPE, & 
                          & FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,I,SOURCE_FIELD_PPE, & 
                          & FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,I,ERR,ERROR,*999)
                    END DO
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Dimension of FITTED and PPE equations set is not consistent.",ERR,ERROR,*999)
                END IF
                !Use calculated values to update mesh
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET_PPE%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_PPE%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_PPE%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_PPE%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_PPE%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)
              ELSE  
                CALL FLAG_ERROR("Mesh update is not defined for non-dynamic problems.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a PPE equation fluid type of a fluid mechanics problem class."
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
    CALL EXITS("POISSON_PRE_SOLVE_UPDATE_PPE_SOURCE")
    RETURN
999 CALL ERRORS("POISSON_PRE_SOLVE_UPDATE_PPE_SOURCE",ERR,ERROR)
    CALL EXITS("POISSON_PRE_SOLVE_UPDATE_PPE_SOURCE")
    RETURN 1
  END SUBROUTINE POISSON_PRE_SOLVE_UPDATE_PPE_SOURCE

  !
  !================================================================================================================================
  !


  !>Sets up the standard Poisson equation for linear sources.
  SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_LINEAR_SOURCE_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS, GEOMETRIC_MESH_COMPONENT,SOURCE_FIELD_NUMBER_OF_COMPONENTS,I, &
      & SOURCE_FIELD_NUMBER_OF_VARIABLES
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: analytic_field,DEPENDENT_FIELD,geometric_field
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_EQUATION_SET_LINEAR_SOURCE_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE.OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                  !Constant source. Materials field components are 1 for each dimension and 1 for the constant source
                  !i.e., k and c in div(k.grad(u(x)))=c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  !Linear source. Materials field components are 1 for each dimension and 2 for the linear source
                  !i.e., k and a and c in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2
                ENDIF
                !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                !Default the source materials components to the first component geometric interpolation with constant interpolation
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !components_idx
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+2, &
                    & ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF

          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                  !Constant source. Materials field components are 1 for each dimension and 1 for the constant source
                  !i.e., k and c in div(k.grad(u(x)))=c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  !Linear source. Materials field components are 1 for each dimension and 2 for the linear source
                  !i.e., k and a and c in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                !Now set the source values to 1.0
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          CASE(EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
              SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
               !Set start action
                CASE(EQUATIONS_SET_SETUP_START_ACTION)
                  IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                    !Create the auto created source field
                    !start field creation with name 'SOURCE_FIELD'
                    CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                      & EQUATIONS_SET%SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                    !start creation of a new field
                    CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    !label the field
                    CALL FIELD_LABEL_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,"Source Field",ERR,ERROR, & 
                      & *999)
                    !define new created field to be source
                    CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                      & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                    !look for decomposition rule already defined
                    CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                      & ERR,ERROR,*999)
                    !apply decomposition rule found on new created field
                    CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                      & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                    !point new field to geometric field
                    CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,EQUATIONS_SET% & 
                      & GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                    !set number of variables to 1 (1 for U)
                    SOURCE_FIELD_NUMBER_OF_VARIABLES=1
                    CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                    & SOURCE_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                      & (/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                    CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                    !calculate number of components with one component for each dimension
                    SOURCE_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                    CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,SOURCE_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                      & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    !Default to the geometric interpolation setup
                    DO I=1,SOURCE_FIELD_NUMBER_OF_COMPONENTS
                      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                        & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    END DO
                    SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                      !Specify fem solution method
                      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                        DO I=1,SOURCE_FIELD_NUMBER_OF_COMPONENTS
                          CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                          & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                        END DO
                        CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                          & ERR,ERROR,*999)
                        CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,GEOMETRIC_SCALING_TYPE, &
                          & ERR,ERROR,*999)
                        !Other solutions not defined yet
                      CASE DEFAULT
                        LOCAL_ERROR="The solution method of " &
                          & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT 
                  ELSE
                    !Check the user specified field
                    CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                      & ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                    !calculate number of components with one component for each dimension and one for pressure
                    SOURCE_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                      & SOURCE_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                    SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                        CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                          & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                        CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                          & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                          &"*",ERR,ERROR))//" is invalid."
                         CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ENDIF    
                !Specify finish action
                CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                  IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                    CALL FIELD_CREATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                  ENDIF
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a linear Poisson subtype"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                & " is invalid for a linear Poisson equation."
               CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 2D
                      IF(NUMBER_OF_DIMENSIONS/=2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analtyic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an exponential source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1)
                    !Check that we have an exponential source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D
                      IF(NUMBER_OF_DIMENSIONS/=3) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be 3 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an exponential source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_NEUMANN)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS>2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_NEUMANN
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_NEUMANN_WITHOUT_SOURCE)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS>2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_NEUMANN_WITHOUT_SOURCE
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_DIRICHLET)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS>2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_DIRICHLET
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_1)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS>2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_MIXED_NEUMANN_DIRICHLET_1
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_2)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS>2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_MIXED_NEUMANN_DIRICHLET_2
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_3)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS>2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_MIXED_NEUMANN_DIRICHLET_3
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_4)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS>2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_MIXED_NEUMANN_DIRICHLET_4
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_6)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS<2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_MIXED_NEUMANN_DIRICHLET_6
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_MIXED_NEUMANN_DIRICHLET_7)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D or 2D
                      IF(NUMBER_OF_DIMENSIONS<3) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 3 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_MIXED_NEUMANN_DIRICHLET_7
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(TEST_CASE_NEUMANN_CUBIC)
                    !Check that we have an constant source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 2D
                      IF(NUMBER_OF_DIMENSIONS<2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be at least 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=TEST_CASE_NEUMANN_CUBIC
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an constant source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a linear source Poisson equation."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                ENDIF
             ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Create the equations
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear source Poisson equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not a linear source Poisson equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_LINEAR_SOURCE_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_EQUATIONS_SET_LINEAR_SOURCE_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_LINEAR_SOURCE_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_LINEAR_SOURCE_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the standard Poisson equation for nonlinear sources.
  SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_NONLINEAR_SOURCE_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: analytic_field,DEPENDENT_FIELD,geometric_field
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_EQUATION_SET_NONLINEAR_SOURCE_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE.OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL POISSON_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE) THEN
                  !Quadratic source. Materials field components are 1 for each dimension and 3 for the quadratic source
                  !i.e., k and a, b and c in div(k.grad(u(x)))=a(x)+b(x)u(x)+c(x)u^2(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ELSE
                  !Exponential source. Matierals field components are 1 for each dimension and 3 for the exponential source
                  !i.e., k, a, b and c in div(k.grad(u(x)))=a(x)+b(x)e^[c(x)u(x)]
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ENDIF
                !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                !Default the source materials components to the first component geometric interpolation with constant interpolation
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999) 
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !components_idx
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+3, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+3, &
                    & ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE) THEN
                  !Quadratic source. Materials field components are 1 for each dimension and 3 for the quadratic source
                  !i.e., k and a, b and c in div(k.grad(u(x)))=a(x)+b(x)u(x)+c(x)u^2(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ELSE
                  !Exponential source. Matierals field components are 1 for each dimension and 3 for the exponential source
                  !i.e., k, a, b and c in div(k.grad(u(x)))=a(x)+b(x)e^[c(x)u(x)]
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                !Now set the source values to 1.0
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
             IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1)
                    !Check that we have an exponential source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 2D
                      IF(NUMBER_OF_DIMENSIONS/=2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an exponential source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1)
                    !Check that we have an exponential source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE) THEN
                      !Check that we are in 3D
                      IF(NUMBER_OF_DIMENSIONS/=3) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be 3 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analytic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an exponential source Poisson equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a nonlinear source Poisson equation."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                ENDIF
             ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Start the creation of the equations
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              ! Use the analytic Jacobian calculation
              CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                & ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,MATRIX_BLOCK_STORAGE_TYPE, &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a nonlinear source Poisson equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not a nonlinear source Poisson equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_NONLINEAR_SOURCE_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_EQUATIONS_SET_NONLINEAR_SOURCE_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_EQUATIONS_SET_NONLINEAR_SOURCE_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_EQUATIONS_SET_NONLINEAR_SOURCE_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the Poisson problem.
  SUBROUTINE POISSON_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Poisson equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_PROBLEM_LINEAR_SOURCE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
        & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE,PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_PROBLEM_PRESSURE_POISSON_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
        CALL POISSON_EQUATION_PROBLEM_NONLINEAR_SOURCE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Poisson equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_PROBLEM_SETUP
  
   !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Poisson equation finite element equations set.
  SUBROUTINE POISSON_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,mi,ms,nh,nhs,ni,ns,I,J,K,L,H,element_node_identity !,k1sum,k2sum
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3),SUM2,DXI_DX(3,3),DELTA_T,DXI_DX_DX(3,3),PHINS
    REAL(DP) :: U_VALUE(3),U_DERIV(3,3),RHO_PARAM,MU_PARAM,U_OLD(3),U_SECOND(3,3,3),X(3),B(3),P_DERIV(3),W_VALUE(3)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS,SOURCE_BASIS,INDEPENDENT_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,SOURCE_FIELD,MATERIALS_FIELD,INDEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    LOGICAL :: INSIDE,BETWEEN
    REAL(DP), POINTER :: INPUT_LABEL(:)
    REAL(DP) :: DIFF_COEFF1,DIFF_COEFF2
    
    CALL ENTERS("POISSON_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems.
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns ???
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
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
                          SUM=SUM+PGMSI(mi)*PGNSI(ni)*EQUATIONS%INTERPOLATION% &
                            & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%GU(mi,ni)
                        ENDDO !ni
                      ENDDO !mi
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                  SUM=0.0_DP
                  !Additional terms for analytic solution
                  IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN

                    SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                    CASE(TEST_CASE_NEUMANN,TEST_CASE_DIRICHLET,TEST_CASE_MIXED_NEUMANN_DIRICHLET_1,&
                      & TEST_CASE_MIXED_NEUMANN_DIRICHLET_2,TEST_CASE_MIXED_NEUMANN_DIRICHLET_3,&
                      & TEST_CASE_MIXED_NEUMANN_DIRICHLET_4,TEST_CASE_NEUMANN_CUBIC)
                      SUM=-8.0_DP !For 2D TEST_CASE_NEUMANN and TEST_CASE_DIRICHLET and TEST_CASE_MIXED_NEUMANN_DIRICHLET
                    CASE DEFAULT
                     SUM=0.0_DP
                    END SELECT

                  ELSE !Use materials field value
                    SUM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                  ENDIF
                  SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ & !This needs to changed to the source vector
                     & SUM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*RWG
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng

        CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????

          INSIDE=.TRUE.
          BETWEEN=.FALSE.

          SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          SOURCE_BASIS=>SOURCE_FIELD%DECOMPOSITION%DOMAIN(SOURCE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE) THEN
            INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD
            INDEPENDENT_BASIS=>INDEPENDENT_FIELD%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD% &
              & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_MESH_VELOCITY_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & INDEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          ENDIF

          DIFF_COEFF1=1.0_DP
          DIFF_COEFF2=1.0_DP

          !Determine inside outside nodes/elements
          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_INPUT_LABEL_SET_TYPE,INPUT_LABEL,ERR,ERROR,*999)
          DO I=1,GEOMETRIC_FIELD%decomposition%domain(1)%ptr%topology%elements%maximum_number_of_element_parameters
            element_node_identity=GEOMETRIC_FIELD%decomposition%domain(1)%ptr%topology% &
              & elements%elements(ELEMENT_NUMBER)%element_nodes(I)
            !this needs to be changed even at the right INPUT_LABEL
            IF(INPUT_LABEL(element_node_identity)<0.5_DP) THEN
              !labelling if !!!ELEMENT!!! is inside or outside
              !so even if only one node is outside node, the whole element is outside element.
              INSIDE=.FALSE.
            ENDIF
          ENDDO
          DO I=1,GEOMETRIC_FIELD%decomposition%domain(1)%ptr%topology%elements%maximum_number_of_element_parameters
            element_node_identity=GEOMETRIC_FIELD%decomposition%domain(1)%ptr%topology% &
              & elements%elements(ELEMENT_NUMBER)%element_nodes(I)
            !this needs to be changed even at the right INPUT_LABEL
            IF(.NOT.INSIDE) THEN
              IF(INPUT_LABEL(element_node_identity)>0.5_DP) THEN
                !labelling if !!!ELEMENT!!! is inside or outside
                !so even if only one node is outside node, the whole element is outside element.
                BETWEEN=.TRUE.
              ENDIF
            ENDIF
          ENDDO

          IF(INSIDE) THEN
            DIFF_COEFF1=1.0_DP
            DIFF_COEFF2=1.0_DP
          ELSE IF(BETWEEN) THEN
!          set to "zero" if "between" nodes should be outside!   
!             DIFF_COEFF1=0.0_DP
            DIFF_COEFF1=1.0_DP/1000.0_DP
            DIFF_COEFF2=0.0_DP
          ELSE
            DIFF_COEFF1=0.0_DP
            DIFF_COEFF2=0.0_DP
          ENDIF
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              W_VALUE(1)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              W_VALUE(2)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
              IF(INDEPENDENT_BASIS%NUMBER_OF_XI==3) THEN
                W_VALUE(3)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
              END IF 
            ELSE
              W_VALUE=0.0_DP
            END IF

! ! !               W_VALUE=0.0_DP

            U_VALUE=0.0_DP
            U_OLD=0.0_DP
            DXI_DX=0.0_DP
            DXI_DX_DX=0.0_DP
            DELTA_T=1.0_DP !need to be dependent on example file
            U_DERIV=0.0_DP
            U_SECOND=0.0_DP
            IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_INPUT_DATA1_SET_TYPE,ELEMENT_NUMBER, & 
                & EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              U_VALUE(1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              U_VALUE(2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
              U_DERIV(1,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S1)
              U_DERIV(1,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S2)
              U_DERIV(2,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S1)
              U_DERIV(2,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S2)
              U_SECOND(1,1,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S1_S1)
              U_SECOND(1,1,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S1_S2)
              U_SECOND(1,2,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S1_S2)
              U_SECOND(1,2,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S2_S2)
              U_SECOND(2,1,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S1_S1)
              U_SECOND(2,1,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S1_S2)
              U_SECOND(2,2,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S1_S2)
              U_SECOND(2,2,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S2_S2)
              IF(DEPENDENT_BASIS%NUMBER_OF_XI==3) THEN
                U_VALUE(3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
                U_DERIV(1,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S3)
                U_DERIV(2,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S3) 
                U_DERIV(3,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S1)
                U_DERIV(3,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S2)
                U_DERIV(3,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S3)
                U_SECOND(1,1,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S1_S3)
                U_SECOND(1,2,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S2_S3)
                U_SECOND(1,3,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S1_S3)
                U_SECOND(1,3,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S2_S3)
                U_SECOND(1,3,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,PART_DERIV_S3_S3)
                U_SECOND(2,1,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S1_S3)
                U_SECOND(2,2,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S2_S3)
                U_SECOND(2,3,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S1_S3)
                U_SECOND(2,3,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S2_S3)
                U_SECOND(2,3,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,PART_DERIV_S3_S3)
                U_SECOND(3,1,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S1_S1)
                U_SECOND(3,1,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S1_S2)
                U_SECOND(3,1,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S1_S3)
                U_SECOND(3,2,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S1_S2)
                U_SECOND(3,2,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S2_S2)
                U_SECOND(3,2,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S2_S3)
                U_SECOND(3,3,1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S1_S3)
                U_SECOND(3,3,2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S2_S3)
                U_SECOND(3,3,3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,PART_DERIV_S3_S3)
              ENDIF
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_INPUT_DATA2_SET_TYPE,ELEMENT_NUMBER, & 
                & EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              U_OLD(1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              U_OLD(2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
              P_DERIV(1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S1)
              P_DERIV(2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S2)
              IF(DEPENDENT_BASIS%NUMBER_OF_XI==3) THEN
                U_OLD(3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
                P_DERIV(3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S3)
              ENDIF
            ENDIF
            !Define RHO_PARAM, density=1
            MU_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
            RHO_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns ???
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      PHINS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          DXI_DX(mi,ni)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & DXI_DX(mi,ni)
                        END DO
                      ENDDO !ni
                      SUM=0.0_DP

                      DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          SUM=SUM+PGMSI(mi)*PGNSI(ni)*EQUATIONS%INTERPOLATION% &
                            & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%GU(mi,ni)*DIFF_COEFF1
                        ENDDO !ni
                      ENDDO !mi
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE) THEN
                  SUM=0.0_DP
                  SUM2=0.0_DP
                  B=0.0_DP
                  IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                    DO K=1,3
                      !this is the temporal term (should be zero for static)
                      IF(INSIDE) THEN
                        B(K)=-RHO_PARAM*((U_VALUE(K)-U_OLD(K))/DELTA_T)
! !                       this is the first viscous term
!                       DO K=1,3
!                         DO L=1,3
!                           B(I)=B(I)+U_DERIV(I,L)*DXI_DX_DX(L,K)
!                         ENDDO  
!                       ENDDO  
!                       this is the second viscous term
                        DO L=1,3
                          DO I=1,3
                            ! this is the ALE term
                            B(K)=B(K)+RHO_PARAM*W_VALUE(I)*U_DERIV(K,L)*DXI_DX(L,I)
                            DO H=1,3
                              B(K)=B(K)+MU_PARAM*U_SECOND(K,L,H)*DXI_DX(L,I)*DXI_DX(H,I)
                            ENDDO  
                          ENDDO  
                        ENDDO  
                      ELSE
                        B(K)=0.0_DP
                      ENDIF
                      !eventually it gets combined to the RHS (source) term
                      DO J=1,3
                        SUM=SUM+B(K)*PGMSI(J)*DXI_DX(J,K)
                      ENDDO
                    ENDDO
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+SUM*RWG
                  ENDIF
                ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE.OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
                  SUM=0.0_DP
                  SUM2=0.0_DP
                  B=0.0_DP
                  IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                    DO K=1,3
                      IF(INSIDE) THEN
                        !this is the temporal term (should be zero for static)
                        B(K)=-RHO_PARAM*((U_VALUE(K)-U_OLD(K))/DELTA_T)
                        !this is the nonlinear term (should be zero for Stokes)
                        DO L=1,3
                          DO I=1,3 
                            !this is the nonlinear/ALE trem
                            B(K)=B(K)-RHO_PARAM*(U_VALUE(I)-W_VALUE(I))*U_DERIV(K,L)*DXI_DX(L,I)
                            DO H=1,3
                              !this is the viscous terms    
                              B(K)=B(K)+MU_PARAM*U_SECOND(K,L,H)*DXI_DX(L,I)*DXI_DX(H,I)
                            ENDDO
                          ENDDO
                        ENDDO
                      ELSE
                        B(K)=0.0_DP
                      ENDIF
                      !now bring the test function in
                      DO J=1,3
                        SUM=SUM+B(K)*PGMSI(J)*DXI_DX(J,K)*DIFF_COEFF2
                      ENDDO
                    ENDDO
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                      IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1) THEN
                        X(1) = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1)
                        X(2) = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,1)
                        X(3) = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,1)
                        SUM2=-4.0_DP*PI*PI/100.0*(3.0_DP*sin(2.0_DP*PI*X(1)/10.0_DP)*sin(2.0_DP*PI*X(2)/10.0)* &
                          & sin(2.0_DP*PI*X(3)/10.0_DP)-6.0_DP*RHO_PARAM*cos(2.0_DP*PI*X(1)/10.0_DP)**2+ &
                          & 8.0_DP*RHO_PARAM*cos(2.0_DP*PI*X(1)/10.0_DP)**2*cos(2.0_DP*PI*X(3)/10.0_DP)**2- &
                          & 2.0_DP*RHO_PARAM*cos(2.0_DP*PI*X(3)/10.0_DP)**2+2.0_DP*RHO_PARAM*cos(2.0_DP*PI*X(1)/10.0_DP)**2* &
                          & cos(2.0_DP*PI*X(2)/10.0_DP)**2+4.0_DP*RHO_PARAM*cos(2.0_DP*PI*X(2)/10.0_DP)**2- &
                          & 2.0_DP*RHO_PARAM*cos(2.0_DP*PI*X(2)/10.0_DP)**2*cos(2.0_DP*PI*X(3)/10.0_DP)**2)
                      ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2) THEN
                        X(1) = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1)
                        X(2) = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,1)
                        X(3) = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,1)
                        SUM2=-12.0_DP* sin(2.0_DP*PI*X(1)/10.0_DP)*PI*PI/100.0_DP*sin(2.0_DP*PI*X(2)/10.0_DP)* &
                          & sin(2.0_DP*PI*X(3)/10.0_DP)
                        SUM=0.0_DP
                      ELSE 
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+SUM*RWG &
                      & -SUM2*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*RWG
                  ENDIF
                ELSE
                  SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                IF(SOURCE_VECTOR%UPDATE_VECTOR) SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR% &
                  & ELEMENT_VECTOR%VECTOR(mhs)*EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)% &
                  & PTR%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
          CALL FLAG_ERROR("Can not calculate finite element stiffness matrices for a nonlinear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
          CALL FLAG_ERROR("Can not calculate finite element stiffness matrices for a nonlinear source.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Poisson equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("POISSON_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices and RHS for a Poisson equation finite element equations set.
  SUBROUTINE POISSON_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element Jacobian evaluation on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nh,nhs,ns
    REAL(DP) :: B_PARAM,C_PARAM,RWG,U_VALUE,VALUE
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
   
    CALL ENTERS("POISSON_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a Jacobian for a Poisson equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a Jacobian for a Poisson equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
          IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
            !Store all these in equations matrices/somewhere else?????
            DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
            GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
            MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            DEPENDENT_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
            GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Loop over gauss points
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
                & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              !Calculate RWG.
              RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
              !Find material parameters and u value at this Gauss point
              C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+3,NO_PART_DERIV)
              U_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              !Loop over field components
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  !Loop over element columns
                  DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      VALUE=-2.0_DP*C_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                        & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*U_VALUE
                      JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+VALUE*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh              
            ENDDO !ng
          ENDIF
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
          IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
            !Store all these in equations matrices/somewhere else?????
            DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
            GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
            MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            DEPENDENT_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
            GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Loop over gauss points
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
                & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              !Calculate RWG.
              RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
              !Find material parameter and u value at this Gauss point
              B_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
              C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+3,NO_PART_DERIV)
              U_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              !Loop over field components
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  !Loop over element columns
                  DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      VALUE=-B_PARAM*C_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                        & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*EXP(C_PARAM*U_VALUE)
                      JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+VALUE*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh              
            ENDDO !ng
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Poisson equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Poisson equation finite element equations set.
  SUBROUTINE POISSON_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nj,nh,nhs,ni,ns
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,K_PARAM,RWG,SUM1,SUM2,PGMJ(3),PGNJ(3),U_VALUE,WG
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("POISSON_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a Poisson equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a Poisson equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)                 
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          DEPENDENT_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems.
            WG=QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN*WG
            !Loop over field components
            IF(EQUATIONS_MATRIX%FIRST_ASSEMBLY) THEN
              IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                B_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                  & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
                mhs=0
                DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    nhs=0
                    !Loop over element columns
                    DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        SUM1=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,GEOMETRIC_BASIS%NUMBER_OF_XI
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(nj,NO_PART_DERIV)
                          SUM1=SUM1+K_PARAM*PGMJ(nj)*PGNJ(nj)
                        ENDDO !nj
                        SUM2=B_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)                        
                        EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & (SUM1+SUM2)*RWG
                      ENDDO !ns
                    ENDDO !nh
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(RHS_VECTOR%FIRST_ASSEMBLY) THEN
              IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                A_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                  & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+1,NO_PART_DERIV)
                mhs=0
                DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                       & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*A_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
              C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+3,NO_PART_DERIV)
              U_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(1,NO_PART_DERIV)
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)- &
                    & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*U_VALUE**2*RWG
                ENDDO !ms
              ENDDO !mh
            ENDIF
          ENDDO !ng
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)                 
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          DEPENDENT_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
           QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            IF(EQUATIONS_MATRIX%FIRST_ASSEMBLY) THEN
              IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                mhs=0
                DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    nhs=0
                    !Loop over element columns
                    DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        SUM1=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,GEOMETRIC_BASIS%NUMBER_OF_XI
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(nj,NO_PART_DERIV)
                          SUM1=SUM1+K_PARAM*PGMJ(nj)*PGNJ(nj)
                        ENDDO !nj
                        EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM1*RWG
                      ENDDO !ns
                    ENDDO !nh
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(RHS_VECTOR%FIRST_ASSEMBLY) THEN
              IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                A_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                  & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+1,NO_PART_DERIV)
                mhs=0
                DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                       & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*A_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
              B_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
              C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+3,NO_PART_DERIV)
              U_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(1,NO_PART_DERIV)
!!TODO: Handle floating point exceptions better
              IF((C_PARAM*U_VALUE)>20000.0_DP) THEN
                LOCAL_ERROR="The value of "//TRIM(NUMBER_TO_VSTRING(C_PARAM*U_VALUE,"*",ERR,ERROR))// &
                  & " is out of range for an exponential function."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)- &
                    & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*B_PARAM*EXP(C_PARAM*U_VALUE)*RWG
                ENDDO !ms
              ENDDO !mh
            ENDIF
          ENDDO !ng
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Poisson equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("POISSON_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a Poisson equation type.
  SUBROUTINE POISSON_EQUATION_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("POISSON_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_POISSON_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE     
      CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_POISSON_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE     
      CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_POISSON_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE
      CASE(PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_POISSON_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE
      CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_POISSON_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE
      CASE(PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_POISSON_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Poisson equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the linear source Poisson equations problem.
  SUBROUTINE POISSON_EQUATION_PROBLEM_LINEAR_SOURCE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT,SUB_LOOP,CONTROL_LOOP1,CONTROL_LOOP2
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
   
    CALL ENTERS("POISSON_EQUATION_PROBLEM_LINEAR_SOURCE_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
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
            !Set the solver to be a linear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
             !Start the linear solver creation
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear source Poisson equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a linear source Poisson equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_PROBLEM_LINEAR_SOURCE_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_PROBLEM_LINEAR_SOURCE_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_PROBLEM_LINEAR_SOURCE_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_PROBLEM_LINEAR_SOURCE_SETUP


  !
  !================================================================================================================================
  !

  !>Sets up the Pressure Poisson equations problem.
  SUBROUTINE POISSON_EQUATION_PROBLEM_PRESSURE_POISSON_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT,SUB_LOOP,CONTROL_LOOP1,CONTROL_LOOP2
    TYPE(SOLVER_TYPE), POINTER :: SOLVER,FITTING_SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,FITTING_SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
   
    CALL ENTERS("POISSON_EQUATION_PROBLEM_PRESSURE_POISSON_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(FITTING_SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(FITTING_SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE.OR.PROBLEM%SUBTYPE== &
        & PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE.OR.PROBLEM%SUBTYPE== &
        & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            NULLIFY(CONTROL_LOOP_ROOT)
            NULLIFY(CONTROL_LOOP1)
            NULLIFY(CONTROL_LOOP2)
            NULLIFY(SUB_LOOP)
            !Set up a simple control loop with 2 subloops
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(SUB_LOOP,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
            !Set the number of iterations for the while loop to 3 for now
! WRITE(*,*)'BE CAREFUL HERE ! MAKE IT MORE GENERAL'
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(SUB_LOOP,1,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          NULLIFY(SUB_LOOP)
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_ROOT,1,sub_LOOP,ERR,ERROR,*999)
          CALL CONTROL_LOOP_GET(sub_LOOP,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)

          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !Set the solver to be a dynamic solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
             !Start the linear solver creation
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(SUB_LOOP)
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_ROOT,1,sub_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_GET(sub_LOOP,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(SUB_LOOP)
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_ROOT,1,sub_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_GET(sub_LOOP,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear source Poisson equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE IF(PROBLEM%SUBTYPE==PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
          CASE(PROBLEM_SETUP_INITIAL_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Do nothing????
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Do nothing????
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a fitted PPE problem."
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
                  & " is invalid for a Pressure Poisson problem."
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
                !Set the first solver to be a linear solver for the velocity field fitting
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,FITTING_SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(FITTING_SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
                !Set solver defaults
                CALL SOLVER_LIBRARY_TYPE_SET(FITTING_SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                !Set the solver to be a quasi-static solvers
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Get the solvers
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                !Finish the solvers creation
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                 & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                 & " is invalid for a fitted PPE problem."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Get the control loop
                CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                !Get the solver
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,FITTING_SOLVER,ERR,ERROR,*999)
                !Create the solver equations
                CALL SOLVER_EQUATIONS_CREATE_START(FITTING_SOLVER,FITTING_SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(FITTING_SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(FITTING_SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(FITTING_SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                !Create the solver equations
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Get the control loop
                CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                !Get the solver equations
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,FITTING_SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(FITTING_SOLVER,FITTING_SOLVER_EQUATIONS,ERR,ERROR,*999)
                !Finish the solver equations creation
                CALL SOLVER_EQUATIONS_CREATE_FINISH(FITTING_SOLVER_EQUATIONS,ERR,ERROR,*999)             
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                !Finish the solver equations creation
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a fitted PPE problem."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a fitted PPE problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a linear source Poisson equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_PROBLEM_PRESSURE_POISSON_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_PROBLEM_PRESSURE_POISSON_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_PROBLEM_PRESSURE_POISSON_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_PROBLEM_PRESSURE_POISSON_SETUP

  
  !
  !================================================================================================================================
  !

  !>Sets up the nonlinear source Poisson equations problem.
  SUBROUTINE POISSON_EQUATION_PROBLEM_NONLINEAR_SOURCE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
    CALL ENTERS("POISSON_EQUATION_PROBLEM_NONLINEAR_SOURCE_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE) THEN
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
              & " is invalid for a nonlinear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation."
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
            !Set the solver to be nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear source Poisson equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a nonlinear source Poisson equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a nonlinear source Poisson equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("POISSON_EQUATION_PROBLEM_NONLINEAR_SOURCE_SETUP")
    RETURN
999 CALL ERRORS("POISSON_EQUATION_PROBLEM_NONLINEAR_SOURCE_SETUP",ERR,ERROR)
    CALL EXITS("POISSON_EQUATION_PROBLEM_NONLINEAR_SOURCE_SETUP")
    RETURN 1
  END SUBROUTINE POISSON_EQUATION_PROBLEM_NONLINEAR_SOURCE_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the Poisson problem post solve.
  SUBROUTINE POISSON_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solver
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("POISSON_POST_SOLVE",ERR,ERROR,*999)
    NULLIFY(SOLVER2)
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
              & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
              CALL POISSON_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
              !Post solve for the linear solver
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !do nothing for fitting problems
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Transferring fitting results to PPE... ",ERR,ERROR,*999)
              ELSE IF (SOLVER%GLOBAL_NUMBER==2) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"PPE post solve... ",ERR,ERROR,*999)
                CALL POISSON_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a Poisson type of a classical field problem class."
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

    CALL EXITS("POISSON_POST_SOLVE")
    RETURN
999 CALL ERRORS("POISSON_POST_SOLVE",ERR,ERROR)
    CALL EXITS("POISSON_POST_SOLVE")
    RETURN 1
  END SUBROUTINE POISSON_POST_SOLVE

  !
  !================================================================================================================================
  !


  !>Sets up the Poisson problem pre solve.
  SUBROUTINE POISSON_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solvers
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set

    CALL ENTERS("POISSON_PRE_SOLVE",ERR,ERROR,*999)
    NULLIFY(SOLVER2)
 
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
              IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1)THEN
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2) THEN
                          !do nothing
                        ELSE
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
                          !Update indpendent data fields
                          CALL POISSON_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"While loop... ",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
                        !Update indpendent data fields
                        CALL POISSON_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"While loop... ",ERR,ERROR,*999)
                      ENDIF 
                    ENDIF
                  ENDIF
                ENDIF
              ELSE
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"While loop... ",ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
              !Pre solve for the linear solver: fitting problem
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
                !Update independent data fields
                CALL POISSON_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Solving fitting problem... ",ERR,ERROR,*999)
              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Update source field ... ",ERR,ERROR,*999)
                  !First update source field for PPE calculation
                CALL POISSON_PRE_SOLVE_UPDATE_PPE_SOURCE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                  !Read in the lable information
!                 CALL POISSON_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Solving PPE problem... ",ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solver global number not associated for PPE problem.",ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
              IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1)THEN
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2) THEN
                          !do nothing
                        ELSE
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
                          !Update indpendent data fields
                          CALL POISSON_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"While loop... ",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
                        !Update indpendent data fields
                        CALL POISSON_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"While loop... ",ERR,ERROR,*999)
                      ENDIF 
                    ENDIF
                  ENDIF
                ENDIF
              !Update mesh
              CALL POISSON_PRE_SOLVE_UPDATE_PPE_MESH(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a Poisson type of a classical field problem class."
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

    CALL EXITS("POISSON_PRE_SOLVE")
    RETURN
999 CALL ERRORS("POISSON_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("POISSON_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE POISSON_PRE_SOLVE
  !
  !================================================================================================================================
  !
 

  !>Update boundary conditions for Poisson pre solve
  SUBROUTINE POISSON_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

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
! !     TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
! !     TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,CURRENT_LOOP_ITERATION
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_VEL_NEW_DATA(:),INPUT_VEL_OLD_DATA(:)
    REAL(DP), POINTER :: INPUT_VEL_LABEL_DATA(:),INPUT_VEL_U_DATA(:),INPUT_VEL_V_DATA(:),INPUT_VEL_W_DATA(:)
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control loop to solve.

    LOGICAL :: BOUNDARY_UPDATE

    BOUNDARY_UPDATE=.FALSE.

    CALL ENTERS("POISSON_PRE_SOLVE_UPDATE_INPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !only read in data if it is the fitting problem
!                 CONTROL_TIME_LOOP=>CONTROL_LOOP%PARENT_LOOP
                CONTROL_TIME_LOOP=>CONTROL_LOOP
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
!                 CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read input data xxxxxxxxxx... ",ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                     !Set pressure field to zero     
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                      CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                      !this is the current time step
!\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
                      INPUT_TYPE=1
                      INPUT_OPTION=1
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_VALUES_SET_TYPE,INPUT_VEL_NEW_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_NEW_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
                  END IF                
                ELSE
                  CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
                END IF 
              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                !this is the interior flag
!                 CONTROL_TIME_LOOP=>CONTROL_LOOP%PARENT_LOOP
                CONTROL_TIME_LOOP=>CONTROL_LOOP
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
!                 CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read input data xxxxxxxxxx... ",ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                     !Set pressure field to zero     
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                      CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                      !this is the current time step
!\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
                      INPUT_TYPE=1
                      INPUT_OPTION=3
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_LABEL_SET_TYPE,INPUT_VEL_LABEL_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_LABEL_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the reference U velocity
                      INPUT_TYPE=1
                      INPUT_OPTION=4
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_VEL1_SET_TYPE,INPUT_VEL_U_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_U_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the reference V velocity
                      INPUT_TYPE=1
                      INPUT_OPTION=5
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_VEL2_SET_TYPE,INPUT_VEL_V_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_V_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the reference W velocity
                      INPUT_TYPE=1
                      INPUT_OPTION=6
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_VEL3_SET_TYPE,INPUT_VEL_W_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_W_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
              & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
                CONTROL_TIME_LOOP=>CONTROL_LOOP%PARENT_LOOP
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read input data... ",ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
!                      !Set pressure field to zero     
!                       CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                         & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
!                       CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
!                         & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
!                       !Set pressure reference point in order to solve for absolute pressure
!                       IF(BOUNDARY_UPDATE) THEN
! !\todo: Let user determine when to change reference value
!                         IF(CURRENT_TIME<0.35_DP) THEN
!                           CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
!                             & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,10000.0_DP,ERR,ERROR,*999)
!                         ELSE
!                           CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
!                             & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
!                         ENDIF
!                           CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                             & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                           CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                             & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                       ENDIF
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                      CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                      !this is the current time step
!\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
                      INPUT_TYPE=1
                      INPUT_OPTION=1
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_DATA1_SET_TYPE,INPUT_VEL_NEW_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_NEW_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the previous time step
!\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
                      INPUT_TYPE=1
                      INPUT_OPTION=2
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_DATA2_SET_TYPE,INPUT_VEL_OLD_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_OLD_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the interior flag
                      INPUT_TYPE=1
                      INPUT_OPTION=3
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_LABEL_SET_TYPE,INPUT_VEL_LABEL_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_LABEL_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the reference U velocity
                      INPUT_TYPE=1
                      INPUT_OPTION=4
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_VEL1_SET_TYPE,INPUT_VEL_U_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_U_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the reference V velocity
                      INPUT_TYPE=1
                      INPUT_OPTION=5
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_VEL2_SET_TYPE,INPUT_VEL_V_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_V_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      !this is the reference W velocity
                      INPUT_TYPE=1
                      INPUT_OPTION=6
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_INPUT_VEL3_SET_TYPE,INPUT_VEL_W_DATA,ERR,ERROR,*999)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_W_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
                  END IF                
                ELSE
                  CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
                END IF 
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a Poisson type of a classical field problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
               & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
               & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
               & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
               & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)
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
    CALL EXITS("POISSON_PRE_SOLVE_UPDATE_INPUT_DATA")
    RETURN
999 CALL ERRORS("POISSON_PRE_SOLVE_UPDATE_INPUT_DATA",ERR,ERROR)
    CALL EXITS("POISSON_PRE_SOLVE_UPDATE_INPUT_DATA")
    RETURN 1
  END SUBROUTINE POISSON_PRE_SOLVE_UPDATE_INPUT_DATA

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE POISSON_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control loop to solve.
    LOGICAL :: EXPORT_FIELD
    TYPE(VARYING_STRING) :: METHOD!,FILE
    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    CALL ENTERS("POISSON_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!       write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
!               do nothing
            CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
              & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
              IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS)THEN
                CONTROL_TIME_LOOP=>CONTROL_LOOP%PARENT_LOOP
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    !Make sure the equations sets are up to date
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                      OUTPUT_ITERATION_NUMBER=CONTROL_TIME_LOOP%TIME_LOOP%OUTPUT_NUMBER
                      IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                        IF(CONTROL_TIME_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_TIME_LOOP%TIME_LOOP%STOP_TIME) THEN
                          IF(CURRENT_LOOP_ITERATION<10) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                          ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                          ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                          ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                            WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                          END IF
                          FILE=OUTPUT_FILE
  !          FILE="TRANSIENT_OUTPUT"
                          METHOD="FORTRAN"
                          EXPORT_FIELD=.TRUE.
                          IF(EXPORT_FIELD) THEN          
                            IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",ERR,ERROR,*999)
                              CALL FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER, &
                                & OUTPUT_FILE,ERR,ERROR,*999)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                            ENDIF
                          ENDIF 
                        ENDIF 
                      ENDIF
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF
            CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
              CONTROL_TIME_LOOP=>CONTROL_LOOP
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_TIME_LOOP%TIME_LOOP%OUTPUT_NUMBER
                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_TIME_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_TIME_LOOP%TIME_LOOP%STOP_TIME) THEN
                      IF(CURRENT_LOOP_ITERATION<10) THEN
                         WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                        END IF
                        FILE=OUTPUT_FILE
!          FILE="TRANSIENT_OUTPUT"
                        METHOD="FORTRAN"
                        EXPORT_FIELD=.TRUE.
                        IF(EXPORT_FIELD) THEN          
                          IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",ERR,ERROR,*999)
                            CALL FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER, &
                              & OUTPUT_FILE,ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                          ENDIF
                        ENDIF 
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a Poisson equation fluid type of a fluid mechanics problem class."
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
    CALL EXITS("POISSON_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 CALL ERRORS("POISSON_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    CALL EXITS("POISSON_POST_SOLVE_OUTPUT_DATA")
    RETURN 1
  END SUBROUTINE POISSON_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

END MODULE POISSON_EQUATIONS_ROUTINES
