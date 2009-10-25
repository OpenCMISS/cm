!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all linear elasticity routines.
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

!>This module handles all linear elasticity routines.
MODULE LINEAR_ELASTICITY_ROUTINES

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
  USE MATHS

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE,LINEAR_ELASTICITY_EQUATIONS_SET_SETUP, &
    & LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET,LINEAR_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET, &
    & LINEAR_ELASTICITY_PROBLEM_SUBTYPE_SET,LINEAR_ELASTICITY_PROBLEM_SETUP

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type,BC_X_counter
    REAL(DP) :: ANALYTIC_VALUE,BC_VALUE,X(3)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    LOGICAL :: SET_BC
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR    
    
    CALL ENTERS("LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR,*999)
    
    !!TODO: Use Geometric/Material Field values to prescribe values in analytic solution, currently hardcodded geometry & material properties
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
            !
            ! IDENTIFY BOUNDARY CONDITION NODES
            !
            BC_X_counter = 0
            DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
              variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
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
                            IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                              !!TODO \todo We should interpolate the geometric field here and the node position.
                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
                                X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                              ENDDO !dim_idx
                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES

                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !If we are a boundary node
                                  SET_BC = .FALSE.
                                  SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                  !
                                  ! TWO DIMENSIONAL LINEAR ELASTICITY
                                  !
                                  CASE(EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_TWO_DIM_PLANE_STRESS_1)
                                    SELECT CASE(component_idx)
                                    CASE(1) !u component
                                      SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                       SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          IF (X(1)==0.0_DP) THEN
                                            BC_X_counter = BC_X_counter + 1
                                          ENDIF
                                        END SELECT
                                      END SELECT
                                    CASE(2) !v component
                                      SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                       SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      END SELECT
                                    END SELECT

                                  !
                                  ! THREE DIMENSIONAL LINEAR ELASTICITY
                                  !
                                  CASE(EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_THREE_DIM_1)
                                    SELECT CASE(component_idx)
                                    CASE(1) !u component
                                      SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                       SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          IF (X(1)==0.0_DP) THEN
                                            BC_X_counter = BC_X_counter + 1
                                          ENDIF
                                        END SELECT
                                      END SELECT
                                    CASE(2) !v component
                                      SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                       SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      END SELECT
                                    CASE(3) !w component
                                      SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                       SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !pass
                                        END SELECT
                                      END SELECT
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The analytic function type of "// &
                                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                ENDIF
                              ENDDO !deriv_idx
                            ENDIF
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
              ELSE
                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !variable_idx

            !
            ! SET BOUNDARY CONDITIONS & ANALYTIC SOLUTION VALUES
            !
            NULLIFY(BOUNDARY_CONDITIONS)
            CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
              variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
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
                              local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
                              X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                            ENDDO !dim_idx
                            !Loop over the derivatives
                            DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                              SET_BC = .FALSE.
                              SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                              !
                              ! TWO DIMENSIONAL LINEAR ELASTICITY
                              !
                              CASE(EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_TWO_DIM_PLANE_STRESS_1)
                                SELECT CASE(component_idx)
                                CASE(1) !u component
                                !u=Sigmax*x/E
                                  SELECT CASE(variable_type)
                                  !!TODO set material parameters from material field
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=X(1)/10E3_DP
                                      IF (X(1)==0.0_DP) THEN
                                        SET_BC = .TRUE.
                                        BC_VALUE=0.0_DP
                                      ENDIF
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=1.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=1.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=1.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      IF (X(1)==0.0_DP) THEN
                                        ANALYTIC_VALUE=-100.0_DP/BC_X_counter
                                      ELSE!IF (X(1)==20.0_DP) THEN
                                        ANALYTIC_VALUE=0.0_DP
                                      ENDIF
                                      IF (X(1)==20.0_DP) THEN
                                        SET_BC = .TRUE.
                                        BC_VALUE=100.0_DP/BC_X_counter
                                      ENDIF
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=1.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=1.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=1.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE(2) !v component
                                !v=Sigmay*y/E
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=-(0.3_DP*X(2))/10E3_DP
                                      IF (X(2)==0.0_DP) THEN
                                        SET_BC = .TRUE.
                                        BC_VALUE=0.0_DP
                                      ENDIF
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE DEFAULT
                                  LOCAL_ERROR="The component index of "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                                    & " is invalid."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                END SELECT
                              !
                              ! THREE DIMENSIONAL LINEAR ELASTICITY
                              !
                              CASE(EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_THREE_DIM_1)
                                SELECT CASE(component_idx)
                                CASE(1) !u component
                                !u=Sigmax*x/E
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=X(1)*100.0_DP/(2500.0_DP*30E6_DP)
                                      IF (X(1)==0.0_DP) THEN
                                        SET_BC = .TRUE.
                                        BC_VALUE=0.0_DP
                                      ENDIF
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      IF (X(1)==0.0_DP) THEN
                                        ANALYTIC_VALUE=-100.0_DP/BC_X_counter
                                      ELSE
                                        ANALYTIC_VALUE=0.0_DP
                                      ENDIF
                                      IF (X(1)==100.0_DP) THEN
                                        SET_BC = .TRUE.
                                        BC_VALUE=100.0_DP/BC_X_counter
                                      ENDIF
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT                                
                                CASE(2) !v component
                                !v=Sigmax*x/E
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=-X(2)*0.25_DP*100.0_DP/(2500.0_DP*30E6_DP)
                                      IF (X(2)==0.0_DP) THEN
                                        SET_BC = .TRUE.
                                        BC_VALUE=0.0_DP
                                      ENDIF
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT                                
                                CASE(3) !w component
                                !w=Sigmax*x/E
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=-X(3)*0.25_DP*100.0_DP/(2500.0_DP*30E6_DP)
                                      IF (X(3)==0.0_DP) THEN
                                        SET_BC = .TRUE.
                                        BC_VALUE=0.0_DP
                                      ENDIF
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                    CASE(NO_GLOBAL_DERIV)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      ANALYTIC_VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE DEFAULT
                                  LOCAL_ERROR="The component index of "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                                    & " is invalid."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                END SELECT
                              CASE DEFAULT
                                LOCAL_ERROR="The analytic function type of "// &
                                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                  & " is invalid."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              END SELECT
                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,ANALYTIC_VALUE,ERR,ERROR,*999)
                              IF (SET_BC) THEN
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                                CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
                                  & BOUNDARY_CONDITION_FIXED,BC_VALUE,ERR,ERROR,*999)
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
            ENDDO !variable_idx
            CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
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
    
    CALL EXITS("LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a linear elasticity finite element equations set.
  SUBROUTINE LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element 
    !<calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    INTEGER(INTG) :: ng,ni,xi,mi,ns,nhs,ms,mhs,TOTAL_DEPENDENT_BASIS_EP,MESH_COMPONENT
    INTEGER(INTG) :: NUMBER_OF_DEPENDENT_COMPONENTS,NUMBER_OF_XI !,NUMBER_OF_GEOMETRIC_COMPONENTS
    INTEGER(INTG) :: OFF_DIAG_COMP(3),OFF_DIAG_DEP_VAR(2,2,3),DIAG_SUB_MAT_LOC(3),OFF_DIAG_SUB_MAT_LOC(2,3)
    INTEGER(INTG) :: DEPENDENT_BASES_EP(3) !,GEOMETRIC_BASES_EP(:)
    REAL(DP) :: JRWG,C(6,6),JRWG_DIAG_C(3,3),JRWG_OFF_DIAG_C(2,3)
    REAL(DP):: SF(64*3)

    !LOGICAL :: SAME_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,FIBRE_FIELD
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(BASIS_PTR_TYPE) :: DEPENDENT_BASES(3) !,GEOMETRIC_BASES(:)
    TYPE(QUADRATURE_SCHEME_PTR_TYPE) :: QUADRATURE_SCHEMES(3)
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: DEPENDENT_INTERPOLATION_PARAMETERS
    !TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FIBRE_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: MATERIALS_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: MATERIALS_INTERP_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERP_POINT_METRICS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE DPHI_DX_COMP_TYPE !A type to store DPHIDX for each mesh component
    REAL(DP) :: DPHIDX(64,3)
    END TYPE DPHI_DX_COMP_TYPE
    TYPE(DPHI_DX_COMP_TYPE) :: DPHIDX_COMP(3)

    CALL ENTERS("LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)
    !!Have a look at XPES40.f in the old CMISS code.
    !!Q - CPB: Need to think about anisotropic materials with fibre fields.
    !!Q - CPB: why store this DPHIDX(ns,xi) as opposed to just using it directly? A - to minimize operations - Otherwise it would be calculated many more times than
    !!         necessary within the loops below 
    !!Q - TPBG: Need to be able to use different Quadrature schemes with different bases? A - No use highest quadrature scheme for all directions
    !!TODO:: Check whether quadrature scheme being used is suffient to interpolate highest order basis function
    !!Q - TPBG: Need to be able to use different Interpolation for Geometric & Dependent field?
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
        GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
        !FIBRE_FIELD=>EQUATIONS%INTERPOLATION%FIBRE_FIELD
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        EQUATIONS_MATRIX=>EQUATIONS_MATRICES%LINEAR_MATRICES%MATRICES(1)%PTR
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        !Note that the dimension/number of components for FIELD_U_VARIABLE_TYPE has to be the same as the FIELD_DELUDELN_VARIABLE_TYPE
        !& Number of Geometric field components = number of dependent field component
        NUMBER_OF_DEPENDENT_COMPONENTS=GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(1)%ptr%NUMBER_OF_COMPONENTS
        !NUMBER_OF_GEOMETRIC_COMPONENTS=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(1)%ptr%NUMBER_OF_COMPONENTS
        !!TODO:: Use highest interpolation scheme's guass points. Warn if Gauss Points insufficient
        !Create an array of Bases with each component 
        DO ns=1,NUMBER_OF_DEPENDENT_COMPONENTS
          MESH_COMPONENT=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(1)%ptr%COMPONENTS(ns)%mesh_component_number
          DEPENDENT_BASES(ns)%PTR=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASES_EP(ns)=DEPENDENT_BASES(ns)%PTR%NUMBER_OF_ELEMENT_PARAMETERS
          QUADRATURE_SCHEMES(ns)%PTR=>DEPENDENT_BASES(ns)%PTR%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        ENDDO
        NUMBER_OF_XI = DEPENDENT_BASES(1)%PTR%NUMBER_OF_XI
        TOTAL_DEPENDENT_BASIS_EP = SUM(DEPENDENT_BASES_EP)
        !DO ns=1,NUMBER_OF_GEOMETRIC_COMPONENTS
        !  MESH_COMPONENT=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(1)%ptr%COMPONENTS(ns)%mesh_component_number
        !  GEOMETRIC_BASES(ns)%PTR=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR% &
        !    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        !  GEOMETRIC_BASES_EP(ns)=GEOMETRIC_BASES(ns)%PTR%NUMBER_OF_ELEMENT_PARAMETERS
        !ENDDO
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        !
        !ONE, TWO & THREE DIMENSIONAL LINEAR ELASTICITY
        !
        CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE,EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
          & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE,EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
          !//
          !Parameters for number of off diagonal stress/strain terms for a given number of xi directions and order of calculation for shear terms
          !These parameters do not change for 1D,2D,3D Linear Elasticity
          OFF_DIAG_COMP = (/0,1,3/)
          OFF_DIAG_DEP_VAR(1,1,:) = (/1,1,2/)
          OFF_DIAG_DEP_VAR(1,2,:) = (/2,3,3/)
          OFF_DIAG_DEP_VAR(2,1,:) = OFF_DIAG_DEP_VAR(1,2,:)
          OFF_DIAG_DEP_VAR(2,2,:) = OFF_DIAG_DEP_VAR(1,1,:)
          !//
          DIAG_SUB_MAT_LOC(:) = (/0,DEPENDENT_BASES_EP(1),SUM(DEPENDENT_BASES_EP(1:2))/)
          OFF_DIAG_SUB_MAT_LOC(1,:) = (/0,0,DEPENDENT_BASES_EP(1)/)
          OFF_DIAG_SUB_MAT_LOC(2,:) = (/DEPENDENT_BASES_EP(1),DIAG_SUB_MAT_LOC(3),DIAG_SUB_MAT_LOC(3)/)

          IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
            !Get the geometric, fibre & material field interpolation parameters
            GEOMETRIC_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS
            !FIBRE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS
            MATERIALS_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            !CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            !  & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & MATERIALS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            GEOMETRIC_INTERP_POINT_METRICS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS
            MATERIALS_INTERP_POINT=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT
            !Loop over gauss points & integrate upper triangular portion of Stiffness matrix
            DO ng=1,QUADRATURE_SCHEMES(1)%PTR%NUMBER_OF_GAUSS !Gauss point index
              !Interpolate geometric, fibre and material fields at gauss points & calculate geometric field metrics
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
              !CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              !  & FIBRE_INTERP_POINT,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,MATERIALS_INTERP_POINT,ERR,ERROR,*999)
              !!TODO:: Add option to only evaluate required metrics
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(NUMBER_OF_XI,GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
              !Calculate JRWG.
              JRWG=QUADRATURE_SCHEMES(1)%PTR%GAUSS_WEIGHTS(ng)*GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN
              DO xi=1,NUMBER_OF_XI
                DPHIDX_COMP(xi)%DPHIDX = 0.0_DP
                DO ns=1,DEPENDENT_BASES_EP(xi)
                  DO mi=1,NUMBER_OF_XI
                    DO ni=1,NUMBER_OF_XI
                      DPHIDX_COMP(xi)%DPHIDX(ns,mi) = DPHIDX_COMP(xi)%DPHIDX(ns,mi)+GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(mi,ni)* &
                        & QUADRATURE_SCHEMES(xi)%PTR%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                    ENDDO !ni
                  ENDDO !mi
                ENDDO !ns
              ENDDO !xi
              !CALL COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT, &
              !  EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT,DXDNU,ERR,ERROR,*999)
              !Create Linear Elasticity Tensor C
              CALL LINEAR_ELASTICITY_TENSOR(EQUATIONS_SET%SUBTYPE,MATERIALS_INTERP_POINT,C,ERR,ERROR,*999)
              !Store Elasticity Tensor diagonal & off diagonal stress coefficients
              JRWG_DIAG_C(3,:) = (/JRWG*C(4,4),JRWG*C(5,5),JRWG*C(3,3)/)
              JRWG_DIAG_C(2,:) = (/JRWG*C(6,6),JRWG*C(2,2),JRWG_DIAG_C(3,2)/)
              JRWG_DIAG_C(1,:) = (/JRWG*C(1,1),JRWG_DIAG_C(2,1),JRWG_DIAG_C(3,1)/)
              JRWG_OFF_DIAG_C(1,:) = (/JRWG*C(1,2),JRWG*C(1,3),JRWG*C(2,3)/)
              JRWG_OFF_DIAG_C(2,:) = (/JRWG*C(6,6),JRWG*C(4,4),JRWG*C(5,5)/)
              !Construct Element Matrix Diagonal Terms
              DO xi=1,NUMBER_OF_XI
                DO ns=1,DEPENDENT_BASES_EP(xi)
                  DO ms=ns,DEPENDENT_BASES_EP(xi)
                    EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(DIAG_SUB_MAT_LOC(xi)+ns,DIAG_SUB_MAT_LOC(xi)+ms) = &
                      & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(DIAG_SUB_MAT_LOC(xi)+ns,DIAG_SUB_MAT_LOC(xi)+ms) + &
                      & DOT_PRODUCT(DPHIDX_COMP(xi)%DPHIDX(ns,:)*DPHIDX_COMP(xi)%DPHIDX(ms,:),JRWG_DIAG_C(xi,1:NUMBER_OF_XI))
                  ENDDO !ms
                ENDDO !ns
              ENDDO !xi
              !Construct Element Matrix Off-Diagonal Terms
              DO xi=1,OFF_DIAG_COMP(NUMBER_OF_XI)
                DO ns=1,DEPENDENT_BASES_EP(OFF_DIAG_DEP_VAR(1,1,xi))
                  DO ms=1,DEPENDENT_BASES_EP(OFF_DIAG_DEP_VAR(1,2,xi))
                    EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(OFF_DIAG_SUB_MAT_LOC(1,xi)+ns,OFF_DIAG_SUB_MAT_LOC(2,xi)+ms) = &
                      & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(OFF_DIAG_SUB_MAT_LOC(1,xi)+ns,OFF_DIAG_SUB_MAT_LOC(2,xi)+ms)+ &
                      & DOT_PRODUCT(DPHIDX_COMP(OFF_DIAG_DEP_VAR(1,1,xi))%DPHIDX(ns,OFF_DIAG_DEP_VAR(1,:,xi))* &
                      &  DPHIDX_COMP(OFF_DIAG_DEP_VAR(1,2,xi))%DPHIDX(ms,OFF_DIAG_DEP_VAR(2,:,xi)),JRWG_OFF_DIAG_C(:,xi))
                  ENDDO !ms
                ENDDO !ns
              ENDDO !xi

             !Below is the full form of constructing the off-Diagonal terms. This will be documented in the linear elasticity equation set page on doxygen for clarity

             !Expanding the DOT_PRODUCT terms

             !OFF_DIAG_DEP_VAR(1,:) = (/1,1,2/)
             !OFF_DIAG_DEP_VAR(2,:) = (/2,3,3/)
             ! DO xi=1,OFF_DIAG_COMP(NUMBER_OF_XI)
             !   DO ns=1,DEPENDENT_BASES_EP(OFF_DIAG_DEP_VAR(1,1,xi))
             !     DO ms=1,DEPENDENT_BASES_EP(OFF_DIAG_DEP_VAR(1,2,xi))
             !       EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(OFF_DIAG_SUB_MAT_LOC(1,xi)+ns,OFF_DIAG_SUB_MAT_LOC(2,xi)+ms) = &
             !         & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(OFF_DIAG_SUB_MAT_LOC(1,xi)+ns,OFF_DIAG_SUB_MAT_LOC(2,xi)+ms)+ &
             !         & JRWG_OFF_DIAG_C(1,xi)*DPHIDX_COMP(OFF_DIAG_DEP_VAR(1,xi))%DPHIDX(ns,OFF_DIAG_DEP_VAR(1,xi))* &
             !         & DPHIDX_COMP(OFF_DIAG_DEP_VAR(2,xi))%DPHIDX(ms,OFF_DIAG_DEP_VAR(2,xi)) + &
             !         & JRWG_OFF_DIAG_C(2,xi)*DPHIDX_COMP(OFF_DIAG_DEP_VAR(1,xi))%DPHIDX(ns,OFF_DIAG_DEP_VAR(2,xi))* &
             !         & DPHIDX_COMP(OFF_DIAG_DEP_VAR(2,xi))%DPHIDX(ms,OFF_DIAG_DEP_VAR(1,xi))
             !     ENDDO !ms
             !   ENDDO !ns
             ! ENDDO !xi

             ! Expanding the xi loop above

!             DO ns=1,DEPENDENT_BASES_EP(1)
!               DO ms=1,DEPENDENT_BASES_EP(2)
!                 EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(ns,ms+DEPENDENT_BASES_EP(1)) = &
!                   & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(ns,ms+DEPENDENT_BASES_EP(1)) + &
!                   & JRWG*C(1,2)*DPHIDX_COMP(1)%DPHIDX(ns,1)*DPHIDX_COMP(2)%DPHIDX(ms,2) + &
!                   & JRWG*C(6,6)*DPHIDX_COMP(1)%DPHIDX(ns,2)*DPHIDX_COMP(2)%DPHIDX(ms,1)
!               ENDDO !ns
!             ENDDO !ms
!             DO ns=1,DEPENDENT_BASES_EP(1)
!               DO ms=1,DEPENDENT_BASES_EP(3)
!                 EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(ns,ms+DEPENDENT_BASES_EP(1)+DEPENDENT_BASES_EP(2)) = &
!                   & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(ns,ms+DEPENDENT_BASES_EP(1)+DEPENDENT_BASES_EP(2)) + &
!                   & JRWG*C(1,3)*DPHIDX_COMP(1)%DPHIDX(ns,1)*DPHIDX_COMP(3)%DPHIDX(ms,3) + &
!                   & JRWG*C(4,4)*DPHIDX_COMP(1)%DPHIDX(ns,3)*DPHIDX_COMP(3)%DPHIDX(ms,1)
!               ENDDO !ns
!             ENDDO !ms
!             DO ns=1,DEPENDENT_BASES_EP(2)
!               DO ms=1,DEPENDENT_BASES_EP(3)
!                 EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(ns+DEPENDENT_BASES_EP(1),ms+DEPENDENT_BASES_EP(1)+DEPENDENT_BASES_EP(2)) = &
!                   & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(ns+DEPENDENT_BASES_EP(1),ms+DEPENDENT_BASES_EP(1)+DEPENDENT_BASES_EP(2)) + &
!                   & JRWG*C(2,3)*DPHIDX_COMP(2)%DPHIDX(ns,2)*DPHIDX_COMP(3)%DPHIDX(ms,3) + &
!                   & JRWG*C(5,5)*DPHIDX_COMP(2)%DPHIDX(ns,3)*DPHIDX_COMP(3)%DPHIDX(ms,2)
!               ENDDO !ns
!             ENDDO !ms
            ENDDO !ng
            !If Plane Stress/Strain problem multiply equation matrix by thickness
            IF(EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE) THEN
              DO mhs=1,TOTAL_DEPENDENT_BASIS_EP
                DO nhs=mhs,TOTAL_DEPENDENT_BASIS_EP
                  !!TODO::Bring 2D plane stress/strain element thickness in through a field - element constant when it can be exported by field i/o. Currently brought in through material field (Temporary)
                  EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                    & MATERIALS_INTERP_POINT%values(3,1)
                ENDDO !nhs
              ENDDO !mhs
            ENDIF
          ENDIF

          !!TODO:: Is this RHS Vector update required? find out/check - RHS not used - BC are prescribed during assembling eg update RHS only when BC change - stiffness matrix should be the same
          IF(RHS_VECTOR%UPDATE_VECTOR) THEN
            RHS_VECTOR%ELEMENT_VECTOR%VECTOR=0.0_DP
          ENDIF

          !Scale factor adjustment, Application of Scale factors is symmetric
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,DEPENDENT_INTERPOLATION_PARAMETERS, &
              & ERR,ERROR,*999)
            DO xi=1,NUMBER_OF_XI
              SF(DIAG_SUB_MAT_LOC(xi)+1:SUM(DEPENDENT_BASES_EP(1:xi)))=DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(:,xi)
            ENDDO !xi
            DO mhs=1,TOTAL_DEPENDENT_BASIS_EP
              IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                DO nhs=mhs,TOTAL_DEPENDENT_BASIS_EP
                  EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)*SF(mhs)*SF(nhs)
                ENDDO !nhs
              ENDIF
              !!TODO:: Check if RHS update required for Linear Elasticity ie is the RHS the force terms but they are set during assembling and not here?
              IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)*SF(mhs)
            ENDDO !mhs
          ENDIF

          IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
            !Transpose upper triangular portion of Stiffness matrix to give lower triangular portion. Has to be done after scale factors are applied
            !!TODO:: Use symmetric linear equation solver or alternatively traspose to give full matrix when asemmbling or when creating solver matrices
            !!TODO:: Better to use SIZE(EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX,1) as apposed to TOTAL_DEPENDENT_BASIS_EP? Is the size re-calculated at end of every loop?
            DO mhs=2,TOTAL_DEPENDENT_BASIS_EP
              DO nhs=1,mhs-1
                EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) = EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(nhs,mhs)
              ENDDO !nhs
            ENDDO !mhs
          ENDIF

        CASE(EQUATIONS_SET_PLATE_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_SHELL_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Linear Elasticity equation type of a Elasticty equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the linear elasticity tensor
  SUBROUTINE LINEAR_ELASTICITY_TENSOR(EQUATIONS_SET_SUBTYPE,MATERIALS_INTERPOLATED_POINT,ELASTICITY_TENSOR,ERR,ERROR,*)

    !Argument variables    
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The subtype of the particular equation set being used
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: MATERIALS_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: ELASTICITY_TENSOR(:,:) !<The Linear Elasticity Tensor C
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    REAL(DP) :: E1,E2,E3,v13,v23,v12,v31,v32,v21,gama
    REAL(DP) :: C11,C22,C33,C12,C13,C23,C21,C31,C32,C44,C55,C66
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LINEAR_ELASTICITY_TENSOR",ERR,ERROR,*999)

    SELECT CASE(EQUATIONS_SET_SUBTYPE)
    !Note: Fortran uses column major format for arrays.
    CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
      !General Orthotropic 3D Linear Elasticity Tensor
      E1 = MATERIALS_INTERPOLATED_POINT%values(1,1)
      E2 = MATERIALS_INTERPOLATED_POINT%values(2,1)
      E3 = MATERIALS_INTERPOLATED_POINT%values(3,1)
      v13 = MATERIALS_INTERPOLATED_POINT%values(4,1)
      v23 = MATERIALS_INTERPOLATED_POINT%values(5,1)
      v12 = MATERIALS_INTERPOLATED_POINT%values(6,1)
      v31 = v13
      v32 = v23
      v21 = v12
      gama = 1.0_DP/(1.0_DP-v12*v21-v23*v32-v31*v13-2.0_DP*v21*v32*v13)
      C11 = E1*(1.0_DP-v23*v32)*gama
      C22 = E2*(1.0_DP-v13*v31)*gama
      C33 = E3*(1.0_DP-v12*v21)*gama
      C12 = E1*(v21+v31*v23)*gama ! = E2*(v12+v32*v13)*gama
      C13 = E1*(v31+v21*v32)*gama ! = E3*(v13+v12*v23)*gama
      C23 = E2*(v32+v12*v31)*gama ! = E3*(v23+v21*v13)*gama
      C21 = C12
      C31 = C13
      C32 = C23
      C44 = E2/(2.0_DP*(1.0_DP+v23)) != G23
      C55 = E1/(2.0_DP*(1.0_DP+v13)) != G13
      C66 = E3/(2.0_DP*(1.0_DP+v12)) != G12
      ELASTICITY_TENSOR(1:6,1)=(/C11,C21,C31,0.0_DP,0.0_DP,0.0_DP/)
      ELASTICITY_TENSOR(1:6,2)=(/C12,C22,C32,0.0_DP,0.0_DP,0.0_DP/)
      ELASTICITY_TENSOR(1:6,3)=(/C13,C23,C33,0.0_DP,0.0_DP,0.0_DP/)
      ELASTICITY_TENSOR(1:6,4)=(/0.0_DP,0.0_DP,0.0_DP,C44,0.0_DP,0.0_DP/)
      ELASTICITY_TENSOR(1:6,5)=(/0.0_DP,0.0_DP,0.0_DP,0.0_DP,C55,0.0_DP/)
      ELASTICITY_TENSOR(1:6,6)=(/0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,C66/)
    CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE)
      !Plane Stress Isotropic Elasticity Tensor
      E1 = MATERIALS_INTERPOLATED_POINT%values(1,1)
      v12 = MATERIALS_INTERPOLATED_POINT%values(2,1)
      v21 = v12
      gama = 1.0_DP/(1.0_DP-v12*v21)
      C11 = E1*gama
      C22 = C11
      C12 = C11*v21
      C21 = C12
      C66 = E1/(2.0_DP*(1.0_DP+v12)) != G12
      ELASTICITY_TENSOR=0.0_DP
      ELASTICITY_TENSOR(1,1)=C11
      ELASTICITY_TENSOR(1,2)=C21
      ELASTICITY_TENSOR(2,1)=C21
      ELASTICITY_TENSOR(2,2)=C22
      ELASTICITY_TENSOR(6,6)=C66
    CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
      !Plane Strain Isotropic Linear Elasticity Tensor
      E1 = MATERIALS_INTERPOLATED_POINT%values(1,1)
      E2 = E1
      v12 = MATERIALS_INTERPOLATED_POINT%values(2,1)
      v21 = v12
      gama = 1.0_DP/(1.0_DP-v12-v21)
      C11 = E1*gama*(1.0_DP-v12)/(1.0_DP+v12)
      C22 = E2*gama*(1.0_DP-v21)/(1.0_DP+v21)
      C12 = C22*v12
      C21 = C12
      C66 = E1/(2.0_DP*(1.0_DP+v12)) != G12
      ELASTICITY_TENSOR=0.0_DP
      ELASTICITY_TENSOR(1,1)=C11
      ELASTICITY_TENSOR(1,2)=C21
      ELASTICITY_TENSOR(2,1)=C21
      ELASTICITY_TENSOR(2,2)=C22
      ELASTICITY_TENSOR(6,6)=C66
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
      !!TODO:: Set this up correctly
      !Plane Strain Isotropic Linear Elasticity Tensor
      E1 = MATERIALS_INTERPOLATED_POINT%values(1,1)
      v12 = MATERIALS_INTERPOLATED_POINT%values(2,1)
      gama = 1.0_DP/(1.0_DP-v12-v21)
      C11 = E1*gama*(1.0_DP-v12)/(1.0_DP+v12)
      C66 = E1/(2.0_DP*(1.0_DP+v12)) != G12
      ELASTICITY_TENSOR=0.0_DP
      ELASTICITY_TENSOR(1,1)=C11
      ELASTICITY_TENSOR(6,6)=C66
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Linear Elasticity equation type of a Elasticty equations set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    CALL EXITS("LINEAR_ELASTICITY_TENSOR")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_TENSOR",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_TENSOR")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_TENSOR

  !
  !================================================================================================================================
  !

  !>Sets up the Linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LINEAR_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a linear elasticity equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE, &
      & NUMBER_OF_COMPONENTS,NUMBER_OF_DIMENSIONS
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LINEAR_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN

      !!TODO:: Update all these so there is a default setup that is valid for 1D,2D & 3D Linear Elasticity
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      !
      ! THREE DIMENSIONAL ELASTICITY
      !
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Default to FEM solution
            CALL LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticity equation."
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
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              DO component_idx=1,NUMBER_OF_DIMENSIONS
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              ENDDO !component_idx
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
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
              & " is invalid for a linear elasticity equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Default to the general 3D orthotropic material
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
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 6,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)                
                DO component_idx=1,6
                  !Default to to the first geometric component with constant interpolation
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,6,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                DO component_idx=1,3
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,30.0E6_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                DO component_idx=4,6
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,0.25_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticity equation."
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
              & " is invalid for a linear elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
!        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
!          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
!          CASE(EQUATIONS_SET_SETUP_START_ACTION)
!            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
!              !Do nothing
!            ELSE
!              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
!            ENDIF
!          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!            !Do nothing
!            !? Maybe set finished flag????
!          CASE DEFAULT
!            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
!              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
!              & " is invalid for a linear elasticity equation."
!            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
!          END SELECT

        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  !List 3 Dimensional Analytic function types currently implemented
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_THREE_DIM_1)
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
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_THREE_DIM_1
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a standard Linear Elasticity equation."
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
          CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
                  CALL LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set analtyic has not been finished.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set analtyic is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Linear Elasticity equation."
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
              CALL FLAG_ERROR("Equations set dependent field has not bee finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
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
              & " is invalid for a linear elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
            IF(ASSOCIATED(EQUATIONS)) THEN
              IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set equations has not been finished.",ERR,ERROR,*999)               
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticiy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear elasticity equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !
      ! TWO DIMENSIONAL PLANE STRESS & PLANE STRAIN ELASTICITY
      !
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE,EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Default to FEM solution
            CALL LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticity equation."
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
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              DO component_idx=1,NUMBER_OF_DIMENSIONS
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              ENDDO !component_idx
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
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
              & " is invalid for a linear elasticity equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Default to the general 3D orthotropic material
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
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 6,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)                
                DO component_idx=1,6
                  !Default to to the first geometric component with constant interpolation
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
                !2 Components for 2D Isotropic Linear Elasticity
                !TODO:: Temporarily set to 3 to allow thickness to passed in. Remove once a thickness, element constant field is defined and can be exported/viewed by cmgui
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,3,ERR,ERROR,*999) 
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,30.0E6_DP,ERR,ERROR,*999) !Young's Modulus
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,2,0.25_DP,ERR,ERROR,*999) !Poisson's Ratio
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticity equation."
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
              & " is invalid for a linear elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
!        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
!          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
!          CASE(EQUATIONS_SET_SETUP_START_ACTION)
!            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
!              !Do nothing
!            ELSE
!              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
!            ENDIF
!          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!            !Do nothing
!            !? Maybe set finished flag????
!          CASE DEFAULT
!            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
!              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
!              & " is invalid for a linear elasticity equation."
!            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
!          END SELECT

        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  !List 3 Dimensional Analytic function types currently implemented
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_TWO_DIM_PLANE_STRESS_1)
                    !Check that we are in 2D
                    !!TODO:: This check may have been done before
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LINEAR_ELASTICITY_EQUATION_TWO_DIM_PLANE_STRESS_1
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a standard Linear Elasticity equation."
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
          CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
                  CALL LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set analtyic has not been finished.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set analtyic is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Linear Elasticity equation."
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
              CALL FLAG_ERROR("Equations set dependent field has not bee finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
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
              & " is invalid for a linear elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
            IF(ASSOCIATED(EQUATIONS)) THEN
              IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set equations has not been finished.",ERR,ERROR,*999)               
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticiy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear elasticity equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !
      ! ONE DIMENSIONAL ELASTICITY
      !
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Default to FEM solution
            CALL LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticity equation."
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
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              DO component_idx=1,NUMBER_OF_DIMENSIONS
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              ENDDO !component_idx
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
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
              & " is invalid for a linear elasticity equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Default to the general 3D orthotropic material
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
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 6,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)                
                DO component_idx=1,6
                  !Default to to the first geometric component with constant interpolation
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
                !2 Components for 2D Isotropic Linear Elasticity
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2,ERR,ERROR,*999) 
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,30.0E6_DP,ERR,ERROR,*999) !Young's Modulus
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,2,0.25_DP,ERR,ERROR,*999) !Poisson's Ratio
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticity equation."
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
              & " is invalid for a linear elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Do nothing
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticity equation."
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
              CALL FLAG_ERROR("Equations set dependent field has not bee finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
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
              & " is invalid for a linear elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
            IF(ASSOCIATED(EQUATIONS)) THEN
              IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set equations has not been finished.",ERR,ERROR,*999)               
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear elasticiy equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear elasticity equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_PLATE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_SHELL_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a linear elasticity equation type of an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LINEAR_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
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
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE)
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
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
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
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
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
      CASE(EQUATIONS_SET_PLATE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_SHELL_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a linear elasticity equation type of an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    CALL EXITS("LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LINEAR_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_FIELD_USER_NUMBER
    TYPE(FIELD_TYPE), POINTER :: DUMMY_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LINEAR_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DUMMY_FIELD_USER_NUMBER=0
      NULLIFY(DUMMY_FIELD)
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)        
        EQUATIONS_SET%CLASS=EQUATIONS_SET_ELASTICITY_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_LINEAR_ELASTICITY_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE  
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_ELASTICITY_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_LINEAR_ELASTICITY_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_ELASTICITY_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_LINEAR_ELASTICITY_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_ELASTICITY_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_LINEAR_ELASTICITY_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE
      CASE(EQUATIONS_SET_PLATE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_SHELL_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a linear elasticity equation type of an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LINEAR_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the linear elasticity problem.
  SUBROUTINE LINEAR_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a linear elasticity equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LINEAR_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE)
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
              & " is invalid for a linear elasticity problem."
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
              & " is invalid for a linear elasticity problem."
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
              & " is invalid for a linear elasticity problem."
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
              & " is invalid for a linear elasticity problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear elasticity problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a linear elasticity type of an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LINEAR_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a linear elasticity type .
  SUBROUTINE LINEAR_ELASTICITY_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LINEAR_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_ELASTICITY_CLASS
        PROBLEM%TYPE=PROBLEM_LINEAR_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NO_SUBTYPE      
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a linear elasticity type of an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LINEAR_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("LINEAR_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("LINEAR_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE LINEAR_ELASTICITY_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

END MODULE LINEAR_ELASTICITY_ROUTINES
