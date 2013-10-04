!> \file
!> \author Chris Bradley
!> \brief This module handles all diffusion equation routines.
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

!>This module handles all advection-diffusion equation routines.
MODULE ADVECTION_DIFFUSION_EQUATION_ROUTINES

  USE ANALYTIC_ANALYSIS_ROUTINES
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
! temporary input for setting velocity field
  USE FLUID_MECHANICS_IO_ROUTINES


  IMPLICIT NONE

  PRIVATE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_SETUP,ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & ADVECTION_DIFFUSION_EQUATION_ANALYTIC_CALCULATE, &
    & AdvectionDiffusion_EquationEquationsSetSpecificationSet,ADVECTION_DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE, &
    & AdvectionDiffusionEquation_ProblemSpecificationSet,ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP, &
    & ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA, ADVECTION_DIFFUSION_PRE_SOLVE,ADVECTION_DIFFUSION_POST_SOLVE, &
    & ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE,ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN
  
CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !>For the advection-diffusion analytic example it is required that the advective velocity
  !>and the source field are set to a particular analytic value, which is performed within this subroutine.
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type
    REAL(DP) :: VALUE,X(3),VALUE_SOURCE,VALUE_INDEPENDENT,VALUE_MATERIAL
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD,SOURCE_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP) :: alpha, phi, Peclet,tanphi

    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR,*999)

    NULLIFY(GEOMETRIC_VARIABLE)
    NULLIFY(GEOMETRIC_PARAMETERS)

    alpha = 1.0_DP
    phi = 0.2_DP
    tanphi = TAN(phi)
    Peclet= 10.0_DP

    !>Set the analytic boundary conditions 
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
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
                                !Default to version 1 of each node derivative
                                local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                  & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                              ENDDO !dim_idx
                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                                  !This is a steady-state solution of advection-diffusion equation
                                  !Velocity field takes form v(x,y)=(sin(6y),cos(6x))
                                  !Solution is u(x,y)=tanh(1 - alpha.(x.tan(Phi) - y))
                                 SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=TANH(1.0-alpha*(X(1)*tanphi-X(2)))
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
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
                                IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
                                  IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                    !If we are a boundary node then set the analytic value on the boundary
                                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                      & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
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

              ENDDO !variable_idx
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
    
    !>Set the independent field (i.e. the advective velocity) to a specified analytical function 
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
        IF(ASSOCIATED(INDEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)
            DO variable_idx=1,INDEPENDENT_FIELD%NUMBER_OF_VARIABLES
              variable_type=INDEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>INDEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
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
                              local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                              X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                            ENDDO !dim_idx
                            !Loop over the derivatives
                            DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                              SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                              CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                                !Velocity field takes form v(x,y)=(sin(6y),cos(6x))
                                IF(component_idx==1) THEN
                                  VALUE_INDEPENDENT=SIN(6*X(1))
                                ELSE
                                  VALUE_INDEPENDENT=COS(6*X(2))           
                                ENDIF
                              CASE DEFAULT
                                LOCAL_ERROR="The analytic function type of "// &
                                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                  & " is invalid."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              END SELECT
                              !Default to version 1 of each node derivative
                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,local_ny,VALUE_INDEPENDENT,ERR,ERROR,*999)
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
                CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
              ENDIF

            ENDDO !variable_idx
            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
          ENDIF            
        ELSE
          CALL FLAG_ERROR("Equations set independent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

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
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)
            DO variable_idx=1,SOURCE_FIELD%NUMBER_OF_VARIABLES
              variable_type=SOURCE_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
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
                              local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                              X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                            ENDDO !dim_idx
                            !Loop over the derivatives
                            DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                              SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                              CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                               VALUE_SOURCE= (1.0/Peclet)*(2.0*TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))*(1.0-(TANH(-0.1E1+alpha*( &
                                  & tanphi*X(1)-X(2)))**2))*alpha*alpha*tanphi*tanphi+2.0*TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)) &
                                  & )*(1.0-(TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))**2))*alpha*alpha-Peclet*(-SIN(6.0*X(2) &
                                  & )*(1.0-(TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))**2))*alpha*tanphi+COS(6.0*X(1))*(1.0- &
                                  & (TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))**2))*alpha))
                                CASE DEFAULT
                                LOCAL_ERROR="The analytic function type of "// &
                                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
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
            ENDDO !variable_idx
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

    !>Set the material field to a specified analytical value
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)
            DO variable_idx=1,MATERIALS_FIELD%NUMBER_OF_VARIABLES
              variable_type=MATERIALS_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>MATERIALS_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                   SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                   CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                     VALUE_MATERIAL= (1.0/Peclet)
                   CASE DEFAULT
                     LOCAL_ERROR="The analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                   END SELECT
                      CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                       & FIELD_VALUES_SET_TYPE,component_idx,VALUE_MATERIAL,ERR,ERROR,*999)
                ENDDO !component_idx
                CALL FIELD_PARAMETER_SET_UPDATE_START(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !variable_idx
            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
          ENDIF            
        ELSE
          CALL FLAG_ERROR("Equations set material field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF


    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_ANALYTIC_CALCULATE


  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equation type of a classical field equations set class.
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a diffusion equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for an advection-diffusion equation type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a diffusion equation type of an classical field equations set class.
  SUBROUTINE ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
         & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
         & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
         & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, & 
         & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)        
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
          & " is not valid for an advection-diffusion equation type of an classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a diffusion equation type of a classical field equations set class.
  SUBROUTINE AdvectionDiffusion_EquationEquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    CALL Enters("AdvectionDiffusion_EquationEquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for an advection-diffusion type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for an advection-diffusion equation type of a classical field equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    CALL Exits("AdvectionDiffusion_EquationEquationsSetSpecificationSet")
    RETURN
999 CALL Errors("AdvectionDiffusion_EquationEquationsSetSpecificationSet",err,error)
    CALL Exits("AdvectionDiffusion_EquationEquationsSetSpecificationSet")
    RETURN 1
  END SUBROUTINE AdvectionDiffusion_EquationEquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the linear advection-diffusion equation.
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS, NUMBER_OF_SOURCE_COMPONENTS, NUMBER_OF_INDEPENDENT_COMPONENTS,imy_matrix,Ncompartments,&
      & GEOMETRIC_COMPONENT_NUMBER,NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS,NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SOURCE
    TYPE(EQUATIONS_SET_INDEPENDENT_TYPE), POINTER :: EQUATIONS_INDEPENDENT
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,EQUATIONS_SET_FIELD_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: num_var,num_var_count,NUMBER_OF_MATERIALS_COUPLING_COMPONENTS    
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS    
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:),VARIABLE_U_TYPES(:),COUPLING_MATRIX_STORAGE_TYPE(:), &
      & COUPLING_MATRIX_STRUCTURE_TYPE(:)
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE
    
    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_EQUATION_SET_LINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL ADVEC_DIFF_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET, &
              & EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
             IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES, &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              ENDIF
             ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
             IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1, 1_INTG, ERR, ERROR, *999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 2, 1_INTG, ERR, ERROR, *999)
              ENDIF
             ENDIF
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
            CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
              SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
              CASE(EQUATIONS_SET_SETUP_START_ACTION)
                EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
                EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
                IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                    & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,& 
                    & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)                
                  DO component_idx = 1, EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                ELSE
                  !Do nothing
                ENDIF
              CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                ! do nothing
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a linear advection-diffusion equation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
               !Do nothing???
            CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
               CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE, ERR, ERROR, *999)
               CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE, ERR, ERROR, *999)
            END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
               IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                 !Create the auto created dependent field
                 CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                   & DEPENDENT_FIELD,ERR,ERROR,*999)
                 CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
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
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                   & FIELD_DELUDELN_VARIABLE_TYPE,1,ERR,ERROR,*999)
                 !Default to the geometric interpolation setup
                 CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
                   & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
                !uses number of compartments to check that appropriate number and type of variables have been set on the dependent field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)                 
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2*Ncompartments,ERR,ERROR,*999)
                  !Create & populate array storing all of the relevant variable types against which to check the field variables
                  ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                  DO num_var=1,Ncompartments
                    VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                    VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  ENDDO
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,ERR,ERROR,*999)

                  DO num_var=1,2*Ncompartments
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var), & 
                       & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),1, &
                       & ERR,ERROR,*999)
                  ENDDO
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                      component_idx=1
                    DO num_var=1,2*Ncompartments
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),component_idx, &
                        & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    ENDDO
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
                 !Check the user specified field
                 CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
                 CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                   & ERR,ERROR,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                   & ERR,ERROR,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                   & ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                   & ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                   & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
!                   DO component_idx=1,NUMBER_OF_DIMENSIONS
                     component_idx=1
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                       & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                       & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
!                   ENDDO !component_idx
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
                END SELECT 
               ENDIF
            CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)

               IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                 !Create the auto created dependent field
                 CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                   & DEPENDENT_FIELD,ERR,ERROR,*999)
                 CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
                 CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                 CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                 CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                   & ERR,ERROR,*999)
                 CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                   & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & 1,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & 1,ERR,ERROR,*999)
                 !Default to the geometric interpolation setup
                 CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_DELVDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
                 CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE,&
                  & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE/),ERR,ERROR,*999)

                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                   & ERR,ERROR,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                   & ERR,ERROR,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                   & ERR,ERROR,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                   & ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                   & 1,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                   & ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, & 
                   & 1,ERR,ERROR,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
 !                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                     component_idx=1
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                       & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                       & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                       & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, & 
                       & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
!                   ENDDO !component_idx
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
              END SELECT
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear advection-diffusion equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                    & MATERIALS_FIELD,ERR,ERROR,*999)
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",ERR,ERROR,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & ERR,ERROR,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,2,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE/), &
                    & ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                   !Set the number of materials components
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                    Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                    NUMBER_OF_MATERIALS_COUPLING_COMPONENTS=Ncompartments
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COUPLING_COMPONENTS,ERR,ERROR,*999)
                  !Default the k materials components to the geometric interpolation setup with constant interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  ENDDO !component_idx
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  DO component_idx=1,NUMBER_OF_MATERIALS_COUPLING_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  ENDDO 
                    !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
               ELSE !standard materials field
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",ERR,ERROR,*999)
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
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                  !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                ENDIF
                 !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)   
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDIF
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
               ENDIF
              ELSE
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE, &
                     & FIELD_V_VARIABLE_TYPE/),ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                    & ERR,ERROR,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                    Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                    NUMBER_OF_MATERIALS_COUPLING_COMPONENTS=Ncompartments
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COUPLING_COMPONENTS,ERR,ERROR,*999)
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
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                    & ERR,ERROR,*999)
                ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                    & ERR,ERROR,*999)
                ENDIF
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
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS             
                ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                  !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  !WRITE(*,'("Setting materials components values :")')
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                    Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  DO component_idx=1,Ncompartments
                   CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                     & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,ERR,ERROR,*999)
                   ENDDO !component_idx
                ENDIF
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Now set the linear source values to 1.0
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Create the auto created source field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SOURCE% &
                  & SOURCE_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,"Source Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                NUMBER_OF_SOURCE_COMPONENTS=1
                !Set the number of source components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_SOURCE_COMPONENTS,ERR,ERROR,*999)
                !Default the source components to the geometric interpolation setup with constant interpolation
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &  
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDIF
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Finish creating the source field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                !Set the default values for the source field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  NUMBER_OF_SOURCE_COMPONENTS=1
                ELSE
                  NUMBER_OF_SOURCE_COMPONENTS=0
                ENDIF
                !Now set the source values to 1.0
                DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          !Setup the equations set for the advective velocity field
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_INDEPENDENT=>EQUATIONS_SET%INDEPENDENT
            IF(ASSOCIATED(EQUATIONS_INDEPENDENT)) THEN
             IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              IF(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_INDEPENDENT% &
                  & INDEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,2,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS,ERR,ERROR,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS=Ncompartments-1
                 !Set the number of independent components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS,ERR,ERROR,*999)
                !Default the k independent components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
                DO component_idx=1,NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE/),&
                   & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                   & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                   & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                  & ERR,ERROR,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS=Ncompartments-1
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS,ERR,ERROR,*999)
                DO component_idx=1,NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO
                DO component_idx=1,NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO
              ENDIF
             ELSE
              IF(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_INDEPENDENT% &
                  & INDEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  NUMBER_OF_INDEPENDENT_COMPONENTS=NUMBER_OF_DIMENSIONS
                 !Set the number of independent components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_COMPONENTS,ERR,ERROR,*999)
                !Default the k independent components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
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
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                  & ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
              ENDIF
             ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set independent is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_INDEPENDENT=>EQUATIONS_SET%INDEPENDENT
            IF(ASSOCIATED(EQUATIONS_INDEPENDENT)) THEN
              IF(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Finish creating the independent field
                CALL FIELD_CREATE_FINISH(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                !Set the default values for the independent field
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)

                 IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                  & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
!                   CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
!                      & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)
                 ENDIF
!                 CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
!                   NUMBER_OF_INDEPENDENT_COMPONENTS=NUMBER_OF_DIMENSIONS             
!                 !First set the k values to 1.0
!                 DO component_idx=1,NUMBER_OF_INDEPENDENT_COMPONENTS
!                   CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
!                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set independent is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c   f i e l d
        !-----------------------------------------------------------------
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
                  CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                      IF(NUMBER_OF_DIMENSIONS/=2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                       & EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a linear advection-diffusion equation."
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
              & " is invalid for a linear advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
              ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set subtype not valid.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
              !Finish the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              SELECT CASE(EQUATIONS_SET_SUBTYPE)
              CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN             
                CALL EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              ENDIF
              CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
                 & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)              
                 EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                 CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                 imy_matrix = EQUATIONS_SET_FIELD_DATA(1)
                 Ncompartments = EQUATIONS_SET_FIELD_DATA(2)    
                 CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                 CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
                 CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,Ncompartments-1,ERR,ERROR,*999)

                 ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                 ALLOCATE(VARIABLE_U_TYPES(Ncompartments-1))
                 DO num_var=1,Ncompartments
                   VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                   VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                 ENDDO
                 num_var_count=0
                 DO num_var=1,Ncompartments
                   IF(num_var/=imy_matrix)THEN
                     num_var_count=num_var_count+1
                     VARIABLE_U_TYPES(num_var_count)=VARIABLE_TYPES(2*num_var-1)
                   ENDIF
                 ENDDO
                 CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,VARIABLE_TYPES(2*imy_matrix-1),ERR,ERROR,*999)
                 CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,VARIABLE_U_TYPES,ERR,ERROR,*999)
                 CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,VARIABLE_TYPES(2*imy_matrix),ERR,ERROR,*999)
                 CALL EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              END SELECT
             CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              !Set up matrix storage and structure
              IF(EQUATIONS%LUMPING_TYPE==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET(EQUATIONS_MATRICES, &
                  & (/EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED/),ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                  & (/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE/),ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                  (/EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE/),ERR,ERROR,*999)
              ELSE
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & (/DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE/),ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & (/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                    & ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                    (/EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE/),ERR,ERROR,*999)    
                  IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                    ALLOCATE(COUPLING_MATRIX_STORAGE_TYPE(Ncompartments-1))
                    ALLOCATE(COUPLING_MATRIX_STRUCTURE_TYPE(Ncompartments-1))
                    DO num_var=1,Ncompartments-1
                     COUPLING_MATRIX_STORAGE_TYPE(num_var)=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                     COUPLING_MATRIX_STRUCTURE_TYPE(num_var)=EQUATIONS_MATRIX_FEM_STRUCTURE
                    ENDDO                    
                    CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                      & COUPLING_MATRIX_STORAGE_TYPE, &
                      & ERR,ERROR,*999)      
                    CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                      COUPLING_MATRIX_STRUCTURE_TYPE,ERR,ERROR,*999)
                  ENDIF                          
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
              CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              !Finish the creation of the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN              
                CALL EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              ENDIF
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
            ENDIF
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
              & " is invalid for a linear advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear advection-diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not a linear advection-diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP

  !
  !===============================================================================================================================
  !
!   SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)
!     !Argument variables
!     TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
!     TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
!     INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
! 
!     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!     CALL EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
!     RETURN
! 999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP",ERR,ERROR)
!     CALL EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
!     RETURN 1
! 
!   END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion problem.
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a diffusion equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for an advection-diffusion equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) mh,mhs,ms,ng,nh,nhs,ni,nj,ns,FIELD_VAR_TYPE,my_compartment,Ncompartments,imatrix,num_var_count
    INTEGER(INTG) :: MESH_COMPONENT_1, MESH_COMPONENT_2
    REAL(DP) :: C_PARAM,K_PARAM,RWG,SUM,PGMJ(3),PGNJ(3),ADVEC_VEL,A_PARAM,COUPLING_PARAM,PGM,PGN
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1, DEPENDENT_BASIS_2
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: DAMPING_MATRIX,STIFFNESS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD,INDEPENDENT_FIELD,EQUATIONS_SET_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(FIELD_VARIABLE_PTR_TYPE) :: FIELD_VARIABLES(99)
    TYPE(EQUATIONS_MATRIX_PTR_TYPE) :: COUPLING_MATRICES(99) 
    INTEGER(INTG) :: FIELD_VAR_TYPES(99)
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME_1, QUADRATURE_SCHEME_2
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    LOGICAL :: UPDATE_DAMPING_MATRIX,UPDATE_STIFFNESS_MATRIX,UPDATE_RHS_VECTOR,UPDATE_SOURCE_VECTOR
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE

    UPDATE_DAMPING_MATRIX = .FALSE.
    UPDATE_STIFFNESS_MATRIX = .FALSE.
    UPDATE_RHS_VECTOR = .FALSE.
    UPDATE_SOURCE_VECTOR = .FALSE.

    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
          CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
        END IF
        EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
        SELECT CASE(EQUATIONS_SET_SUBTYPE)
        CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD!Stores the advective velocity field
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD
          ENDIF
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
             STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
             DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
             STIFFNESS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          ELSE

          ENDIF
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN            
             SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          ENDIF
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             IF(ASSOCIATED(DAMPING_MATRIX)) UPDATE_DAMPING_MATRIX=DAMPING_MATRIX%UPDATE_MATRIX
          ENDIF
          IF(ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
          IF(ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN             
             IF(ASSOCIATED(SOURCE_VECTOR)) UPDATE_SOURCE_VECTOR=SOURCE_VECTOR%UPDATE_VECTOR
          ENDIF
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
           & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
           EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
           CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
             & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
            my_compartment = EQUATIONS_SET_FIELD_DATA(1)
            Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            num_var_count=0
            DO imatrix = 1,Ncompartments
             IF(imatrix/=my_compartment)THEN
              num_var_count=num_var_count+1
              COUPLING_MATRICES(num_var_count)%PTR=>LINEAR_MATRICES%MATRICES(num_var_count)%PTR
              FIELD_VARIABLES(num_var_count)%PTR=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(num_var_count)%VARIABLE
              FIELD_VAR_TYPES(num_var_count)=FIELD_VARIABLES(num_var_count)%PTR%VARIABLE_TYPE
              COUPLING_MATRICES(num_var_count)%PTR%ELEMENT_MATRIX%MATRIX=0.0_DP
             ENDIF
            END DO
          ENDIF
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &  
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
             FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
             FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          ELSE

          ENDIF
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
               & MATERIALS_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          ENDIF

          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & SOURCE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          ENDIF
          !the following line has been changed to use fieldinputdata1settype
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & INDEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN  
           DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS=> &
            & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR
           CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
           DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT=> &
            & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR
          ENDIF


          !Select whether using standard Galerkin scheme, or the stabilised streamwise-upwinding Petrov-Galerkin scheme
          SELECT CASE(EQUATIONS_SET_SUBTYPE)
          CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,& 
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE , & 
               & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Interpolate to get the advective velocity
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)            
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)   
            ENDIF
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                 & MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            ENDIF
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components

            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(UPDATE_STIFFNESS_MATRIX .OR. UPDATE_DAMPING_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(UPDATE_STIFFNESS_MATRIX) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          !Advection term is constructed here and then added to SUM for updating the stiffness matrix outside of this loop
                          ADVEC_VEL=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+ADVEC_VEL*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*PGNJ(nj)   
                        ENDDO !nj
                        IF (EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                        ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. & 
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                          A_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS,NO_PART_DERIV)
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                            ! A_PARAM is the material parameter that multiplies the linear source u
                        ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                        !for multi-compartment model must include additional terms into the
                          COUPLING_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR% &
                            & VALUES(my_compartment,NO_PART_DERIV)
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ & 
                            & SUM*RWG + QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG*COUPLING_PARAM
                        ENDIF
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                      IF(UPDATE_DAMPING_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
               ENDDO !ms
              ENDDO !mh
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN 
                IF(UPDATE_SOURCE_VECTOR) THEN
                    C_PARAM=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1, NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                       mhs=mhs+1
                       SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                IF(UPDATE_SOURCE_VECTOR) THEN
                    CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                          & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,ERR,ERROR,*999)
                    C_PARAM=DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                       mhs=mhs+1
                       SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ENDIF
            IF(UPDATE_RHS_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
  
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
            !Calculate the coupling matrices

              !Loop over element rows
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS !field_variable is the variable associated with the equations set under consideration

                MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS_1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%PTR% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME_1 => DEPENDENT_BASIS_1%QUADRATURE% &
                  & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                RWG = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN * &
                  & QUADRATURE_SCHEME_1%GAUSS_WEIGHTS(ng)

                DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1

                  num_var_count=0
                  DO imatrix = 1,Ncompartments
                  IF(imatrix/=my_compartment)THEN
                    num_var_count=num_var_count+1

!need to test for the case where imatrix==mycompartment
!the coupling terms then needs to be added into the stiffness matrix
                    IF(COUPLING_MATRICES(num_var_count)%PTR%UPDATE_MATRIX) THEN

!                       !Loop over element columns
                      nhs=0
! !                       DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      DO nh=1,FIELD_VARIABLES(num_var_count)%PTR%NUMBER_OF_COMPONENTS

                        MESH_COMPONENT_2 = FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS_2 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%PTR% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        !--- We cannot use two different quadrature schemes here !!!
                        QUADRATURE_SCHEME_2 => DEPENDENT_BASIS_2%QUADRATURE% &
                         & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                        !RWG = EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN * &
                        !  & QUADRATURE_SCHEME_2%GAUSS_WEIGHTS(ng)

                        DO ns=1,DEPENDENT_BASIS_2%NUMBER_OF_ELEMENT_PARAMETERS
                          nhs=nhs+1

!                           !-------------------------------------------------------------------------------------------------------------
!                           !concentration test function, concentration trial function
!                           !For now, this is only a dummy implementation - this still has to be properly set up.
!                           IF(mh==nh.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN ! don't need this for diffusion equation

!                             SUM = 0.0_DP

                            PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                            !Get the coupling coefficients 
                              COUPLING_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR% &
                                & VALUES(imatrix,NO_PART_DERIV)

!                              SUM = SUM + COUPLING_PARAM * PGM * PGN
 
                             COUPLING_MATRICES(num_var_count)%PTR%ELEMENT_MATRIX%MATRIX(mhs,nhs) = &
                               & COUPLING_MATRICES(num_var_count)%PTR%ELEMENT_MATRIX%MATRIX(mhs,nhs) + & 
                               & COUPLING_PARAM * PGM * PGN * RWG
!                           ENDIF
 
                        ENDDO !ns
                      ENDDO !nh
                    ENDIF
                   ENDIF
                  ENDDO !imatrix
                ENDDO !ms
              ENDDO !mh

            ENDIF

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
                IF(UPDATE_STIFFNESS_MATRIX .OR. UPDATE_DAMPING_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(UPDATE_STIFFNESS_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                      IF(UPDATE_DAMPING_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(UPDATE_RHS_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                IF(UPDATE_SOURCE_VECTOR) SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
              ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
          CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999) 
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Interpolate to get the advective velocity
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)            
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(UPDATE_STIFFNESS_MATRIX .OR. UPDATE_DAMPING_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(UPDATE_STIFFNESS_MATRIX) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          !Advection term is constructed here and then added to SUM for updating the stiffness matrix outside of this loop
                          ADVEC_VEL=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+ADVEC_VEL*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*PGNJ(nj)   
                        ENDDO !nj
                        IF (EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                        ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                          A_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS,NO_PART_DERIV)
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                            ! A_PARAM is the material parameter that multiplies the linear source u
                        ENDIF
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE) THEN
                      IF(UPDATE_DAMPING_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
               ENDDO !ms
              ENDDO !mh
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
                IF(UPDATE_SOURCE_VECTOR) THEN
                    C_PARAM=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1, NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                       mhs=mhs+1
                       SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ENDIF
            IF(UPDATE_RHS_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
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
                IF(UPDATE_STIFFNESS_MATRIX .OR. UPDATE_DAMPING_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(UPDATE_STIFFNESS_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE) THEN
                      IF(UPDATE_DAMPING_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(UPDATE_RHS_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                IF(UPDATE_SOURCE_VECTOR) SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
              ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
          CASE DEFAULT
            LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for an advection-diffusion equation type of a classical field equations set class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
          CALL FLAG_ERROR("Can not calculate finite element stiffness matrices for a nonlinear source.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for an advection-diffusion equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !INSERT CODE HERE TO DEAL WITH THE NON-LINEAR SOURCE TERMS - DECIDE HOW TO SOLVE THEM, USE JACOBIAN AS FOR NON-LINEAR POISSON EXAMPLE?

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a diffusion equation type.
  SUBROUTINE AdvectionDiffusionEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    CALL Enters("AdvectionDiffusionEquation_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
             & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
             & PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for an advection-diffusion equation type of a classical field problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Advection-diffusion problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    CALL Exits("AdvectionDiffusionEquation_ProblemSpecificationSet")
    RETURN
999 CALL Errors("AdvectionDiffusionEquation_ProblemSpecificationSet",err,error)
    CALL Exits("AdvectionDiffusionEquation_ProblemSpecificationSet")
    RETURN 1
  END SUBROUTINE AdvectionDiffusionEquation_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equations.
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    INTEGER(INTG) :: PROBLEM_SUBTYPE
    
    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
      END IF
      PROBLEM_SUBTYPE=PROBLEM%SPECIFICATION(3)
      IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
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
              & " is invalid for a linear advection-diffusion equation."
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
              & " is invalid for a linear advection-diffusion equation."
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
            !Set the solver to be a first order dynamic solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
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
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear advection-diffusion equation."
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
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
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
              & " is invalid for a linear advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear advection-diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSEIF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
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
              & " is invalid for a linear static advection-diffusion equation."
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
              & " is invalid for a linear static advection-diffusion equation."
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
              & " is invalid for a linear static advection-diffusion equation."
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
              & " is invalid for a linear static advection-diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear static advection-diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a linear advection-diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP
  
  !
  !================================================================================================================================
  !
    !>Sets up the diffusion equations.
!   SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)
! 
!     !Argument variables
!     TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
!     TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
!     INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
! 
!     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999) 
! 
!     CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
!     RETURN
! 999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",ERR,ERROR)
!     CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
!     RETURN 1
! 
!   END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Poisson problem pre solve.
  SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solvers
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: PROBLEM_SUBTYPE

    CALL ENTERS("ADVECTION_DIFFUSION_PRE_SOLVE",ERR,ERROR,*999)
    NULLIFY(SOLVER2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          PROBLEM_SUBTYPE=CONTROL_LOOP%PROBLEM%SPECIFICATION(3)
          IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
            !Update independent data fields
            CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            !CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
          ELSE IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
            !Update independent data fields
            CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            !CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CALL ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
          ELSE IF(PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN         
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          ELSE IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
            !do nothing
          ELSE
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a advection-diffusion type of a classical field problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END IF
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    END IF

    CALL EXITS("ADVECTION_DIFFUSION_PRE_SOLVE")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE
  !   
  !================================================================================================================================
  !
  !>Update mesh position and velocity for ALE advection-diffusion problem
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_DIFFUSION !<A pointer to the solvers
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)

    INTEGER(INTG) :: dof_number,TOTAL_NUMBER_OF_DOFS,NDOFS_TO_PRINT

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_DATA1(:)

    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH",ERR,ERROR,*999)

    NULLIFY(SOLVER_ALE_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
      ELSE IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL>1) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP%PARENT_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
      ENDIF
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
                        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
                          ! do nothing ???
                        CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Advection-diffusion update mesh ... ",ERR,ERROR,*999)
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                            !--- First, read mesh displacement values from file

                           CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)

                           INPUT_TYPE=42
                           INPUT_OPTION=2
                           NULLIFY(INPUT_DATA1)
                           !CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            !& FIELD_VALUES_SET_TYPE,INPUT_DATA1,ERR,ERROR,*999)
                           CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                            & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)

                            NULLIFY(MESH_DISPLACEMENT_VALUES)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                            IF(DIAGNOSTICS1) THEN
                              NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                                & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",3(X,E13.6))','3(3(X,E13.6))', &
                                & ERR,ERROR,*999)
                            ENDIF

                           CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                            & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)

                            TOTAL_NUMBER_OF_DOFS = GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                              & TOTAL_NUMBER_OF_DOFS

                            !--- Second, update geometric field
                            DO dof_number=1,TOTAL_NUMBER_OF_DOFS
                              CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(GEOMETRIC_FIELD, & 
                                & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, & 
                                & MESH_DISPLACEMENT_VALUES(dof_number), &
                                & ERR,ERROR,*999)
                            END DO
                            CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

                            !--- Third, use displacement values to calculate velocity values
                            ALPHA=1.0_DP/TIME_INCREMENT
                            CALL FIELD_PARAMETER_SETS_COPY(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,ERR,ERROR,*999)
                            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Geometric field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="Equations set subtype " &
                            & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
                            & " is not valid for an advection-diffusion equation type of a classical field problem class."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
                ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
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

    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH",ERR,ERROR)
    CALL EXITS("ADVECTION_ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH
  !   
  !================================================================================================================================
  !


  SUBROUTINE ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ADVECTION_DIFFUSION !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_ADVECTION_DIFFUSION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_ADVECTION_DIFFUSION !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_ADVECTION_DIFFUSION !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_ADVECTION_DIFFUSION !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION
    INTEGER(INTG) :: I

    CALL ENTERS("ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_ADVECTION_DIFFUSION)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the advection-diffusion equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Store value of advection-diffusion & 
                  & (dependent field - U variable_type) at time, t ... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
                SOLVER_EQUATIONS_ADVECTION_DIFFUSION=>SOLVER_ADVECTION_DIFFUSION%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)) THEN
                  SOLVER_MAPPING_ADVECTION_DIFFUSION=>SOLVER_EQUATIONS_ADVECTION_DIFFUSION%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_ADVECTION_DIFFUSION)) THEN
                    EQUATIONS_SET_ADVECTION_DIFFUSION=>SOLVER_MAPPING_ADVECTION_DIFFUSION%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_ADVECTION_DIFFUSION)) THEN
                      DEPENDENT_FIELD_ADVECTION_DIFFUSION=>EQUATIONS_SET_ADVECTION_DIFFUSION%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_ADVECTION_DIFFUSION)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_ADVECTION_DIFFUSION, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_ADVECTION_DIFFUSIONE is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Advection-diffusion equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Advection-diffusion solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Advection-diffusion solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,I,ERR,ERROR,*999)
                  END DO

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',ERR,ERROR,*999)
!                   CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
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



    CALL EXITS("ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN")
    RETURN
999 CALL ERRORS("ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN",ERR,ERROR)
    CALL EXITS("ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN")
    RETURN 1
  END SUBROUTINE ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN    
  !
  !================================================================================================================================
  !
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ADVECTION_DIFFUSION, SOLVER_DIFFUSION  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION, SOURCE_FIELD_ADVECTION_DIFFUSION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_ADVECTION_DIFFUSION, SOLVER_EQUATIONS_DIFFUSION  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_ADVECTION_DIFFUSION, SOLVER_MAPPING_DIFFUSION !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_ADVECTION_DIFFUSION, EQUATIONS_SET_DIFFUSION !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION
    INTEGER(INTG) :: I


    CALL ENTERS("ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_ADVECTION_DIFFUSION)
      NULLIFY(SOLVER_DIFFUSION)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Update advection-diffusion source field ... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DIFFUSION,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION=>SOLVER_DIFFUSION%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION)) THEN
                  SOLVER_MAPPING_DIFFUSION=>SOLVER_EQUATIONS_DIFFUSION%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION)) THEN
                    EQUATIONS_SET_DIFFUSION=>SOLVER_MAPPING_DIFFUSION%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION)) THEN
                      DEPENDENT_FIELD_DIFFUSION=>EQUATIONS_SET_DIFFUSION%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DIFFUSION, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_DIFFUSION is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Diffusion equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Diffusion solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Diffusion solver equations are not associated.",ERR,ERROR,*999)
                END IF


                !--- Get the source field for the advection-diffusion equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
                SOLVER_EQUATIONS_ADVECTION_DIFFUSION=>SOLVER_ADVECTION_DIFFUSION%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)) THEN
                  SOLVER_MAPPING_ADVECTION_DIFFUSION=>SOLVER_EQUATIONS_ADVECTION_DIFFUSION%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_ADVECTION_DIFFUSION)) THEN
                    EQUATIONS_SET_ADVECTION_DIFFUSION=>SOLVER_MAPPING_ADVECTION_DIFFUSION%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_ADVECTION_DIFFUSION)) THEN
                      SOURCE_FIELD_ADVECTION_DIFFUSION=>EQUATIONS_SET_ADVECTION_DIFFUSION%SOURCE%SOURCE_FIELD
                      IF(ASSOCIATED(SOURCE_FIELD_ADVECTION_DIFFUSION)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(SOURCE_FIELD_ADVECTION_DIFFUSION, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("SOURCE_FIELD_ADVECTION_DIFFUSION is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Advection-diffusion equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Advection-diffusion solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Advection-diffusion solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Copy the result from diffusion's dependent field to advection-diffusion's source field
                IF(NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION==NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION) THEN
                  DO I=1,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,SOURCE_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,ERR,ERROR,*999)
                  END DO
                ELSE
                  LOCAL_ERROR="Number of components of diffusion dependent field "// &
                    & "is not consistent with advection-diffusion equation source field."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END IF

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',ERR,ERROR,*999)
!                   CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
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

    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE
  !   
  !================================================================================================================================
  !
  !>Update independent field (velocity) for advection-diffusion pre solve
  SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

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

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_DATA1(:)
    INTEGER(INTG) :: PROBLEM_SUBTYPE

    CALL ENTERS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          PROBLEM_SUBTYPE=CONTROL_LOOP%PROBLEM%SPECIFICATION(3)
          IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & PROBLEM_SUBTYPE==PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read input data... ",ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  INPUT_TYPE=1
                  INPUT_OPTION=1
                  NULLIFY(INPUT_DATA1)
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VALUES_SET_TYPE,INPUT_DATA1,ERR,ERROR,*999)
                  CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                    & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
              END IF
           ELSE
             CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
           END IF
          ELSEIF(PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        END IF

        CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
           & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
        CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
           & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    END IF
    CALL EXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA

  !
  !================================================================================================================================
  !
  !Update the boundary conditions
  SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(CONTROL_LOOP,SOLVER,ERR,ERROR,*)
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
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT

!\todo: Reduce number of variable used
    INTEGER(INTG) :: BOUNDARY_CONDITION_CHECK_VARIABLE,node_idx
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: local_ny,global_ny
    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    INTEGER(INTG), POINTER :: BOUNDARY_NODES(:)

    CALL ENTERS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
      WRITE (*,*) CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE (PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
            & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
            & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                        & ERR,ERROR,*999)
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                         & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                        NULLIFY(BOUNDARY_VALUES)
                        NULLIFY(BOUNDARY_NODES)
                        CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS_ITERATION(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                          & BOUNDARY_NODES,NUMBER_OF_COMPONENTS,BOUNDARY_CONDITION_FIXED,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                          & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER)
                        WRITE(*,*) SIZE(BOUNDARY_VALUES)
                        DO node_idx=1,SIZE(BOUNDARY_VALUES)
                          !Default to version 1 of each node derivative
                          CALL FIELD_COMPONENT_DOF_GET_USER_NODE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & 1,NO_GLOBAL_DERIV,BOUNDARY_NODES(node_idx), &
                            & 1,local_ny,global_ny,ERR,ERROR,*999)
                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                            & CONDITION_TYPES(local_ny)
                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                              & BOUNDARY_VALUES(node_idx),ERR,ERROR,*999)
                          END IF
                        ENDDO
                      ELSE
                        CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Equations set dependent variable is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Boundary conditions are not associated.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
              END IF
            ELSE
              CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
            END IF
          ELSE
            CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
          END IF
          CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for an advection-diffusion equation of a classical field problem class."
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

    CALL EXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC")
    RETURN 1

  END SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC
  !
  !================================================================================================================================
  !
  SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solver
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ADVECTION_DIFFUSION_POST_SOLVE",ERR,ERROR,*999)
    NULLIFY(SOLVER2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
            CALL ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
          CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE,PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
            ! do nothing ???
            !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for an advection-diffusion type of a classical field problem class."
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

    CALL EXITS("ADVECTION_DIFFUSION_POST_SOLVE")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_POST_SOLVE",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_POST_SOLVE")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE
  !
  !================================================================================================================================
  !
  SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    CALL ENTERS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR

                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
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
!                        FILE="TRANSIENT_OUTPUT"
!                         METHOD="FORTRAN"
!                         EXPORT_FIELD=.TRUE.
!                         IF(EXPORT_FIELD) THEN          
!                          IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",ERR,ERROR,*999)
                          CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                              & ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
!                           ENDIF
!                         ENDIF 
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
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

    CALL EXITS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 CALL ERRORS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    CALL EXITS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA")
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA
  !
  !================================================================================================================================
  !
END MODULE ADVECTION_DIFFUSION_EQUATION_ROUTINES

