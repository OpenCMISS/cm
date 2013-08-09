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

!>This module handles all diffusion equation routines.
MODULE DIFFUSION_EQUATION_ROUTINES

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

  USE FLUID_MECHANICS_IO_ROUTINES

  IMPLICIT NONE

  PRIVATE
  
  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE,DIFFUSION_EQUATION_ANALYTIC_CALCULATE

  PUBLIC DIFFUSION_EQUATION_EQUATIONS_SET_SETUP

  PUBLIC DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET
  
  PUBLIC DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET

  PUBLIC DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE

  PUBLIC DIFFUSION_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE
  
  PUBLIC DIFFUSION_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE
  
  PUBLIC DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET

  PUBLIC DIFFUSION_EQUATION_PROBLEM_SETUP

  PUBLIC DIFFUSION_EQUATION_PRE_SOLVE,DIFFUSION_EQUATION_POST_SOLVE

!!TODO: should the following two routines really be public???

  PUBLIC DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE
  
  PUBLIC DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION
 
CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !Calculates a two-dimensional unsteady heat equation solution (diffusion coefficient is 1)
  SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type,version_idx
    REAL(DP) :: VALUE,X(3),INITIAL_VALUE
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:),GEOMETRIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    !TYPE(VARYING_STRING) :: LOCAL_ERROR    
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,imy_matrix
    !THESE ARE TEMPORARY VARIABLES - they need to be replace by constant field values and the current simulation time
    REAL(DP) :: TIME,NORMAL(3),TANGENTS(3,3)
    !CURRENT_TIME = 1.2_DP

    CALL ENTERS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
            ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
            ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            NULLIFY(GEOMETRIC_PARAMETERS)
            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)
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
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)THEN
                !If a multi-comp model, we will use the equations set field information to assign only the appropriate field
                !variable boundary conditions
                TIME=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)
                !Use predetermined mapping from equations set field compartment number to field variable type
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                imy_matrix = EQUATIONS_SET_FIELD_DATA(1)
                DO variable_idx=0,1
                  variable_type=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(imy_matrix-1))+variable_idx
                  !variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
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
                                  GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX
                                  CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET%SUBTYPE,ANALYTIC_FUNCTION_TYPE,&
                                    & X,TANGENTS,NORMAL,TIME,variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                    & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*999)
                                  !Default to version 1 of each node derivative
                                  local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                    & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                  IF(MOD(variable_type,FIELD_NUMBER_OF_VARIABLE_SUBTYPES)==FIELD_U_VARIABLE_TYPE) THEN
                                    IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                      !If we are a boundary node then set the analytic value on the boundary
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                        & FIELD_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
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
                    CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                      & ERR,ERROR,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                      & ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !variable_idx
                CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)            
              ELSE
                !for single physics diffusion problems use standard analytic calculate           
                ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                TIME=EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME
                DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                  variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                  FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                    CALL Field_ParameterSetEnsureCreated(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                      & ERR,ERROR,*999)
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
                                  GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX
                                  CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET%SUBTYPE,ANALYTIC_FUNCTION_TYPE,&
                                    & X,TANGENTS,NORMAL,0.0_DP,variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                    & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,INITIAL_VALUE,ERR,ERROR,*999)
                                  CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET%SUBTYPE,ANALYTIC_FUNCTION_TYPE,&
                                    & X,TANGENTS,NORMAL,TIME,variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                    & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*999)
                                  DO version_idx=1,DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%NUMBER_OF_VERSIONS
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(version_idx)
                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
                                      IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                        !If we are a boundary node then set the analytic value on the boundary
                                        CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                          & local_ny,BOUNDARY_CONDITION_FIXED,INITIAL_VALUE,ERR,ERROR,*999)
                                      ELSE
                                        !Set the initial condition.
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                          & FIELD_VALUES_SET_TYPE,local_ny,INITIAL_VALUE,ERR,ERROR,*999)
                                      ENDIF
                                    ENDIF
                                  ENDDO !version_idx
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
              ENDIF
            ENDIF
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

    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE")
    RETURN 1
    
  END SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_CALCULATE


  !
  !================================================================================================================================
  !
  !>Evaluate the analytic solutions for a diffusion equation
  SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SUBTYPE,ANALYTIC_FUNCTION_TYPE,X, &
    & TANGENTS,NORMAL,TIME,VARIABLE_TYPE,GLOBAL_DERIVATIVE,COMPONENT_NUMBER,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SUBTYPE !<The subtype of equation to evaluate
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: X(:) !<X(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: TANGENTS(:,:) !<TANGENTS(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: NORMAL(:) !<NORMAL(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: TIME !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: GLOBAL_DERIVATIVE !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: ANALYTIC_PARAMETERS(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: MATERIALS_PARAMETERS(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the analytic function value.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: k,phi,A,B,C,D,A1,A2,A3,A4
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,K_PARAM,L_PARAM,CONST_PARAM,BETA_PARAM,LAMBDA_PARAM,MU_PARAM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    !These are parameters for the analytical solution
!     k = 1.0_DP !this is a time decay constant for the exponential term
!     phi = 0.785398163397_DP !pi/4 - this sets the orientation of the solution relative to the axes
    k = 10.0_DP !this is a time decay constant for the exponential term
    phi = 1.0_DP !pi/4 - this sets the orientation of the solution relative to the axes
 
    !Solution parameters for 
    A1 = 0.4_DP
    A2 = 0.3_DP
    A3 = 0.2_DP
    A4 = 0.1_DP

    CALL ENTERS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR,*999)

    SELECT CASE(EQUATIONS_SUBTYPE)
    CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
        !For \del u/\del t = a \del^2 u/\del x^2
        !u(x,t)=A.exp(-4.\mu^2.t)cos(\mu.x+B)+C
        !see http://eqworld.ipmnet.ru/en/solutions/lpde/lpde101.pdf
        !OpenCMISS has \del u/\del t + k \del^2 u/\del x^2 = 0, thereform with \mu=2.\pi/L we have
        !u(x,t)=A.exp(4.\pi^2.k.t/L^2)cos(2.\pi.x/L+B)+C
        K_PARAM=MATERIALS_PARAMETERS(1)
        A_PARAM=ANALYTIC_PARAMETERS(1)
        B_PARAM=ANALYTIC_PARAMETERS(2)
        C_PARAM=ANALYTIC_PARAMETERS(3)
        L_PARAM=ANALYTIC_PARAMETERS(4)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A_PARAM*EXP(4.0_DP*PI**2*K_PARAM*TIME/L_PARAM**2)*COS(2.0_DP*PI*X(1)/L_PARAM+B_PARAM)+C_PARAM
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
        !u=exp(-kt)*sin(sqrt(k)*(x*cos(phi)+y*sin(phi)))
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=EXP(-k*TIME)*SIN((SQRT(k))*(X(1)*COS(phi)+X(2)*SIN(phi)))!Need to specify time, k and phi!
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
        !u=A1*exp(-t)*(x^2+y^2+z^2)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A1*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1)
        !u=exp(-kt)*sin(sqrt(k)*(x*cos(phi)+y*sin(phi)))
        !These are parameters for the 3D analytical solution with a linear source
        A = -0.25_DP
        B = 0.5_DP   
        C = 0.5_DP
        D = 0.5_DP
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=EXP(A*TIME)*EXP(B*X(1))*EXP(C*X(2))*EXP(D*X(3))!Need to specify time, k and phi!
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
     SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
        !For del u/del t = del^2 u/del x^2 + a + bu + cu^m with a = 0 then
        !u(x,t) = [+/-\beta + C.exp(\lamba.t+/\mu.x)]^(2/1-m) where
        !\beta=\sqrt(-c/b); \lamba=b(1-m)(m+3)/(2(m+1)); \mu = \sqrt((b(1-m)^2)/(2.(m+1))
        !see http://eqworld.ipmnet.ru/en/solutions/npde/npde1104.pdf
        A_PARAM=MATERIALS_PARAMETERS(1)
        B_PARAM=MATERIALS_PARAMETERS(2)
        C_PARAM=MATERIALS_PARAMETERS(3)
        BETA_PARAM=SQRT(-C_PARAM/B_PARAM)
        LAMBDA_PARAM=-5.0_DP*B_PARAM/6.0_DP
        MU_PARAM=SQRT(B_PARAM/6.0_DP)
        CONST_PARAM=1.0_DP
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=1.0_DP/(BETA_PARAM+CONST_PARAM*EXP(LAMBDA_PARAM*TIME+MU_PARAM*X(1)))**2
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP                                 
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
        !For del u/del t = del^2 u/del x^2 + a + b.e^c.u
        !u(x,t) = -2/c.ln[+/-\beta + C.exp(+/-\mu.x-a.c.t/2)] where
        !\beta=\sqrt(-b/a); \mu = \sqrt(a.c/2)
        !see http://eqworld.ipmnet.ru/en/solutions/npde/npde1105.pdf
        A_PARAM=MATERIALS_PARAMETERS(1)
        B_PARAM=MATERIALS_PARAMETERS(2)
        C_PARAM=MATERIALS_PARAMETERS(3)
        CONST_PARAM=1.0_DP
        BETA_PARAM=SQRT(-B_PARAM/A_PARAM)
        MU_PARAM=SQRT(A_PARAM*C_PARAM/2.0_DP)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=-2.0_DP/C_PARAM*LOG(BETA_PARAM+CONST_PARAM*EXP(MU_PARAM*X(1)-A_PARAM*C_PARAM*TIME/2.0_DP))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP 
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A1*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_V_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A2*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELVDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A1*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_V_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A2*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELVDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A1*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_V_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A2*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELVDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_U1_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A3*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELU1DELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_U2_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A4*EXP(-1.0_DP*TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELU2DELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
          CASE(GLOBAL_DERIV_S1)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SUBTYPE,"*",ERR,ERROR))// &
        & " is invalid."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
     
    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equation type of a classical field equations set class.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a diffusion equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a diffusion equation type of an classical field equations set class.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,& 
         & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE,EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE,EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)        
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
          & " is not valid for a diffusion equation type of an classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a diffusion equation type of a classical field equations set class.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the linear diffusion equation.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
      & NUMBER_OF_DIMENSIONS, NUMBER_OF_MATERIALS_COMPONENTS, NUMBER_OF_SOURCE_COMPONENTS,imy_matrix,Ncompartments, &
      & GEOMETRIC_COMPONENT_NUMBER
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SOURCE
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,EQUATIONS_SET_FIELD_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: num_var,num_var_count,NUMBER_OF_MATERIALS_COUPLING_COMPONENTS    
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS    
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:),VARIABLE_U_TYPES(:),COUPLING_MATRIX_STORAGE_TYPE(:), &
      & COUPLING_MATRIX_STRUCTURE_TYPE(:)

    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
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
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Equations",ERR,ERROR,*999)
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
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
            !do nothing 
          CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
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
              !Do nothing
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a linear diffusion equation."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
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
            SELECT CASE(EQUATIONS_SET%SUBTYPE)
            CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
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
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & "del U/del n",ERR,ERROR,*999)
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
                SELECT CASE(EQUATIONS_SET%SUBTYPE)
                CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                  !Check the field created by advection-diffusion routines for the coupled problem
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, & 
                    & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
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
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                    & ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    component_idx=1
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, & 
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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

                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                  !uses number of compartments to check that appropriate number and type of variables have been set on the
                  !dependent field
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
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
                  & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
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
                    DO component_idx=1,NUMBER_OF_DIMENSIONS
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                        & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                        & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    ENDDO !component_idx
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
            END SELECT
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation"
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
                !Create the auto created materials field
                SELECT CASE(EQUATIONS_SET%SUBTYPE)
                CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
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
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "Materials",ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                  ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                    !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                  ELSE
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2
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
                  !Default the source materials components to the first component geometric interpolation with constant
                  !interpolation
                  IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                    ENDDO !components_idx
                  ENDIF
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
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
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "Materials",ERR,ERROR,*999)
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
                END SELECT
              ELSE
                !Check the user specified field
                SELECT CASE(EQUATIONS_SET%SUBTYPE)
                CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)

                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                      & ERR,ERROR,*999)
                  ELSE
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                      & ERR,ERROR,*999)
                  ENDIF
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
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
                END SELECT
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
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS             
                ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                  !Constant source. Materials field components are 1 for each dimension and 1 for the constant source
                  !i.e., k and c in div(k.grad(u(x)))=c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  DO component_idx=1,Ncompartments
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDIF
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
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
              & " is invalid for a linear diffusion equation."
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
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Source",ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                NUMBER_OF_SOURCE_COMPONENTS=1
                !Set the number of source components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_SOURCE_COMPONENTS,ERR,ERROR,*999)
                !Default the source components to the geometric interpolation setup with constant interpolation
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
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
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
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
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c  T y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                  IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                    EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                    IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                      IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                          & ERR,ERROR,*999)
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)=0.0_DP
                        SELECT CASE(EQUATIONS_SET%SUBTYPE)
                        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
                            IF(NUMBER_OF_DIMENSIONS/=1) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a no source diffusion equation requires that there be 1 geometric dimension."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=4
                            !Set analytic function type
                            EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1
                          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
                            !Check that domain is 2D
                            IF(NUMBER_OF_DIMENSIONS/=2) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a no source diffusion equation requires that there be 2 geometric dimensions."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a no source diffusion equation."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a constant source diffusion equation requires that there be 3 geometric dimensions."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a constant source diffusion equation."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a linear source diffusion equation requires that there be 3 geometric dimensions."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a linear source diffusion equation."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM)
                            !Check that domain is 2D
                            IF(NUMBER_OF_DIMENSIONS/=2) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a multi-compartment diffusion equation requires that there be 2 geometric dimensions."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= &
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a multi-compartment diffusion equation requires that there be 3 geometric dimensions."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a multi-compartment diffusion equation requires that there be 3 geometric dimensions."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_FOUR_COMP_THREE_DIM)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " for a multi-compartment diffusion requires that there be 3 geometric dimensions."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_FOUR_COMP_THREE_DIM
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a multi-compartment diffusion equation."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        CASE DEFAULT
                          LOCAL_ERROR="The equation set subtype of "// &
                            & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                            & " is invalid for an analytical diffusion equation."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                        !Create analytic field if required
                        IF(NUMBER_OF_ANALYTIC_COMPONENTS>=1) THEN
                          IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                            !Create the auto created source field
                            CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                              & EQUATIONS_ANALYTIC%ANALYTIC_FIELD,ERR,ERROR,*999)
                            CALL FIELD_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,"Analytic Field",ERR,ERROR,*999)
                            CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                            CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_INDEPENDENT_TYPE, &
                              & ERR,ERROR,*999)
                            CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                              & ERR,ERROR,*999)
                            CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                              & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                            CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,EQUATIONS_SET%GEOMETRY% &
                              & GEOMETRIC_FIELD,ERR,ERROR,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,1,ERR,ERROR,*999)
                            CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,[FIELD_U_VARIABLE_TYPE], &
                              & ERR,ERROR,*999)
                            CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & "Analytic",ERR,ERROR,*999)
                            CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                            CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_DP_TYPE,ERR,ERROR,*999)
                            !Set the number of analytic components
                            CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,ERR,ERROR,*999)
                            !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
                            CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                            DO component_idx=1,NUMBER_OF_ANALYTIC_COMPONENTS
                              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                              CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            ENDDO !component_idx
                            !Default the field scaling to that of the geometric field
                            CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & ERR,ERROR,*999)
                            CALL FIELD_SCALING_TYPE_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                          ELSE
                            !Check the user specified field
                            CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                            CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                            CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                            IF(NUMBER_OF_ANALYTIC_COMPONENTS==1) THEN
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                            ELSE
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                            ENDIF
                            CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                              & ERR,ERROR,*999)
                            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,ERR,ERROR,*999)
                          ENDIF
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations set materials has not been finished.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  !Finish creating the analytic field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,ERR,ERROR,*999)
                  !Set the default values for the analytic field
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET%SUBTYPE)
                  CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
                    SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
                      !Set A
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1.0_DP,ERR,ERROR,*999)
                      !Set B
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,2,1.0_DP,ERR,ERROR,*999)
                      !Set C
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,3,1.0_DP,ERR,ERROR,*999)
                      !Set L
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,4,1.0_DP,ERR,ERROR,*999)                      
                    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
                      !Do nothing
                    CASE DEFAULT
                      LOCAL_ERROR="The specified analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for a no source diffusion equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                     END SELECT
                  CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE DEFAULT
                    LOCAL_ERROR="The equation set subtype of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                      & " is invalid for an analytical linear diffusion equation."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
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
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SUBTYPE)      
              CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_V_VARIABLE_TYPE,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELVDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                imy_matrix = EQUATIONS_SET_FIELD_DATA(1)
                Ncompartments = EQUATIONS_SET_FIELD_DATA(2)    
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
              CASE DEFAULT
                CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              END SELECT
              IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN              
                CALL EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              ENDIF
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              !Set up matrix storage and structure
              IF(EQUATIONS%LUMPING_TYPE==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET(EQUATIONS_MATRICES, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                  [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],ERR,ERROR,*999)
              ELSE
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                  IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)THEN
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not a linear diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP")
    RETURN 1
    
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the non-linear diffusion equation.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
      & NUMBER_OF_DIMENSIONS,NUMBER_OF_MATERIALS_COMPONENTS
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE.OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)            
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
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & "U",ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & "del U/del n",ERR,ERROR,*999)
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
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, & 
                & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
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
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
              & " is invalid for a nonlinear diffusion equation"
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
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Materials",ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) THEN
                  !Quadratic source. Materials field components are 1 for each dimension and 3 for the quadratic source
                  !i.e., k and a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)u(x)+c(x)u^2(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ELSE
                  !Exponential source. Matierals field components are 1 for each dimension and 3 for the exponential source
                  !i.e., k, a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)e^[c(x)u(x)]
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
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
                !Default the source materials components to the first component geometric interpolation with constant
                !interpolation                
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
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
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2                                
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
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
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) THEN
                  !Quadratic source. Materials field components are 1 for each dimension and 3 for the quadratic source
                  !i.e., k and a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)u(x)+c(x)u^2(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ELSE
                  !Exponential source. Matierals field components are 1 for each dimension and 3 for the exponential source
                  !i.e., k, a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)e^[c(x)u(x)]
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                !Set the source values to 1.0
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing put the constant source directly into the RHS
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing put the constant source directly into the RHS
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                  IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                    IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                      GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                      IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                          & ERR,ERROR,*999)
                        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) THEN
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
                            !Check that domain is 1D
                            IF(NUMBER_OF_DIMENSIONS/=1) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " requires that there be 1 geometric dimension."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 3,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            !Check that the a parameter is zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,1,A_PARAM,ERR,ERROR,*999)
                            IF(ABS(A_PARAM)>ZERO_TOLERANCE)  &
                              & CALL FLAG_ERROR("The 1st material component must be zero.",ERR,ERROR,*999)
                            !Check that the b parameter is not zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,2,B_PARAM,ERR,ERROR,*999)
                            IF(B_PARAM<ZERO_TOLERANCE)  &
                              & CALL FLAG_ERROR("The 2nd material component must be greater than zero.",ERR,ERROR,*999)
                            !Check to ensure we get real solutions
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,2,B_PARAM,ERR,ERROR,*999)
                            IF((B_PARAM*C_PARAM)>ZERO_TOLERANCE) &
                              & CALL FLAG_ERROR("The product of the 2nd and 3rd material components must not be positive.", &
                              & ERR,ERROR,*999)
                            !Set the number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=1
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= &
                              & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a nonlinear diffusion equation with a quadratic source."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        ELSE
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
                            !Check that domain is 1D
                            IF(NUMBER_OF_DIMENSIONS/=1) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " requires that there be 1 geometric dimension."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 3,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            !Check that the a parameter is not zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,1,A_PARAM,ERR,ERROR,*999)
                            IF(ABS(A_PARAM)<ZERO_TOLERANCE)  &
                              & CALL FLAG_ERROR("The 1st material component must not be zero.",ERR,ERROR,*999)
                            !Check that the c parameter is not zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,3,C_PARAM,ERR,ERROR,*999)
                            IF(ABS(C_PARAM)<ZERO_TOLERANCE)  &
                              & CALL FLAG_ERROR("The 3rd material component must not be zero.",ERR,ERROR,*999)
                            !Check to ensure we get real solutions
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,2,B_PARAM,ERR,ERROR,*999)
                            IF((A_PARAM*B_PARAM)>ZERO_TOLERANCE) &
                              & CALL FLAG_ERROR("The product of the 1st and 2nd material components must not be positive.", &
                              & ERR,ERROR,*999)
                            IF((A_PARAM*C_PARAM)<ZERO_TOLERANCE) &
                              & CALL FLAG_ERROR("The product of the 1st and 3rd material components must not be negative.", &
                              & ERR,ERROR,*999)
                            !Set the number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= &
                              & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a nonlinear diffusion equation with an exponential source."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        ENDIF
                        !Create analytic field if required
                        IF(NUMBER_OF_ANALYTIC_COMPONENTS>=1) THEN
                          IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                            !Create the auto created source field
                            CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                              & EQUATIONS_ANALYTIC%ANALYTIC_FIELD,ERR,ERROR,*999)
                            CALL FIELD_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,"Analytic Field",ERR,ERROR,*999)
                            CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                            CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_INDEPENDENT_TYPE, &
                              & ERR,ERROR,*999)
                            CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                              & ERR,ERROR,*999)
                            CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                              & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                            CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,EQUATIONS_SET%GEOMETRY% &
                              & GEOMETRIC_FIELD,ERR,ERROR,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,1,ERR,ERROR,*999)
                            CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,[FIELD_U_VARIABLE_TYPE], &
                              & ERR,ERROR,*999)
                            CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & "Analytic",ERR,ERROR,*999)
                            CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                            CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_DP_TYPE,ERR,ERROR,*999)
                            !Set the number of analytic components
                            CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,ERR,ERROR,*999)
                            !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
                            CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                            DO component_idx=1,NUMBER_OF_ANALYTIC_COMPONENTS
                              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                              CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            ENDDO !component_idx
                            !Default the field scaling to that of the geometric field
                            CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & ERR,ERROR,*999)
                            CALL FIELD_SCALING_TYPE_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                          ELSE
                            !Check the user specified field
                            CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                            CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                            CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                            IF(NUMBER_OF_ANALYTIC_COMPONENTS==1) THEN
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                            ELSE
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                            ENDIF
                            CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                              & ERR,ERROR,*999)
                            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,ERR,ERROR,*999)
                          ENDIF
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations set materials is not finished.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  !Finish creating the analytic field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,ERR,ERROR,*999)
                  !Set the default values for the analytic field
                  SELECT CASE(EQUATIONS_SET%SUBTYPE)
                  CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE DEFAULT
                    LOCAL_ERROR="The equation set subtype of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                      & " is invalid for an analytical nonlinear diffusion equation."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
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
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              ! Use the analytic Jacobian calculation
              CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                & ERR,ERROR,*999)
              !Set up matrix storage and structure
              IF(EQUATIONS%LUMPING_TYPE==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET(EQUATIONS_MATRICES, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],ERR,ERROR,*999)
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                    & ERR,ERROR,*999)
                 CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                    EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
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
              & " is invalid for a nonlinear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a nonlinear diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not a nonlinear diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
    RETURN 1
    
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion problem.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a diffusion equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the diffusion problem pre-solve.
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local variables
    LOGICAL :: UPDATE_MATERIALS
    LOGICAL :: UPDATE_BOUNDARY_CONDITIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    UPDATE_MATERIALS = .FALSE.
    UPDATE_BOUNDARY_CONDITIONS = .TRUE.

    IF( UPDATE_MATERIALS ) THEN
      !CALL DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_MATERIALS_FIELD(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
    ENDIF
    
    !     IF( UPDATE_BOUNDARY_CONDITIONS ) THEN
    !       CALL DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
    !     ENDIF

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
            & PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
            ! do nothing ???
            CALL DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_ANALYTIC_VALUES(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
          CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
            & PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE diffusion pre solve... ",ERR,ERROR,*999)
            IF(SOLVER%DYNAMIC_SOLVER%ALE) THEN
              !First update mesh and calculate boundary velocity values
              CALL DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              !Then apply both normal and moving mesh boundary conditions
              !CALL DIFFUSION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            ELSE  
              CALL FLAG_ERROR("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a diffusion equation type of a classical field problem class."
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

    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !
  !>Within the diffusion pre-solve, update the boundary conditions
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
!     TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
! !    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field
!     TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
!     TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
!     TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
!     TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
!     TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
!     TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
!     TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
! !    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
!     TYPE(VARYING_STRING) :: LOCAL_ERROR
! !    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
! !    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
! !    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
!     REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
!     INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE
! 
!     REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
!     REAL(DP) :: VALUE,X(3) !<The value to add
! !     REAL(DP) :: k_xx, k_yy, k_zz
!     INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,variable_idx
!     INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
!     INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE
!     INTEGER(INTG) :: GLOBAL_DERIV_INDEX
! !    INTEGER(INTG) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
! !    INTEGER(INTG) :: DERIVATIVE_NUMBER !<The node derivative number
! !    INTEGER(INTG) :: COMPONENT_NUMBER !<The field variable component number
! !    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of (geometry) nodes
! !    INTEGER(INTG) :: LOCAL_NODE_NUMBER
! !    INTEGER(INTG) :: EQUATIONS_SET_IDX
! !    INTEGER(INTG) :: equations_row_number

    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",ERR,ERROR,*999)

    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)!This routine previously set analytic BCs, but this has been moved. Needs rewriting to set
    !boundary conditions from file, time varying if appropriate.

!     IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
!        !write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!        !write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
!       IF(ASSOCIATED(SOLVER)) THEN
!         IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
!           SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
!             CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
!                 SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
!                 IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
!                   SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
!                   EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
!                   IF(ASSOCIATED(EQUATIONS)) THEN
!                     EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
!                     IF(ASSOCIATED(EQUATIONS_SET)) THEN
!                       IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
!                         DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
!                         IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
!                           GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
!                           IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
!                             CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
!                               & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
!                             NULLIFY(GEOMETRIC_VARIABLE)
!                             NULLIFY(GEOMETRIC_PARAMETERS)
!                             CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
!                             CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,& 
!                               & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
!                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
!                               variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
!                               FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
!                               IF(ASSOCIATED(FIELD_VARIABLE)) THEN
!                                 DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!                                   IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
!                                     & FIELD_NODE_BASED_INTERPOLATION) THEN
!                                     DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
!                                     IF(ASSOCIATED(DOMAIN)) THEN
!                                       IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
!                                         DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
!                                         IF(ASSOCIATED(DOMAIN_NODES)) THEN
!                                           !Loop over the local nodes excluding the ghosts.
!                                           DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!                                             !!TODO \todo We should interpolate the geometric field here and the node position.
!                                             DO dim_idx=1,NUMBER_OF_DIMENSIONS
!                                               local_ny= & 
!                                           & GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
!                                               X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
!                                             ENDDO !dim_idx
!                                             !Loop over the derivatives
!                                             DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
!                                               ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
!                                               GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx)
!                                               CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X, & 
!                                                 & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX, &
!                                                 & ANALYTIC_FUNCTION_TYPE,ERR,ERROR,*999)
!                                               local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
!                                                 & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
!                                               CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
!                                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
!                                               BOUNDARY_CONDITION_CHECK_VARIABLE=SOLVER_EQUATIONS%BOUNDARY_CONDITIONS% &
!                                                 & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% & 
!                                                 & CONDITION_TYPES(local_ny)
!                                               IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
!                                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, & 
!                                                  & variable_type,FIELD_VALUES_SET_TYPE,local_ny, & 
!                                                  & VALUE,ERR,ERROR,*999)
!                                               ENDIF
! !                                              IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
! !                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
!                                                   !If we are a boundary node then set the analytic value on the boundary
! !                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
! !                                                    & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
! !                                                ENDIF
! !                                              ENDIF
!                                             ENDDO !deriv_idx
!                                           ENDDO !node_idx
!                                         ELSE
!                                           CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
!                                         ENDIF
!                                       ELSE
!                                         CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
!                                       ENDIF
!                                     ELSE
!                                       CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
!                                     ENDIF
!                                   ELSE
!                                     CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
!                                   ENDIF
!                                 ENDDO !component_idx
!                                 CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
!                                 CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
!                                 CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                                 CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                               ELSE
!                                 CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
!                               ENDIF
! 
!                              ENDDO !variable_idx
!                              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
!                               & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
!                           ELSE
!                             CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
!                           ENDIF            
!                         ELSE
!                           CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
!                         ENDIF
!                       ELSE
!                         !CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
!                       ENDIF
!                     ELSE
!                       CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
!                     ENDIF
!                   ELSE
!                     CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
!                   END IF                
!                 ELSE
!                   CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
!                 END IF  
!                 CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                 CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!             !do nothing?! 
!             CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
!             !do nothing?! 
!             CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
!                 SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
!                 IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
!                   SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
!                   EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
!                   IF(ASSOCIATED(EQUATIONS)) THEN
!                     EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
!                     IF(ASSOCIATED(EQUATIONS_SET)) THEN
!                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
!                         DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
!                         IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
!                           GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
!                           IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
!                             CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
!                               & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
!                             NULLIFY(GEOMETRIC_VARIABLE)
!                             NULLIFY(GEOMETRIC_PARAMETERS)
!                             CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
!                             CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,& 
!                               & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
!                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
!                               variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
!                               FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
!                               IF(ASSOCIATED(FIELD_VARIABLE)) THEN
!                                 DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!                                   IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
!                                     & FIELD_NODE_BASED_INTERPOLATION) THEN
!                                     DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
!                                     IF(ASSOCIATED(DOMAIN)) THEN
!                                       IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
!                                         DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
!                                         IF(ASSOCIATED(DOMAIN_NODES)) THEN
!                                           !Loop over the local nodes excluding the ghosts.
!                                           DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!                                             !!TODO \todo We should interpolate the geometric field here and the node position.
!                                             DO dim_idx=1,NUMBER_OF_DIMENSIONS
!                                               local_ny= & 
!                                           & GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
!                                               X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
!                                             ENDDO !dim_idx
!                                             !Loop over the derivatives
!                                             DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
!                                               ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
!                                               GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx)
!                                               CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X, & 
!                                                 & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX, &
!                                                 & ANALYTIC_FUNCTION_TYPE,ERR,ERROR,*999)
!                                               local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
!                                                 & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
!                                               CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
!                                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
!                                               BOUNDARY_CONDITION_CHECK_VARIABLE=SOLVER_EQUATIONS%BOUNDARY_CONDITIONS% &
!                                                 & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% & 
!                                                 & CONDITION_TYPES(local_ny)
!                                               IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
!                                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, & 
!                                                  & variable_type,FIELD_VALUES_SET_TYPE,local_ny, & 
!                                                  & VALUE,ERR,ERROR,*999)
!                                               ENDIF
! !                                              IF(variable_type==FIELD_U_VARIABLE_TYPE .OR. variable_type==FIELD_V_VARIABLE_TYPE) THEN
! !                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
! !                                                   If we are a boundary node then set the analytic value on the boundary
! !                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
! !                                                    & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
! !                                                ENDIF
! !                                              ENDIF
!                                             ENDDO !deriv_idx
!                                           ENDDO !node_idx
!                                         ELSE
!                                           CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
!                                         ENDIF
!                                       ELSE
!                                         CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
!                                       ENDIF
!                                     ELSE
!                                       CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
!                                     ENDIF
!                                   ELSE
!                                     CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
!                                   ENDIF
!                                 ENDDO !component_idx
!                                 CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
!                                 CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
!                                 CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                                 CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                               ELSE
!                                 CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
!                               ENDIF
! 
!                              ENDDO !variable_idx
!                              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
!                               & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
!                           ELSE
!                             CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
!                           ENDIF            
!                         ELSE
!                           CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
!                         ENDIF
!                       ELSE
!                         !CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
!                       ENDIF
!                     ELSE
!                       CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
!                     ENDIF
!                   ELSE
!                     CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
!                   END IF                
!                 ELSE
!                   CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
!                 END IF  
!                 CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!                 CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
!             CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
!               ! do nothing ???
!             CASE(PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
!               ! do nothing ???
!             CASE(PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
!               ! do nothing ???
!             CASE DEFAULT
!               LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
!                 & " is not valid for a diffusion equation type of a classical field problem class."
!             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
!           END SELECT
!         ELSE
!           CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
!         ENDIF
!       ELSE
!         CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
!       ENDIF
!     ELSE
!       CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
!     ENDIF
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS

  !   
  !================================================================================================================================
  !
  !updates the boundary conditions and source term to the required analytic values
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_ANALYTIC_VALUES(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

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
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SOURCE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    !    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
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
    REAL(DP) :: A1,D1
!    INTEGER(INTG) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
!    INTEGER(INTG) :: DERIVATIVE_NUMBER !<The node derivative number
!    INTEGER(INTG) :: COMPONENT_NUMBER !<The field variable component number
!    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of (geometry) nodes
!    INTEGER(INTG) :: LOCAL_NODE_NUMBER
!    INTEGER(INTG) :: EQUATIONS_SET_IDX
!    INTEGER(INTG) :: equations_row_number

    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_ANALYTIC_VALUES",ERR,ERROR,*999)

    A1 = 0.4_DP
    D1=1.0_DP

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
            & PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
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
                          EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)=CURRENT_TIME
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
                          !                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                          variable_type=DEPENDENT_FIELD%VARIABLES(2*eqnset_idx-1)%VARIABLE_TYPE
                          FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                            BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
                                    & FIELD_NODE_BASED_INTERPOLATION) THEN
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
                                              !If we are a boundary node then set the analytic value on the boundary
                                              !                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
                                              !                                                    & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                              !                                                ENDIF
                                              !                                              ENDIF
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
                              ELSE
                                CALL FLAG_ERROR("Boundary conditions variable is not associated",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Boundary conditions are not associated",ERR,ERROR,*999)
                            ENDIF
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
                IF(CONTROL_LOOP%PROBLEM%SUBTYPE==PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)THEN
                  !>Set the source field to a specified analytical function
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                      EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
                      IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
                        SOURCE_FIELD=>EQUATIONS_SOURCE%SOURCE_FIELD
                        IF(ASSOCIATED(SOURCE_FIELD)) THEN
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                              & ERR,ERROR,*999)
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
                                            local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                              & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                            X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                          ENDDO !dim_idx
                                          !Loop over the derivatives
                                          DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                            SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                            CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
                                              VALUE_SOURCE=-1*A1*EXP(-1*CURRENT_TIME)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)+6)
                                            CASE DEFAULT
                                              LOCAL_ERROR="The analytic function type of "// &
                                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*", &
                                                & ERR,ERROR))//" is invalid."
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
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !eqnset_idx
            ELSE
              CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a diffusion equation type of a classical field problem class."
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
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_ANALYTIC_VALUES")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_ANALYTIC_VALUES",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_ANALYTIC_VALUES")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_ANALYTIC_VALUES

  !
  !================================================================================================================================
  !
  !>Update mesh position and velocity for ALE diffusion problem
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

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

    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH",ERR,ERROR,*999)

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
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      SELECT CASE(EQUATIONS_SET%SUBTYPE)
                        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
                          ! do nothing ???
                        CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Diffusion update mesh ... ",ERR,ERROR,*999)
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
!                            CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, &
!                             & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CURRENT_TIME)

                            NULLIFY(MESH_DISPLACEMENT_VALUES)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                            IF(DIAGNOSTICS1) THEN
                              NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                                & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",3(X,E13.6))','3(3(X,E13.6))', &
                                & ERR,ERROR,*999)
                            ENDIF

!                            CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, &
!                             & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CURRENT_TIME)

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
                            & //TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                            & " is not valid for a diffusion equation type of a classical field problem class."
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
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
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

    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH
  !   
  !================================================================================================================================
  !
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION_ONE !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION_ONE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION_ONE !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_DIFFUSION_ONE !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_DIFFUSION_ONE !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
    INTEGER(INTG) :: I

    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_DIFFUSION_ONE)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion-one equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Store value diffusion-one dependent field at time, t ... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_DIFFUSION_ONE,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      DEPENDENT_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_ONE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DIFFUSION_ONE, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_DIFFUSION_ONE is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Diffusion-one equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Diffusion-one solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Diffusion-one solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_DIFFUSION_ONE, & 
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
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==2) THEN
                !--- Get the dependent field of the diffusion equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Store value of diffusion solution &
                   & (dependent field - V variable_type) at time, t ... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DIFFUSION_ONE,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      DEPENDENT_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_ONE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DIFFUSION_ONE, &
                          & FIELD_V_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_DIFFUSION_ONE is not associated.",ERR,ERROR,*999)
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

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_V_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,I,ERR,ERROR,*999)
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
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
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



    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION    
  !   
  !================================================================================================================================
  !
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION_ONE, SOLVER_DIFFUSION_TWO  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION_TWO, SOURCE_FIELD_DIFFUSION_ONE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION_ONE, SOLVER_EQUATIONS_DIFFUSION_TWO  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_DIFFUSION_ONE, SOLVER_MAPPING_DIFFUSION_TWO !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_DIFFUSION_ONE, EQUATIONS_SET_DIFFUSION_TWO !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE
    INTEGER(INTG) :: I


    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_DIFFUSION_ONE)
      NULLIFY(SOLVER_DIFFUSION_TWO)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion_two equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Update diffusion-one source field ... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DIFFUSION_TWO,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION_TWO=>SOLVER_DIFFUSION_TWO%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_TWO)) THEN
                  SOLVER_MAPPING_DIFFUSION_TWO=>SOLVER_EQUATIONS_DIFFUSION_TWO%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_TWO)) THEN
                    EQUATIONS_SET_DIFFUSION_TWO=>SOLVER_MAPPING_DIFFUSION_TWO%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_TWO)) THEN
                      DEPENDENT_FIELD_DIFFUSION_TWO=>EQUATIONS_SET_DIFFUSION_TWO%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_TWO)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DIFFUSION_TWO, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_DIFFUSION_TWO is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Diffusion-two equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Diffusion-two solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Diffusion-two solver equations are not associated.",ERR,ERROR,*999)
                END IF


                !--- Get the source field for the diffusion_one equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_DIFFUSION_ONE,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      SOURCE_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%SOURCE%SOURCE_FIELD
                      IF(ASSOCIATED(SOURCE_FIELD_DIFFUSION_ONE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(SOURCE_FIELD_DIFFUSION_ONE, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("SOURCE_FIELD_DIFFUSION_ONE is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Diffusion-one equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Diffusion-one solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Diffusion-one solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Copy the result from diffusion-two's dependent field to diffusion-one's source field
                IF(NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE==NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO) THEN
                  DO I=1,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_DIFFUSION_TWO, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,SOURCE_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,ERR,ERROR,*999)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIELDMESHDISPLACEMENTTYPE needs to be changed to appropriate type for this problem
                  END DO
                ELSE
                  LOCAL_ERROR="Number of components of diffusion-two dependent field "// &
                    & "is not consistent with diffusion-one-equation source field."
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
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
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



    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE
  !   
  !================================================================================================================================
  !
  !>Sets up the diffusion problem post solve.
  SUBROUTINE DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solver
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DIFFUSION_EQUATION_POST_SOLVE",ERR,ERROR,*999)
    
    NULLIFY(SOLVER2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
            & PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
            !CALL DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
          CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
            ! do nothing ???
          CASE(PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a diffusion type of a classical field problem class."
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

    CALL EXITS("DIFFUSION_EQUATION_POST_SOLVE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_POST_SOLVE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_POST_SOLVE")
    RETURN 1
    
  END SUBROUTINE DIFFUSION_EQUATION_POST_SOLVE
  
  !   
  !================================================================================================================================
  !
  
  !>Output data post solve
  SUBROUTINE DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

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

    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    CALL ENTERS("DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!       write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
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
  !          FILE="TRANSIENT_OUTPUT"
!!!!!!!!ADAPT THIS TO WORK WITH DIFFUSION AND NOT JUST FLUID MECHANICS
!                         METHOD="FORTRAN"
!                         EXPORT_FIELD=.TRUE.
!                         IF(EXPORT_FIELD) THEN          
!                           IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",ERR,ERROR,*999)
!                             CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
!                               & ERR,ERROR,*999)
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
!                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
!                           ENDIF
!                         ENDIF 

                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1 .OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1) THEN
                            CALL ANALYTIC_ANALYSIS_OUTPUT(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,ERR,ERROR,*999)
                          ENDIF
                        ENDIF
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
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
    CALL EXITS("DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !
  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,mh,mhs,ms,ng,nh,nhs,ni,nj,ns,my_compartment,Ncompartments,imatrix,num_var_count
    INTEGER(INTG) :: MESH_COMPONENT_1, MESH_COMPONENT_2
    REAL(DP) :: C_PARAM,K_PARAM,RWG,SUM,PGMJ(3),PGNJ(3),A_PARAM,COUPLING_PARAM,PGM,PGN
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1, DEPENDENT_BASIS_2
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: DAMPING_MATRIX,STIFFNESS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD,EQUATIONS_SET_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(FIELD_VARIABLE_PTR_TYPE) :: FIELD_VARIABLES(99)
    TYPE(EQUATIONS_MATRIX_PTR_TYPE) :: COUPLING_MATRICES(99) 
    INTEGER(INTG) :: FIELD_VAR_TYPES(99)
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME_1, QUADRATURE_SCHEME_2
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATION_PARAMETERS, &
      & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT, &
      & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
     
    CALL ENTERS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD
          ENDIF
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          ENDIF
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
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
          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & SOURCE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          ENDIF
          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN  
            ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATION_PARAMETERS=> &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT=> &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS=> &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PREVIOUS_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT=> &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR
          ENDIF
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
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
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
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
                        ENDDO !nj
                        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                        ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                          A_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS,NO_PART_DERIV)
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                        ENDIF
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
              IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                C_PARAM=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1, NO_PART_DERIV)
                mhs=0
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                      & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
              IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                !The value of the source term is +0.5*(C_1^{t}+C_1_{t+1}-C_2^{t}) 
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                  & ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT,ERR,ERROR,*999)
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                  & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,ERR,ERROR,*999)
                write(*,*) ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                write(*,*) DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                C_PARAM=0.5_DP*ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)- &
                  & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                !                     C_PARAM_1_T= EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1, NO_PART_DERIV)!<This is the value of the solution from the advection-diffusion equation at time T
                !                     C_PARAM_1_TPLUSONE= !<This is the value of the solution from the advection-diffusion equation at time T+deltaT
                !                     C_PARAM_2_T= EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(1, NO_PART_DERIV)!<This is the value of the solution from the diffusion equation at time T
                !                     C_PARAM=C_PARAM_1_T+C_PARAM_1_TPLUSONE+C_PARAM_2_T
                mhs=0
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                      & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP 
          ENDDO !ng

          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS=> &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PREVIOUS_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT=> &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                !The value of the source term is +0.5*(C_1^{t}+C_1_{t+1}-C_2^{t}) 
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                  & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,ERR,ERROR,*999)
                write(*,*) ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                C_PARAM=0.5_DP*ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                mhs=0
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                      & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh   
              ENDIF
              IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP 
            ENDDO !ng
          ENDIF

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
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                  IF(SOURCE_VECTOR%UPDATE_VECTOR) SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)= & 
                    & SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
!!!!!!!!!!!!!!!MULTI-COMPARTMENT DIFFUSION - PROTOTYPE FOR OTHER MULTI-COMPARTMENT MODELS IN FUTURE
!!!!!!!!!!!!!!!HAS BEEN SEPARATED HERE FOR EASE OF DEVELOPMENT & READABILITY OF THIS NEW FEATURE
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)

          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD
          EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD

          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
          DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING


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


          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
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
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & SOURCE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
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
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
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
                        ENDDO !nj
                        COUPLING_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR% &
                          & VALUES(my_compartment,NO_PART_DERIV)
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ & 
                          & SUM*RWG + QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG*COUPLING_PARAM
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
            IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
              C_PARAM=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1, NO_PART_DERIV)
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                    & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                ENDDO !ms
              ENDDO !mh
            ENDIF
            IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP 
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
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                IF(SOURCE_VECTOR%UPDATE_VECTOR) SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)= & 
                  & SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF

        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE,EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not calculate finite element stiffness matrices for a nonlinear source.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a diffusion equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a diffusion equation type.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE  
      CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE
      CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE
      CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE  
      CASE(PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE
      CASE(PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equations.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
    CALL ENTERS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM%SUBTYPE==PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM%SUBTYPE==PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM%SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
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
              & " is invalid for a linear diffusion equation."
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
              & " is invalid for a linear diffusion equation."
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
            CALL SOLVER_LABEL_SET(SOLVER,"Dynamic solver",ERR,ERROR,*999)
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
              & " is invalid for a linear diffusion equation."
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a linear diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the nonlinear diffusion problem
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
    CALL ENTERS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF( PROBLEM%SUBTYPE==PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
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
              & " is invalid for a nonlinear diffusion problem."
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
              & " is invalid for a nonlinear diffusion problem."
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
            CALL SOLVER_LABEL_SET(SOLVER,"Dynamic solver",ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
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
              & " is invalid for a nonlinear diffusion problem."
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
              & " is invalid for a nonlinear diffusion problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a nonlinear diffusion problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a nonlinear diffusion problem subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP
  
  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices for a diffusion equation finite element equations set.
  SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)
    
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
   
    CALL ENTERS("DIFFUSION_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a diffusion equation with no source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a diffusion equation with a constant source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a diffusion equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
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
              B_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
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
                      VALUE=-2.0_DP*B_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                        & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*U_VALUE
                      JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+VALUE*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh              
            ENDDO !ng
          ENDIF
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)                 
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
        CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for an ALE diffusion equation with no source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for an ALE diffusion equation with a constant source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for an ALE diffusion equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a multi component transport diffusion equation.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a diffusion equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
    
  END SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Diffusion equation finite element equations set.
  SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nj,nh,nhs,ni,ns
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,K_PARAM,RWG,SUM1,SUM2,PGMJ(3),PGNJ(3),U_VALUE
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: STIFFNESS_MATRIX,DAMPING_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("DIFFUSION_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a diffusion equation with no source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a diffusion equation with a constant source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a diffusion equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
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
            IF(STIFFNESS_MATRIX%FIRST_ASSEMBLY.OR.DAMPING_MATRIX%FIRST_ASSEMBLY.OR.RHS_VECTOR%FIRST_ASSEMBLY) THEN
              B_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                    nhs=0
                    !Loop over element columns
                    DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
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
                            K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                              & VALUES(nj,NO_PART_DERIV)
                            SUM1=SUM1+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          ENDDO !nj
                          SUM2=B_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)                        
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & (SUM1+SUM2)*RWG
                        ENDIF
                        IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                          DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                        ENDIF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                ENDDO !ms
              ENDDO !mh
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
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)                 
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
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
            IF(STIFFNESS_MATRIX%FIRST_ASSEMBLY.OR.DAMPING_MATRIX%FIRST_ASSEMBLY.OR.RHS_VECTOR%FIRST_ASSEMBLY) THEN
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                    nhs=0
                    !Loop over element columns
                    DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
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
                            K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                              & VALUES(nj,NO_PART_DERIV)
                            SUM1=SUM1+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          ENDDO !nj
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM1*RWG
                        ENDIF
                        IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                          DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                        ENDIF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                  IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
                ENDDO !ms
              ENDDO !mh
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
        CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for an ALE diffusion equation with no source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for an ALE diffusion equation with a constant source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for an ALE diffusion equation with a linear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not evaluate a residual for a multi component transport diffusion equation.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a diffusion equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE
 
  !
  !================================================================================================================================
  !

END MODULE DIFFUSION_EQUATION_ROUTINES
