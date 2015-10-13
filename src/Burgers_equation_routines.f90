!> \file
!> \author David Ladd
!> \brief This module handles all Burgers equation routines.
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
!> Contributor(s): Chris Bradley
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

!>This module handles all Burgers equation routines.
MODULE BURGERS_EQUATION_ROUTINES

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

#include "macros.h"

  IMPLICIT NONE

  PRIVATE
  
  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Burgers_AnalyticFunctionsEvaluate

  PUBLIC Burgers_BoundaryConditionsAnalyticCalculate

  PUBLIC BURGERS_EQUATION_EQUATIONS_SET_SETUP

  PUBLIC Burgers_EquationsSetSolutionMethodSet
  
  PUBLIC Burgers_EquationsSetSpecificationSet

  PUBLIC Burgers_FiniteElementJacobianEvaluate
  
  PUBLIC Burgers_FiniteElementResidualEvaluate
  
  PUBLIC Burgers_ProblemSpecificationSet

  PUBLIC BURGERS_EQUATION_PROBLEM_SETUP

  PUBLIC BURGERS_EQUATION_PRE_SOLVE,BURGERS_EQUATION_POST_SOLVE
 
CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !Calculates a one-dimensional dynamic solution to the burgers equation
  SUBROUTINE Burgers_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

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
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE
    !THESE ARE TEMPORARY VARIABLES - they need to be replace by constant field values and the current simulation time
    REAL(DP) :: TIME,NORMAL(3),TANGENTS(3,3)
    !CURRENT_TIME = 1.2_DP

    ENTERS("Burgers_BoundaryConditionsAnalyticCalculate",ERR,ERROR,*999)
    
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
            ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
            TIME=EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL Field_ParameterSetEnsureCreated(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
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
                                CALL Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE, &
                                  & X,TANGENTS,NORMAL,0.0_DP,variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                  & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,INITIAL_VALUE,ERR,ERROR,*999)
                                CALL Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE, &
                                  & X,TANGENTS,NORMAL,TIME,variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                  & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*999)
                                DO version_idx=1,DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%numberOfVersions
                                  local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                    & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(version_idx)
                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                  IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
                                    IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                      !If we are a boundary node then set the analytic value on the boundary
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
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
                            CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !component_idx
                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !variable_idx
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("Burgers_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Burgers_BoundaryConditionsAnalyticCalculate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Burgers_BoundaryConditionsAnalyticCalculate


  !
  !================================================================================================================================
  !
  !>Evaluate the analytic solutions for a Burgers equation
  SUBROUTINE Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE,X, &
    & TANGENTS,NORMAL,TIME,VARIABLE_TYPE,GLOBAL_DERIVATIVE,COMPONENT_NUMBER,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<The equations set to evaluate
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
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,D_PARAM,E_PARAM,X0_PARAM
    INTEGER(INTG) :: EQUATIONS_SUBTYPE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("Burgers_AnalyticFunctionsEvaluate",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) THEN
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ELSE
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      ELSE
        EQUATIONS_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
      END IF
    END IF
    SELECT CASE(EQUATIONS_SUBTYPE)
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1)
        !For del[u]/del[t] + u.(del[u]/del[x]) = nu.(del^2[u]/del[x]^2)
        !u(x,t)=1-tanh(x-x_0-t)/(2.nu))   with BCs,
        !u(0,t) = 2, u_{n} = 2.u_{n-1} - u_{n-2}
        !see http://www.cfd-online.com/Wiki/Burgers_equation
        !OpenCMISS has del[u]/del[t] + K.(del^2[u]/del[x]^2) + u.(del[u]/del[x]) = 0,
        !u(x,t)= 1 - tanh(x-x_0 - t)/(2.K)
        B_PARAM=MATERIALS_PARAMETERS(1)  !nu
        X0_PARAM=ANALYTIC_PARAMETERS(1)  !x_0
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=1.0_DP - TANH((X(1)-X0_PARAM-TIME)/(2.0_DP*B_PARAM)) 
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
          & " is invalid for a Burgers equation."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
      !a.del u/del t + b.del^2 u/del x^2 + c.u.del u/del x = 0
      A_PARAM=MATERIALS_PARAMETERS(1)
      B_PARAM=MATERIALS_PARAMETERS(2)
      C_PARAM=MATERIALS_PARAMETERS(3)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_1)
        !Analytic solution is u(x,t)=(D+a.x)/(E+c.t)
        D_PARAM = ANALYTIC_PARAMETERS(1)
        E_PARAM = ANALYTIC_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=(D_PARAM+A_PARAM*X(1))/(E_PARAM+C_PARAM*TIME)
          CASE(GLOBAL_DERIV_S1)
            VALUE=D_PARAM/(E_PARAM+C_PARAM*TIME)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            VALUE=0.0_DP
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2)
        !Analytic_solution=a.D+2.b/c(x-c.D.t+E)
        D_PARAM = ANALYTIC_PARAMETERS(1)
        E_PARAM = ANALYTIC_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A_PARAM*D_PARAM+2.0_DP*B_PARAM/(C_PARAM*(X(1)-C_PARAM*D_PARAM*TIME+E_PARAM))
          CASE(GLOBAL_DERIV_S1)
            VALUE=-2.0_DP*B_PARAM/(C_PARAM*(X(1)-C_PARAM*D_PARAM*TIME+E_PARAM)**2)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            VALUE=0.0_DP
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
          & " is invalid for a generalised Burgers equation."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SUBTYPE,"*",ERR,ERROR))// &
        & " is invalid."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    
    EXITS("Burgers_AnalyticFunctionsEvaluate")
    RETURN
999 ERRORSEXITS("Burgers_AnalyticFunctionsEvaluate",ERR,ERROR)
    RETURN 1
  END SUBROUTINE Burgers_AnalyticFunctionsEvaluate

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a burgers equation type of an fluid mechanics equations set class.
  SUBROUTINE Burgers_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Burgers_EquationsSetSolutionMethodSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
        & EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)     
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
          & " is not valid for a Burgers equation type of an classical field equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Burgers_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Burgers_EquationsSetSolutionMethodSet",ERR,ERROR)
    EXITS("Burgers_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Burgers_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Burgers type of a fluid mechanics equations set.
  SUBROUTINE Burgers_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Burgers_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers equation set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_BURGERS_SUBTYPE, &
          & EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE, &
          & EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
          & EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
          & " is not valid for a Burgers type of fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_BURGERS_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Burgers_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Burgers_EquationsSetSpecificationSet",err,error)
    EXITS("Burgers_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Burgers_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Burgers equation type of a fluid mechanics equations set class.
  SUBROUTINE BURGERS_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a BURGERS equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
      & NUMBER_OF_DEPENDENT_COMPONENTS,NUMBER_OF_GEOMETRIC_COMPONENTS,NUMBER_OF_MATERIALS_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BURGERS_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
        EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Burgers_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear burgers equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_GEOMETRIC_COMPONENTS,ERR,ERROR,*999)
                NUMBER_OF_DEPENDENT_COMPONENTS=NUMBER_OF_GEOMETRIC_COMPONENTS
              ELSE
                NUMBER_OF_DEPENDENT_COMPONENTS=1
              ENDIF
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DEPENDENT_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_DEPENDENT_COMPONENTS,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              DO component_idx=1,NUMBER_OF_DEPENDENT_COMPONENTS
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              ENDDO !component_idx
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_DEPENDENT_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, & 
                & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_GEOMETRIC_COMPONENTS,ERR,ERROR,*999)
                NUMBER_OF_DEPENDENT_COMPONENTS=NUMBER_OF_GEOMETRIC_COMPONENTS
              ELSE
                NUMBER_OF_DEPENDENT_COMPONENTS=1
              ENDIF
              IF(NUMBER_OF_DEPENDENT_COMPONENTS==1) THEN
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, & 
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
              ENDIF
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DEPENDENT_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & NUMBER_OF_DEPENDENT_COMPONENTS,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_DEPENDENT_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear Burgers equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                !Not an inviscid Burgers equation
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
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    !1 materials field component
                    !i.e., k = viscosity*(-1) in du/dt + k*(d^2u/dx^2)+ u*(du/dx) = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=1
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) 
                    !3 materials field components
                    !i.e., a.du/dt + b.(d^2u/dx^2) + c.u*(du/dx)  = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=3
                  CASE DEFAULT
                    LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                      & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a nonlinear Burgers equation."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                  !Set the number of materials components
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                  !Default the materials components to the 1st geometric component interpolation setup with constant interpolation
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  DO component_idx=1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
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
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    !1 materials field component
                    !i.e., k = viscosity*(-1) in du/dt + k*(d^2u/dx^2)+ u*(du/dx) = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=1
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                      & ERR,ERROR,*999)
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) 
                    !3 materials field components
                    !i.e., a.du/dt + b.(d^2u/dx^2) + c.u*(du/dx)  = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=3
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                      & ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                      & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a nonlinear Burgers equation."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                !Not an inviscid Burgers equation
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  !Set the default values for the materials field
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    !1 materials field component. Default to
                    !du/dt - d^2u/dx^2 + u*(du/dx) = 0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,-1.0_DP,ERR,ERROR,*999)
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) 
                    !3 materials field components. Default to
                    !du/dt - d^2u/dx^2 + u*(du/dx) = 0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,1.0_DP,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,2,-1.0_DP,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,3,1.0_DP,ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                      & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a nonlinear Burgers equation."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
             ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear Burgers equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing 
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing 
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear Burgers equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_GEOMETRIC_COMPONENTS,ERR,ERROR,*999)
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1)
                            !Check that domain is 1D
                            IF(NUMBER_OF_GEOMETRIC_COMPONENTS/=1) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GEOMETRIC_COMPONENTS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " requires that there be 1 geometric dimension."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1
                            NUMBER_OF_ANALYTIC_COMPONENTS=1
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a Burgers equation."
                            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_1, &
                            & EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2)
                            !Check that domain is 1D
                            IF(NUMBER_OF_GEOMETRIC_COMPONENTS/=1) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GEOMETRIC_COMPONENTS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " requires that there be 1 geometric dimension."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 3,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            NUMBER_OF_ANALYTIC_COMPONENTS=2
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for a generalised Burgers equation."
                            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                          CALL FlagError("Not implemented.",ERR,ERROR,*999)
                        CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
                          CALL FlagError("Not implemented.",ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The equation set subtype of "// &
                            & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
                            & " is invalid for an analytical nonlinear Burgers equation."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
                        CALL FlagError("Equations set materials is not finished.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set geometric field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations analytic is not associated.",ERR,ERROR,*999)
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
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
                    !Default the analytic parameter value to 0.0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
                    !Default the analytic parameter values to 1.0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,1.0_DP,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,2,1.0_DP,ERR,ERROR,*999)                    
                  CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The equation set subtype of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
                      & " is invalid for an analytical nonlinear Burgers equation."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear Burgers equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE (EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",ERR,ERROR,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                  CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.FALSE.,ERR,ERROR,*999)
                ELSE
                  CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
                ENDIF
                CALL EquationsMapping_ResidualVariableTypesSet(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                !Set up matrix storage and structure
                IF(EQUATIONS%LUMPING_TYPE==EQUATIONS_LUMPED_MATRICES) THEN
                  !Set up lumping
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                    CALL EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET(EQUATIONS_MATRICES, &
                      & [EQUATIONS_MATRIX_LUMPED],ERR,ERROR,*999)
                  ELSE
                    CALL EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET(EQUATIONS_MATRICES, &
                      & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],ERR,ERROR,*999)
                  ENDIF
                  SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                  CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],ERR,ERROR,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                        [EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],ERR,ERROR,*999)
                    ELSE
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],ERR,ERROR,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                        [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],ERR,ERROR,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                      & ERR,ERROR,*999)                    
                  CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],ERR,ERROR,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                        & [EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],ERR,ERROR,*999)
                    ELSE
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],ERR,ERROR,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                        & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],ERR,ERROR,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES, &
                      & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                    CALL EquationsMatrices_NonlinearStructureTypeSet(EQUATIONS_MATRICES,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                      & ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The equations matrices sparsity type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                  CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],ERR,ERROR,*999)
                    ELSE
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],ERR,ERROR,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                      & ERR,ERROR,*999)
                  CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                        & ERR,ERROR,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                        [EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                    ELSE
                      CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                        & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                        & ERR,ERROR,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                        [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES, &
                      & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                    CALL EquationsMatrices_NonlinearStructureTypeSet(EQUATIONS_MATRICES, &
                      EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The equations matrices sparsity type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
                ! Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a nonlinear Burgers equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",ERR,ERROR,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES,MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, & 
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                ! Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a nonlinear Burgers equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is invalid for an analytical nonlinear Burgers equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a nonlinear Burgers equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not a nonlinear Burgers equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BURGERS_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BURGERS_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the BURGERS problem pre-solve.
  SUBROUTINE BURGERS_EQUATION_PRE_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BURGERS_EQUATION_PRE_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification array is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
              DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
              IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
                IF(DYNAMIC_SOLVER%SOLVER_INITIALISED) &
                  & CALL Burgers_PreSolveUpdateAnalyticValues(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Solver dynamic solver is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a Burgers equation type of a fluid mechanics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("BURGERS_EQUATION_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_PRE_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BURGERS_EQUATION_PRE_SOLVE
      

  !   
  !================================================================================================================================
  !
  !updates the boundary conditions and source term to the required analytic values
  SUBROUTINE Burgers_PreSolveUpdateAnalyticValues(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:),GEOMETRIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP) :: NORMAL(3),TANGENTS(3,3),VALUE,X(3) !<The value to add
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,variable_idx,eqnset_idx
    INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX

    ENTERS("Burgers_PreSolveUpdateAnalyticValues",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STATIC_BURGERS_SUBTYPE,PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              !Loop over all the equation sets and set the appropriate field variable type BCs and
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
                          DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
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
                                          CALL Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET, &
                                            & ANALYTIC_FUNCTION_TYPE,X,TANGENTS,NORMAL,CURRENT_TIME,variable_type, &
                                            & GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS, &
                                            & VALUE,ERR,ERROR,*999)
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                            & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                            & FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                            BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                              & CONDITION_TYPES(local_ny)
                                            IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, & 
                                                & variable_type,FIELD_VALUES_SET_TYPE,local_ny, & 
                                                & VALUE,ERR,ERROR,*999)
                                            ENDIF
                                          ELSE
                                            CALL FlagError("Boundary conditions variable is not associated",ERR,ERROR,*999)
                                          ENDIF
                                        ENDDO !deriv_idx
                                      ENDDO !node_idx
                                    ELSE
                                      CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
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
                            CALL FlagError("Field variable is not associated.",ERR,ERROR,*999)
                          ENDIF
                          
                        ENDDO !variable_idx
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
                            & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Equations set geometric field is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
              ENDDO !eqnset_idx
            ELSE
              CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a BURGERS equation type of a fluid mechanics problem class."
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
    
    EXITS("Burgers_PreSolveUpdateAnalyticValues")
    RETURN
999 ERRORS("Burgers_PreSolveUpdateAnalyticValues",ERR,ERROR)
    EXITS("Burgers_PreSolveUpdateAnalyticValues")
    RETURN 1
    
  END SUBROUTINE Burgers_PreSolveUpdateAnalyticValues


  !   
  !================================================================================================================================
  !
  SUBROUTINE Burgers_PreSolveStoreCurrentSolution(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("Burgers_PreSolveStoreCurrentSolution",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
              ! do nothing ???
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a Burgers equation type of a fluid mechanics problem class."
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

    EXITS("Burgers_PreSolveStoreCurrentSolution")
    RETURN
999 ERRORS("Burgers_PreSolveStoreCurrentSolution",ERR,ERROR)
    EXITS("Burgers_PreSolveStoreCurrentSolution")
    RETURN 1
    
  END SUBROUTINE Burgers_PreSolveStoreCurrentSolution
  
  !   
  !================================================================================================================================
  !
    
  !>Sets up the Burgers problem post solve.
  SUBROUTINE BURGERS_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BURGERS_EQUATION_POST_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STATIC_BURGERS_SUBTYPE,PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
            !CALL BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a BURGERS type of a fluid mechanics problem class."
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

    EXITS("BURGERS_EQUATION_POST_SOLVE")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_POST_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BURGERS_EQUATION_POST_SOLVE
  
  !   
  !================================================================================================================================
  !
  
  !>Output data post solve
  SUBROUTINE BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP 
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    ENTERS("BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN        
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification array is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_BURGERS_SUBTYPE,PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
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
 
                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,ERR,ERROR,*999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a BURGERS equation type of a fluid mechanics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver solvers is not associated.",ERR,ERROR,*999)
    ENDIF
  ELSE
    CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
  ENDIF
  
    EXITS("BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Burgers problem.
  SUBROUTINE Burgers_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("Burgers_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STATIC_BURGERS_SUBTYPE, &
            & PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Burgers type of a fluid mechanics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_BURGERS_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Burgers equation problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("Burgers_ProblemSpecificationSet")
    RETURN
999 ERRORS("Burgers_ProblemSpecificationSet",err,error)
    EXITS("Burgers_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Burgers_ProblemSpecificationSet

  !
  !================================================================================================================================
  !
 
  !>Sets up the Burgers problem.
  SUBROUTINE BURGERS_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Burgers problem on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    
    ENTERS("BURGERS_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
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
              & " is invalid for a Burgers problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
                & " is invalid for a Burgers problem."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
            !Set the solver to be a static nonlinear solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Nonlinear solver",ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Burgers problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
              & " is invalid for a Burgers problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a Burgers problem."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
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
              & " is invalid for a Burgers problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
              & " is invalid for a Burgers problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
            !Set the solver to be a static nonlinear solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Nonlinear dynamic solver",ERR,ERROR,*999)
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
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a nonlinear burgers problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC, &
              & ERR,ERROR,*999)
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
              & " is invalid for a Burgers problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a Burgers problem."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a Burgers problem subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BURGERS_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BURGERS_EQUATION_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices for a BURGERS equation finite element equations set.
  SUBROUTINE Burgers_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,nj,ms,nh,nhs,ni,ns,MESH_COMPONENT1,MESH_COMPONENT2
    REAL(DP) :: C_PARAM,JGW,SUM1,SUM2,SUM3,DXI_DX,DPHINS_DXI,PHIMS,PHINS,U_VALUE,U_DERIV,VALUE
    LOGICAL :: EVALUATE_JACOBIAN,UPDATE_JACOBIAN_MATRIX
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME1,QUADRATURE_SCHEME2

    ENTERS("Burgers_FiniteElementJacobianEvaluate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
        UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
        EVALUATE_JACOBIAN=UPDATE_JACOBIAN_MATRIX
        IF(EVALUATE_JACOBIAN) THEN
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) &
            & CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          C_PARAM=1.0_DP
          !Loop over all Gauss points 
          DO ng=1,QUADRATURE_SCHEME1%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
            ENDIF
            !Loop over rows
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
              DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                !Loop over element columns
                DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                    &(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                  DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                    nhs=nhs+1
                    PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                    !Loop over xi directions
                    SUM1=0.0_DP
                    DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      U_DERIV=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR% &
                        & VALUES(mh,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni))
                      DXI_DX=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nh)
                      SUM1=SUM1+U_DERIV*DXI_DX
                    ENDDO !ni
                    SUM1=SUM1*PHINS
                    SUM2=0.0_DP
                    IF(nh==mh) THEN
                      !Loop over spatial directions
                      DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                        U_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(nj,NO_PART_DERIV)
                        SUM3=0.0_DP
                        DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                          DPHINS_DXI=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          DXI_DX=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                          SUM3=SUM3+DPHINS_DXI*DXI_DX
                        ENDDO !ni
                        SUM2=SUM2+U_VALUE*SUM3
                      ENDDO !nj
                    ENDIF
                    VALUE=C_PARAM*(SUM1+SUM2)*PHIMS
                    JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+VALUE*JGW
                  ENDDO !ns    
                ENDDO !nh
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("Burgers_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("Burgers_FiniteElementJacobianEvaluate",ERR,ERROR)
    EXITS("Burgers_FiniteElementJacobianEvaluate")
    RETURN 1
    
  END SUBROUTINE Burgers_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Burgers equation finite element equations set.
  SUBROUTINE Burgers_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nj,nh,nhs,ni,mi,ns
    INTEGER(INTG) MESH_COMPONENT1,MESH_COMPONENT2
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,JGW,SUM,SUM1,PHIMS,PHINS,DPHIMS_DXI(3),DPHINS_DXI(3),U_VALUE
    LOGICAL :: EVALUATE_ANY,EVALUATE_DAMPING,EVALUATE_LINEAR_DYNAMIC,EVALUATE_RESIDUAL,EVALUATE_RHS, &
      & EVALUATE_STIFFNESS,FIRST_DAMPING,FIRST_RHS,FIRST_STIFFNESS,UPDATE_STIFFNESS,UPDATE_DAMPING,UPDATE_RHS,UPDATE_RESIDUAL
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: STIFFNESS_MATRIX,DAMPING_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: GEOMETRIC_VARIABLE,FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    ENTERS("Burgers_FiniteElementResidualEvaluate",ERR,ERROR,*999)

    NULLIFY(DAMPING_MATRIX)
    NULLIFY(DYNAMIC_MATRICES)
    NULLIFY(STIFFNESS_MATRIX)

    UPDATE_STIFFNESS=.FALSE.
    FIRST_STIFFNESS=.FALSE.
    UPDATE_DAMPING=.FALSE.
    FIRST_DAMPING=.FALSE.
    UPDATE_RHS=.FALSE.
    FIRST_RHS=.FALSE.
    UPDATE_RESIDUAL=.FALSE.

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
        DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
        GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
        GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
        MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
        GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          STIFFNESS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a BURGERS equation type of a fluid mechanics equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        IF(ASSOCIATED(STIFFNESS_MATRIX)) THEN
          UPDATE_STIFFNESS=STIFFNESS_MATRIX%UPDATE_MATRIX
          FIRST_STIFFNESS=STIFFNESS_MATRIX%FIRST_ASSEMBLY
        ENDIF
        IF(ASSOCIATED(DAMPING_MATRIX)) THEN
          UPDATE_DAMPING=DAMPING_MATRIX%UPDATE_MATRIX
          FIRST_DAMPING=DAMPING_MATRIX%FIRST_ASSEMBLY
        ENDIF
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          UPDATE_RHS=RHS_VECTOR%UPDATE_VECTOR
          FIRST_RHS=RHS_VECTOR%FIRST_ASSEMBLY
        ENDIF
        IF(ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
        EVALUATE_RESIDUAL=UPDATE_RESIDUAL
        EVALUATE_STIFFNESS=FIRST_STIFFNESS.OR.UPDATE_STIFFNESS
        EVALUATE_DAMPING=FIRST_DAMPING.OR.UPDATE_DAMPING
        EVALUATE_RHS=FIRST_RHS.OR.UPDATE_RHS
        EVALUATE_LINEAR_DYNAMIC=EVALUATE_STIFFNESS.OR.EVALUATE_DAMPING.OR.EVALUATE_RHS
        EVALUATE_ANY=EVALUATE_LINEAR_DYNAMIC.OR.UPDATE_RESIDUAL
        IF(EVALUATE_ANY) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) &
            & CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          A_PARAM=1.0_DP
          B_PARAM=1.0_DP
          C_PARAM=1.0_DP
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & MATERIALS_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
                A_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                B_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
                C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
              ELSE
                B_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              ENDIF
            ENDIF
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                IF(EVALUATE_LINEAR_DYNAMIC) THEN
                  IF(EVALUATE_STIFFNESS.OR.EVALUATE_DAMPING) THEN
                    nhs=0
                    !Loop over element columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                        &(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        !Diffusion matrix
                        IF(EVALUATE_STIFFNESS) THEN
                          DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            DPHIMS_DXI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                            DPHINS_DXI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          END DO !ni
                          SUM=0.0_DP
                          !Calculate SUM 
                          DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                              SUM=SUM+DPHINS_DXI(ni)*DPHIMS_DXI(mi)*EQUATIONS%INTERPOLATION% &
                                & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%GU(mi,ni)
                            ENDDO !ni
                          ENDDO !mi
                          STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)- &
                            & B_PARAM*SUM*JGW
                        ENDIF
                        !Mass matrix
                        IF(EVALUATE_DAMPING) THEN
                          PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                          DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & A_PARAM*PHIMS*PHINS*JGW
                        ENDIF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF !Stiffness or Damping
                  !Calculate RHS
                  IF(EVALUATE_RHS) THEN
                    RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
                  ENDIF
                ENDIF !Evaluate linear dynamic
                !Calculate nonlinear vector
                IF(EVALUATE_RESIDUAL) THEN
                  PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  SUM=0.0_DP
                  DO nj=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    U_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(nj,NO_PART_DERIV)
                    SUM1=0.0_DP
                    !Calculate SUM 
                    DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      SUM1=SUM1+EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(nj, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni))*EQUATIONS%INTERPOLATION% &
                        & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
                    ENDDO !ni
                    SUM=SUM+U_VALUE*SUM1
                  ENDDO !nj
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES% &
                    & ELEMENT_RESIDUAL%VECTOR(mhs)+C_PARAM*SUM*PHIMS*JGW
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(EVALUATE_STIFFNESS.OR.EVALUATE_DAMPING) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(EVALUATE_STIFFNESS) &
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      IF(EVALUATE_DAMPING) &
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(EVALUATE_RHS) &
                  & RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                IF(EVALUATE_RESIDUAL) &
                  & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF
          
        ENDIF !Evaluate any
        
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("Burgers_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("Burgers_FiniteElementResidualEvaluate",ERR,ERROR)
    EXITS("Burgers_FiniteElementResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE Burgers_FiniteElementResidualEvaluate
 
  !
  !================================================================================================================================
  !

END MODULE BURGERS_EQUATION_ROUTINES
