!> \file
!> \author Chris Bradley
!> \brief This module handles all Hamilton-Jacobi equations routines.
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

!>This module handles all Hamilton-Jacobi equations routines.
MODULE HAMILTON_JACOBI_EQUATIONS_ROUTINES

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
  USE MATHS
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

!!MERGE: move

  PUBLIC HJ_BoundaryConditionsAnalyticCalculate
  
  PUBLIC HJEquation_EquationsSetSolutionMethodSet
  
  PUBLIC HJ_EQUATION_EQUATIONS_SET_SETUP
  
  PUBLIC HJEquation_EquationsSetSpecificationSet

  PUBLIC HJ_EQUATION_FINITE_ELEMENT_CALCULATE
  
  PUBLIC HJ_EQUATION_PROBLEM_SETUP
  
  PUBLIC HJEquation_ProblemSpecificationSet
  
  PUBLIC NUMBER_OF_INPUT_NODES,PRE_PROCESS_INFORMATION,SOLVE_PROBLEM_FMM,SOLVE_PROBLEM_GEODESIC
  PUBLIC SOLVE_PROBLEM_GEODESIC_CONNECTIVITY,SOLVE_PROBLEM_FMM_CONNECTIVITY
  PUBLIC FIND_MINIMAX,POST_PROCESS_DATA

CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE HJ_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type
    REAL(DP) :: VALUE,X(3)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR    
    
    ENTERS("HJ_BoundaryConditionsAnalyticCalculate",ERR,ERROR,*999)

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
                                CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_1)
                                  !u=x^2+2.x.y-y^2
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=X(1)*X(1)-2.0_DP*X(1)*X(2)-X(2)*X(2)
                                    CASE(GLOBAL_DERIV_S1)
                                      VALUE=2.0_DP*X(1)+2.0_DP*X(2)
                                    CASE(GLOBAL_DERIV_S2)
                                      VALUE=2.0_DP*X(1)-2.0_DP*X(2)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      VALUE=2.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_2)
                                  !u=cos(x).cosh(y)
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=COS(X(1))*COSH(X(2))
                                    CASE(GLOBAL_DERIV_S1)
                                      VALUE=-SIN(X(1))*COSH(X(2))
                                    CASE(GLOBAL_DERIV_S2)
                                      VALUE=COS(X(1))*SINH(X(2))
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      VALUE=-SIN(X(1))*SINH(X(2))
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
                                    CASE(GLOBAL_DERIV_S1)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_1)
                                  !u=x^2+y^2-2.z^2
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=X(1)*X(1)+X(2)*X(2)-2.0_DP*X(3)*X(3)
                                    CASE(GLOBAL_DERIV_S1)
                                      VALUE=2.0_DP*X(1)
                                    CASE(GLOBAL_DERIV_S2)
                                      VALUE=2.0_DP*X(2)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S3)
                                      VALUE=-4.0_DP*X(3)
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      VALUE=0.0_DP
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S3)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_2)
                                  !u=cos(x).cosh(y).z
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=COS(X(1))*COSH(X(2))*X(3)
                                    CASE(GLOBAL_DERIV_S1)
                                      VALUE=-SIN(X(1))*COSH(X(2))*X(3)
                                    CASE(GLOBAL_DERIV_S2)
                                      VALUE=COS(X(1))*SINH(X(2))*X(3)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      VALUE=-SIN(X(1))*SINH(X(2))*X(3)
                                    CASE(GLOBAL_DERIV_S3)
                                      VALUE=COS(X(1))*COSH(X(2))
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      VALUE=-SIN(X(1))*COSH(X(2))
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      VALUE=COS(X(1))*SINH(X(2))
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      VALUE=-SIN(X(1))*SINH(X(2))
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
                                    CASE(GLOBAL_DERIV_S1)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S3)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      !CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                          & ERR,ERROR))//" is invalid."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                                      & " is invalid."
                                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                CASE DEFAULT
                                  LOCAL_ERROR="The analytic function type of "// &
                                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                    & " is invalid."
                                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",ERR,ERROR,*999)
                ENDIF

              ENDDO !variable_idx
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
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
    
    EXITS("HJ_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("HJ_BoundaryConditionsAnalyticCalculate",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJ_BoundaryConditionsAnalyticCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Hamilton-Jacobi equation finite element equations set.
  SUBROUTINE HJ_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,mi,ms,nh,nhs,ni,ns
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CHARACTER(26) :: CVAR
    INTEGER :: GAUSS_POINT_LOOP_PHASE(2) = (/ 0, 0 /)
    SAVE GAUSS_POINT_LOOP_PHASE
#endif

    ENTERS("HJ_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Hamilton-Jacobi type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Hamilton-Jacobi is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
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
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
#ifdef TAUPROF
              WRITE (CVAR,'(a17,i2)') 'Gauss Point Loop ',ng
              CALL TAU_PHASE_CREATE_DYNAMIC(GAUSS_POINT_LOOP_PHASE,CVAR)
              CALL TAU_PHASE_START(GAUSS_POINT_LOOP_PHASE)
#endif
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
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
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
#ifdef TAUPROF
            CALL TAU_PHASE_STOP(GAUSS_POINT_LOOP_PHASE)
#endif
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
              ENDDO !ms
            ENDDO !mh
          ENDIF       
        CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a Hamilton-Jacobi equation type of a classical field equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HJ_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("HJ_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJ_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Calculates the connectivity and seed values for a Hamilton-Jacobi equation fast marching equations set.
  SUBROUTINE HJ_EQUATION_FAST_MARCHING_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,mi,ms,nh,nhs,ni,ns
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CHARACTER(26) :: CVAR
    INTEGER :: GAUSS_POINT_LOOP_PHASE(2) = (/ 0, 0 /)
    SAVE GAUSS_POINT_LOOP_PHASE
#endif

    ENTERS("HJ_EQUATION_FAST_MARCHING_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Hamilton-Jacobi type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Hamilton-Jacobi is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
!          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD!
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
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
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
#ifdef TAUPROF
              WRITE (CVAR,'(a17,i2)') 'Gauss Point Loop ',ng
              CALL TAU_PHASE_CREATE_DYNAMIC(GAUSS_POINT_LOOP_PHASE,CVAR)
              CALL TAU_PHASE_START(GAUSS_POINT_LOOP_PHASE)
#endif
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
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
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
#ifdef TAUPROF
            CALL TAU_PHASE_STOP(GAUSS_POINT_LOOP_PHASE)
#endif
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
              ENDDO !ms
            ENDDO !mh
          ENDIF       
        CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a Hamilton-Jacobi equation type of a classical field equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HJ_EQUATION_FAST_MARCHING_CALCULATE")
    RETURN
999 ERRORSEXITS("HJ_EQUATION_FAST_MARCHING_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJ_EQUATION_FAST_MARCHING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the Hamilton-Jacobi equation type of a classical field equations set class.
  SUBROUTINE HJ_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Hamilton-Jacobi equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HJ_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Hamilton-Jacobi type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      
      CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
        CALL HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Hamilton-Jacobi equation type of a classical field equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HJ_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("HJ_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJ_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Hamilton-Jacobi equation type of an classical field equations set class.
  SUBROUTINE HJEquation_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HJEquation_EquationsSetSolutionMethodSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Hamilton-Jacobi type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)        
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
!        CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
!          CALL FlagError("Not implemented.",ERR,ERROR,*999)
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
      CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)        
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
!        CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
!          CALL FlagError("Not implemented.",ERR,ERROR,*999)
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
          & " is not valid for a Hamilton-Jacobi equation type of an classical field equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HJEquation_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("HJEquation_EquationsSetSolutionMethodSet",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJEquation_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Hamilton-Jacobi equation type of a classical field equations set class.
  SUBROUTINE HJEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    CALL Enters("HJEquation_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Hamilton-Jacobi type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
        !ok
      CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a Hamilton-Jacobi type of a classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_HJ_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("HJEquation_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("HJEquation_EquationsSetSpecificationSet",err,error)
    EXITS("HJEquation_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE HJEquation_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the standard Hamilton-Jacobi equation.
  SUBROUTINE HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS, GEOMETRIC_MESH_COMPONENT
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS

    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HJ_EQUATION_EQUATION_SET_STANDARD_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Hamilton-Jacobi type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_STANDARD_HJ_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL HJEquation_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Hamilton-Jacobi equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
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
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"Phi",ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
                & ERR,ERROR,*999)
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
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"Phi",ERR,ERROR,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & "del Phi/del n",ERR,ERROR,*999)
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
!              CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
!                CALL FlagError("Not implemented.",ERR,ERROR,*999)
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
!              CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
!                CALL FlagError("Not implemented.",ERR,ERROR,*999)
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
              & " is invalid for a standard Hamilton-Jacobi equation"
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
            
            
            

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
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
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
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+2, &
                    & ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
            
            
            
            
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            

            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE) THEN
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
              CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF            
            
            
            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Hamilton-Jacobi equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Hamilton-Jacobi equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
                  CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_1)
                    !Check that we are in 2D
                    IF(NUMBER_OF_DIMENSIONS/=2) THEN
                      LOCAL_ERROR="The number of geometric dimensions of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that there be 2 geometric dimensions."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    !Create analytic field if required
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_HJ_EQUATION_TWO_DIM_1
                  CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_2)
                    !Check that we are in 2D
                    IF(NUMBER_OF_DIMENSIONS/=2) THEN
                      LOCAL_ERROR="The number of geometric dimensions of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that there be 2 geometric dimensions."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    !Create analytic field if required
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_HJ_EQUATION_TWO_DIM_2
                  CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_1)
                    !Check that we are in 3D
                    IF(NUMBER_OF_DIMENSIONS/=3) THEN
                      LOCAL_ERROR="The number of geometric dimensions of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that there be 3 geometric dimensions."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    !Create analytic field if required
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_HJ_EQUATION_THREE_DIM_1
                  CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_2)
                    !Check that we are in 3D
                    IF(NUMBER_OF_DIMENSIONS/=3) THEN
                      LOCAL_ERROR="The number of geometric dimensions of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that there be 3 geometric dimensions."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    !Create analytic field if required
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_HJ_EQUATION_THREE_DIM_2
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a standard Hamilton-Jacobi equation."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FlagError("Equations set geometric field is not associated.",ERR,ERROR,*999)
                ENDIF
             ELSE
                CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",ERR,ERROR,*999)
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
              CALL FlagError("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Hamilton-Jacobi equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EquationsMapping_LinearMatricesNumberSet(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EquationsMapping_LinearMatricesVariableTypesSet(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
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
                CALL EquationsMatrices_LinearStructureTypeSet(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
!            CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
!              CALL FlagError("Not implemented.",ERR,ERROR,*999)
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
              & " is invalid for a standard Hamilton-Jacobi equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Hamilton-Jacobi equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a standard Hamilton-Jacobi equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN
999 ERRORSEXITS("HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the Hamilton-Jacobi problem.
  SUBROUTINE HJ_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Hamilton-Jacobi equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HJ_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Hamilton-Jacobi problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_STANDARD_HJ_SUBTYPE)
        CALL HJ_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_HJ_SUBTYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Hamilton-Jacobi equation type of a classical field problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HJ_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("HJ_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJ_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Hamilton-Jacobi equation type.
  SUBROUTINE HJEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("HJEquation_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STANDARD_HJ_SUBTYPE)
          !ok
        CASE(PROBLEM_GENERALISED_HJ_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Hamilton-Jacobi type of a classical field problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_HJ_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Hamilton-Jacobi problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("HJEquation_ProblemSpecificationSet")
    RETURN
999 ERRORS("HJEquation_ProblemSpecificationSet",err,error)
    EXITS("HJEquation_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE HJEquation_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the standard Hamilton-Jacobi equations problem.
  SUBROUTINE HJ_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
    ENTERS("HJ_EQUATION_PROBLEM_STANDARD_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Hamilton-Jacobi problem.",err,error,*999)
      END IF
      IF(PROBLEM%SPECIFICATION(3)==PROBLEM_STANDARD_HJ_SUBTYPE) THEN
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
              & " is invalid for a standard Hamilton-Jacobi equation."
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
              & " is invalid for a standard Hamilton-Jacobi equation."
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
                & " is invalid for a standard Hamilton-Jacobi equation."
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
              & " is invalid for a standard Hamilton-Jacobi equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Hamilton-Jacobi equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a standard Hamilton-Jacobi equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HJ_EQUATION_PROBLEM_STANDARD_SETUP")
    RETURN
999 ERRORSEXITS("HJ_EQUATION_PROBLEM_STANDARD_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HJ_EQUATION_PROBLEM_STANDARD_SETUP

  !
  !================================================================================================================================
  !
!
!
!


  !>Calculates to give back the number of nodes from input file.
  SUBROUTINE NUMBER_OF_INPUT_NODES(INPUT_FILE_NAME,INPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,&
  &TOTAL_NUMBER_OF_CONNECTIVITY,Err)

    !subroutine variables
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_NODES
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_CONNECTIVITY
    CHARACTER (LEN=300) :: INPUT_FILE_NAME
    CHARACTER (LEN=10) :: INPUT_FILE_FORMAT
    INTEGER(INTG) :: Err

    !Argument variables
!    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string

    !Local variables
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:):: CONNECTIVITY_LIST
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:):: ELEMENT_LIST
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: CONNECTIVITY_NUMBER
    INTEGER(INTG) :: I,J,K,N,NUMBER_OF_NODES_PER_ELEMENT,THERE_IS_IN_CONNECTIVITY_LIST
    CHARACTER (LEN=300) :: STRING
        
!    ENTERS("GENERATE_STATUS_MASK",Err,Error,*999)
    
! """""""""""""""""""""""""""""""""""INPUT OF TABC FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TABC") THEN 

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".tabc"
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
!      PRINT *, STRING
      TOTAL_NUMBER_OF_NODES=-1
      DO WHILE (STRING .ne. "Connectivity") 
        READ(11,*) STRING
        TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+1
      ENDDO
      
      TOTAL_NUMBER_OF_CONNECTIVITY=0
      DO I=1,TOTAL_NUMBER_OF_NODES
        READ(11,*) STRING,J
        TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+J
      ENDDO
      
      CLOSE (11)
      TOTAL_NUMBER_OF_ELEMENTS = TOTAL_NUMBER_OF_CONNECTIVITY ! SHOULD BE DEFINED
    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".vtk"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_NODES
!      PRINT *, STRING,TOTAL_NUMBER_OF_NODES
      
      DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)
      
        READ(11,*) STRING

      ENDDO
!      READ(11,*) STRING
!      print*,I,INT(TOTAL_NUMBER_OF_NODES/3.0+0.5),STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_ELEMENTS
    
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET1NPL") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".vtk"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_NODES,STRING
!      PRINT *, STRING,TOTAL_NUMBER_OF_NODES,STRING
      
      DO I=1,TOTAL_NUMBER_OF_NODES
      
        READ(11,*) STRING

      ENDDO
!      READ(11,*) STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_ELEMENTS

      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF CARP FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "CARP") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".pts"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_NODES
!      PRINT *, TOTAL_NUMBER_OF_NODES
      CLOSE (11)
      
      STRING = INPUT_FILE_NAME(1:I)//".elem"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF TETGEN FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TETGEN") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".node"
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_NODES
      CLOSE (11)
      
      STRING = INPUT_FILE_NAME(1:I)//".ele"
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE (11)
      
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

    ENDIF


  END SUBROUTINE NUMBER_OF_INPUT_NODES

  !
  !================================================================================================================================
  !


  !>to READ input file.
  SUBROUTINE PRE_PROCESS_INFORMATION(MATERIAL_BEHAVIOUR,INPUT_FILE_NAME,INPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,&
&INPUT_TYPE_FOR_SEED_VALUE,INPUT_TYPE_FOR_SPEED_FUNCTION,SPEED_FUNCTION_ALONG_EIGEN_VECTOR,INPUT_TYPE_FOR_CONDUCTIVITY,&
&STATUS_MASK,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,&
&SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

    !subroutine variables
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LIST(:,:)

    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: COLUMN_INDEX
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: RAW_INDEX
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    REAL(DP) :: SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS,TOTAL_NUMBER_OF_CONNECTIVITY
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NODES_PER_ELEMENT
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SEED_VALUE
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SPEED_FUNCTION
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_CONDUCTIVITY
    CHARACTER (LEN=300) :: INPUT_FILE_NAME
    CHARACTER (LEN=10)  :: INPUT_FILE_FORMAT
    CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    
    !Local variables
    CHARACTER (LEN=300) :: STRING
    INTEGER(INTG) :: I,J,K,N,FIRST_NODE_NUMBER
    INTEGER(INTG) :: TEXT_LENGTH, THERE_IS_IN_CONNECTIVITY_LIST
    REAL(DP) :: A(3),B(3),C(3)
    REAL(DP) :: DOT_PRODUCT_VALUE
        
!INITIALIZE PARAMETERS:
    DO I=1,TOTAL_NUMBER_OF_NODES
      CONNECTIVITY_NUMBER(I)=0
      DO J=1,3
        NODE_LIST(I,J) = 0.0
        SPEED_FUNCTION_TABLE(I,J) = 0.0
      ENDDO
      DO J=1,9
        CONDUCTIVITY_TENSOR(I,J) = 0.0
      ENDDO
      RAW_INDEX(I)=0
    ENDDO
    RAW_INDEX(TOTAL_NUMBER_OF_NODES+1)=0

    DO I=1,TOTAL_NUMBER_OF_CONNECTIVITY
      COLUMN_INDEX(I)=0
      DO J=1,3
        SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(I,J) = 0.0
      ENDDO
      DO J=1,9
        CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(I,J) = 0.0
      ENDDO
    ENDDO
     
! """""""""""""""""""""""""""""""""""INPUT OF TABC FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TABC") THEN 

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".tabc"
      OPEN (11,FILE=STRING)
      READ(11,*) STRING

! SOSIOISOISOISOISOIS      load data for the case material behaves = ISOTROPIC 
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!      input type for *velocity function* = FILE and *seed points* = FILE
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),SEED_VALUE(I)
            
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
           
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

          ENDDO
        ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SEED_VALUE(I)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

          ENDDO
        ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

            IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
              SEED_VALUE(I) = 1000.0_DP
            ENDIF

          ENDDO
        ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

            IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
              SEED_VALUE(I) = 1000.0_DP
            ENDIF

          ENDDO
        ENDIF

      ENDIF

! ANANANAIANSOANAIANI      load data for the case material behaves = ANISOTROPIC 
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 


!      if conductivity format is TENSOR type i.e. three EigenVectors
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN

!      input type for *velocity function* = FILE and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3),SEED_VALUE(I)

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SEED_VALUE(I)
 
            ENDDO
          ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

        ENDIF


!      if conductivity format is VECTOR type i.e. first EigenVectors
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN

!      input type for *velocity function* = FILE and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3),SEED_VALUE(I)

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SEED_VALUE(I)
 
            ENDDO
          ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

          DO I=1,TOTAL_NUMBER_OF_NODES
!            CALL CALCULATE_SECOND_EIGENVECTOR()
            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            CALL CALCULATE_SECOND_EIGENVECTOR()
          ENDDO

        ENDIF

      ENDIF

! CONNSDONCOCNCNOSKCN      load data for the CONNECTIVITY list 
      READ(11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_NODES

        READ(11,*) STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
        
        DO J=1,3
          SPEED_FUNCTION_TABLE(I,J)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
        ENDDO
          
!        PRINT *,STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))

      ENDDO

      CLOSE(11)

    ENDIF


! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"

      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING

      DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-1

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3),(NODE_LIST(3*(I-1)+3,J),J=1,3)

      ENDDO

      I=INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)

      IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 0) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3),(NODE_LIST(3*(I-1)+3,J),J=1,3)

      ENDIF
      IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 1) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3)


      ENDIF
      IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 2) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3)

      ENDIF
      

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 

      READ(11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO
        
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO
        
      ENDDO
      
      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        DO I=1,TOTAL_NUMBER_OF_NODES
          STATUS_MASK(I) = ""
        ENDDO

        OPEN (11,FILE=STRING)
        READ (11,*) N

        DO I=1,N 

          READ(11,*) J,SEED_VALUE(J+1)
!          PRINT*,(NODE_LIST(J+1,K),K=1,3),SEED_VALUE(J+1)      
          STATUS_MASK(J+1) = "SEED POINT"

        ENDDO

        CLOSE(11)

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          IF (STRING .EQ. '#') THEN! begin if

            DO WHILE (STRING .NE. 'fiber')
      
              READ(11,*) STRING
        
            ENDDO
      
!      PRINT *,STRING
      
            DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-1

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDDO

            I=INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)

            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 0) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 1) THEN
    
              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 2) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3)

            ENDIF
      
          
          ELSE
            DO I=1,TOTAL_NUMBER_OF_NODES 
            
              READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3) 
              
            ENDDO               
          ENDIF ! end if
          
          DO I=1,TOTAL_NUMBER_OF_NODES 
          
!            READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)

            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 
        
        DO I=1,TOTAL_NUMBER_OF_NODES
          DO J=1,3
            SPEED_FUNCTION_TABLE(I,J) = 0.0_DP
          ENDDO
        ENDDO
        
          DO I=1,TOTAL_NUMBER_OF_NODES
!            DO J=1,CONNECTIVITY_NUMBER(I)
!              IF (CONNECTIVITY_LIST(I,J) > TOTAL_NUMBER_OF_NODES) THEN
!                PRINT*, I,J,CONNECTIVITY_NUMBER(I),CONNECTIVITY_LIST(I,J),TOTAL_NUMBER_OF_NODES
!              ENDIF
              SPEED_FUNCTION_TABLE(I,1) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
              SPEED_FUNCTION_TABLE(I,2) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)
              SPEED_FUNCTION_TABLE(I,3) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),1) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),2) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),3) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
!            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET1NPL") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"

      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_NODES

        READ(11,*) (NODE_LIST(I,J),J=1,3)

      ENDDO

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 

      READ(11,*) STRING

      DO N=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) N

        DO I=1,N 

          READ(11,*) J,SEED_VALUE(J+1)
       
          STATUS_MASK(J+1) = "SEED POINT"

        ENDDO

        CLOSE(11)

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING
!          print *,STRING
        
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          IF (STRING .EQ. '#') THEN! begin if

            DO WHILE (STRING .NE. 'fiber')
      
              READ(11,*) STRING
        
            ENDDO
      
!      PRINT *,STRING
      
            DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-1

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDDO

            I=INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)

            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 0) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 1) THEN
    
              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 2) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3)

            ENDIF
      
          
          ELSE
            DO I=1,TOTAL_NUMBER_OF_NODES 
            
              READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3) 
              
            ENDDO               
          ENDIF ! end if
          
          DO I=1,TOTAL_NUMBER_OF_NODES 
          
!            READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF CARP FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "CARP") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".pts"

      OPEN (11,FILE=STRING)
      READ (11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_NODES 

        READ(11,*) (NODE_LIST(I,J),J=1,3)

      ENDDO

      CLOSE(11)

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".elem"

      OPEN (11,FILE=STRING)
      READ (11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      DO N=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        DO I=1,TOTAL_NUMBER_OF_NODES 
 
          READ(11,*) STRING,SEED_VALUE(I)

        ENDDO

        CLOSE(11)

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".lon"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR" .OR. INPUT_TYPE_FOR_CONDUCTIVITY .EQ. " ") THEN
          DO I=1,TOTAL_NUMBER_OF_ELEMENTS
            READ(11,*) (CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),1),CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),2), &
             & CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),4),CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),5), &
             & CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF
            
            DO J=2,NUMBER_OF_NODES_PER_ELEMENT
              DO K=1,9
                CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,J),K)=CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),K)
              ENDDO
            ENDDO

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF
    
! """""""""""""""""""""""""""""""""""INPUT OF TETGEN FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TETGEN") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".node"

      OPEN (11,FILE=STRING)
      READ (11,*) STRING

      READ(11,*) FIRST_NODE_NUMBER,(NODE_LIST(1,J),J=1,3)

      DO I=2,TOTAL_NUMBER_OF_NODES 


        READ(11,*) STRING,(NODE_LIST(I,J),J=1,3)


      ENDDO

      CLOSE(11)

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".ele"

      OPEN (11,FILE=STRING)
      READ (11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      DO N=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        IF (FIRST_NODE_NUMBER .EQ. 0) THEN        
          DO J=1,NUMBER_OF_NODES_PER_ELEMENT
            ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
          ENDDO
        ENDIF

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        DO I=1,TOTAL_NUMBER_OF_NODES 
 
          READ(11,*) STRING,SEED_VALUE(I)

        ENDDO

        CLOSE(11)

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE

                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF

    ! to set RAW_INDEX and COLUMN_INDEX
    RAW_INDEX(1) = 0
    DO I=1,TOTAL_NUMBER_OF_NODES
        
      RAW_INDEX(I+1) = RAW_INDEX(I) + CONNECTIVITY_NUMBER(I)
      DO J = 1,CONNECTIVITY_NUMBER(I)
        COLUMN_INDEX(RAW_INDEX(I)+J) = CONNECTIVITY_LIST(I,J)
      ENDDO          

    ENDDO
    
    DO I=1,TOTAL_NUMBER_OF_NODES
        
      DO J = RAW_INDEX(I)+1,RAW_INDEX(I+1)
        DO K = 1,3
          SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,K) = &
          & (SPEED_FUNCTION_TABLE(I,K)+SPEED_FUNCTION_TABLE(COLUMN_INDEX(J),K))/2.0_DP
        ENDDO
        DO K = 1,9
          CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,K) = &
          & (CONDUCTIVITY_TENSOR(I,K)+CONDUCTIVITY_TENSOR(COLUMN_INDEX(J),K))/2.0_DP
        ENDDO
      ENDDO

    ENDDO

!    EXITS("GENERATE_STATUS_MASK")
!    RETURN
999 ERRORSEXITS("GENERATE_STATUS_MASK",ERR,ERROR)
    RETURN

  END SUBROUTINE PRE_PROCESS_INFORMATION


  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_FMM(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,&
  &SEED_VALUE,CONNECTIVITY_NUMBER,CONNECTIVITY_LIST,STATUS_MASK)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER, CHANGED_STATUS, MIN_TRIAL_NODE, TRIAL_STATUS
    REAL(DP), DIMENSION(3)   :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW,CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    !Start Program

    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

!    MAIN LOOP 
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=1,CONNECTIVITY_NUMBER(MIN_TRIAL_NODE)
        TIME_NEW=1000

        DO J=1,CONNECTIVITY_NUMBER(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))
          IF (STATUS_MASK(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=(&
                             &/NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),1)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),2)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),3)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3)/)

!            CONDUCTION_RATIO=SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)/SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)
!            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE = (/1.0_DP,0.0_DP,0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO/), SHAPE = (/3,3/))
            CONDUCTION_RATIO= &
            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)/&
            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)
!            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO/), SHAPE = (/3,3/))


!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),4),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),5),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),6),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),7),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),8),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),9)/&
                               &),SHAPE = (/3,3/))

            CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
            CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)
!            CALL INVERT(FMFT,INV_FMFT,DET,Err,Error,*999)

!	    PRINT *,F(1,1),F(1,2),F(1,3),F(2,1),F(2,2),F(2,3),F(3,1),F(3,2),F(3,3)

            CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J))+SQRT(ABS(VMV))*&
          &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)
!          &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)


            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .lt. SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))) THEN
          SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = TIME_NEW
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "KNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0.AND.STATUS_MASK(I) .EQ. "SEED POINT".AND.SEED_VALUE(I).LT.MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


999 ERRORSEXITS("SOLVE_PROBLEM_FMM",ERR,ERROR)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_FMM

  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_FMM_CONNECTIVITY(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_CONNECTIVITY
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDEX(:)  
    INTEGER(INTG), ALLOCATABLE :: RAW_INDEX(:)    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)

    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER, CHANGED_STATUS, MIN_TRIAL_NODE, TRIAL_STATUS
    REAL(DP), DIMENSION(3)   :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(2)   :: CONDUCTION_RATIO
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    !Start Program

    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

! ::::::: MAIN LOOP :::::::	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=RAW_INDEX(MIN_TRIAL_NODE)+1,RAW_INDEX(MIN_TRIAL_NODE+1)
        TIME_NEW=1000

        DO J=RAW_INDEX(COLUMN_INDEX(I))+1,RAW_INDEX(COLUMN_INDEX(I)+1)
          IF (STATUS_MASK(COLUMN_INDEX(J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=(/NODE_LIST(COLUMN_INDEX(I),1)-NODE_LIST(COLUMN_INDEX(J),1)&
                            &,NODE_LIST(COLUMN_INDEX(I),2)-NODE_LIST(COLUMN_INDEX(J),2)&
                            &,NODE_LIST(COLUMN_INDEX(I),3)-NODE_LIST(COLUMN_INDEX(J),3)/)

            CONDUCTION_RATIO(1)= SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,2)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
            CONDUCTION_RATIO(2)= SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,3)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO(1)*CONDUCTION_RATIO(1),0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO(2)*CONDUCTION_RATIO(2)/), SHAPE = (/3,3/))


!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,1),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,2),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,3),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,4),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,5),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,6),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,7),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,8),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,9)/),SHAPE = (/3,3/))

            CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
            CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)
!            CALL INVERT(FMFT,INV_FMFT,DET,Err,Error,*999)

            CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(COLUMN_INDEX(J))+SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .LT. SEED_VALUE(COLUMN_INDEX(I))) THEN
          SEED_VALUE(COLUMN_INDEX(I)) = TIME_NEW
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "KNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0 .AND. STATUS_MASK(I) .EQ. "SEED POINT" .AND. SEED_VALUE(I) .LT. MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


999 ERRORSEXITS("SOLVE_PROBLEM_FMM_CONNECTIVITY",ERR,ERROR)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_FMM_CONNECTIVITY


  !
  !================================================================================================================================
  !

  SUBROUTINE SOLVE_PROBLEM_GEODESIC_CONNECTIVITY(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK,TRACE_NODE)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_CONNECTIVITY
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDEX(:)
    INTEGER(INTG), ALLOCATABLE :: RAW_INDEX(:)
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: TRACE_NODE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    REAL(DP), ALLOCATABLE :: CONNECTIVITY_WEIGHT(:)
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: CHANGED_STATUS,MIN_TRIAL_NODE,TRIAL_STATUS,TRACE_NODE_NUMBER
    REAL(DP), DIMENSION(3) :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW
    REAL(DP), DIMENSION(2) :: CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    
    !initialize data
    DO I=1,TOTAL_NUMBER_OF_NODES
      TRACE_NODE(I)=0
    ENDDO

    !Start Program
    ALLOCATE(CONNECTIVITY_WEIGHT(TOTAL_NUMBER_OF_CONNECTIVITY),STAT=ERR)
    DO I=1,TOTAL_NUMBER_OF_NODES
      DO J=RAW_INDEX(I)+1,RAW_INDEX(I+1)
        DISTANCE_VECTOR=(/NODE_LIST(I,1)-NODE_LIST(COLUMN_INDEX(J),1)&
                        &,NODE_LIST(I,2)-NODE_LIST(COLUMN_INDEX(J),2)&
                        &,NODE_LIST(I,3)-NODE_LIST(COLUMN_INDEX(J),3)/)

        CONDUCTION_RATIO(1) = SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,2)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
        CONDUCTION_RATIO(2) = SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,3)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

        CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                            &0.0_DP,CONDUCTION_RATIO(1)*CONDUCTION_RATIO(1),0.0_DP,&
                                            &0.0_DP,0.0_DP,CONDUCTION_RATIO(2)*CONDUCTION_RATIO(2)/), SHAPE = (/3,3/))

        F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,1),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,2),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,3),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,4),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,5),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,6),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,7),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,8),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,9)/),SHAPE = (/3,3/))

        CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
        CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
        CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)

        CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
        CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

        CONNECTIVITY_WEIGHT(J)=SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
      ENDDO
    ENDDO
          
    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

! ::::::: MAIN LOOP :::::::	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
!    LOOP_NUMBER = 0
    PRINT *,"Running GEODESIC Solver ...."

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=RAW_INDEX(MIN_TRIAL_NODE)+1,RAW_INDEX(MIN_TRIAL_NODE+1)
        TIME_NEW=1000

        DO J=RAW_INDEX(COLUMN_INDEX(I))+1,RAW_INDEX(COLUMN_INDEX(I)+1)
          IF (STATUS_MASK(COLUMN_INDEX(J)) == "KNOWN") THEN 

            TIME_ITER=SEED_VALUE(COLUMN_INDEX(J))+CONNECTIVITY_WEIGHT(J)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
              TRACE_NODE_NUMBER=COLUMN_INDEX(J)
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .LT. SEED_VALUE(COLUMN_INDEX(I))) THEN
          SEED_VALUE(COLUMN_INDEX(I)) = TIME_NEW
          TRACE_NODE(COLUMN_INDEX(I)) = TRACE_NODE_NUMBER
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "KNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE = I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0 .AND. STATUS_MASK(I) .EQ. "SEED POINT" .AND. SEED_VALUE(I) .LT. MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO

999 ERRORSEXITS("SOLVE_PROBLEM_GEODESIC_CONNECTIVITY",ERR,ERROR)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_GEODESIC_CONNECTIVITY

  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_GEODESIC(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,&
  & SEED_VALUE,CONNECTIVITY_NUMBER,CONNECTIVITY_LIST,STATUS_MASK,TRACE_NODE)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    INTEGER(INTG), ALLOCATABLE :: TRACE_NODE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER,CHANGED_STATUS,MIN_TRIAL_NODE,TRIAL_STATUS,TRACE_NODE_NUMBER
    REAL(DP), DIMENSION(3) :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW,CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0

    !Start Program

    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

!   4444444444444444 MAIN LOOP 55555555555555555	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=1,CONNECTIVITY_NUMBER(MIN_TRIAL_NODE)
        TIME_NEW=1000

        DO J=1,CONNECTIVITY_NUMBER(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))
          IF (STATUS_MASK(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=(&
                             &/NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),1)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),2)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),3)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3)/)

            CONDUCTION_RATIO=SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)/&
                            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO/), SHAPE = (/3,3/))

!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),4),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),5),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),6),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),7),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),8),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),9)/&
                               &),SHAPE = (/3,3/))

            CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
            CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)

            CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J))+&
                      &SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
              TRACE_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .lt. SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))) THEN
          SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = TIME_NEW
          TRACE_NODE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))=TRACE_NODE_NUMBER
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "KNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000
      CHANGED_STATUS = 0

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0.AND.STATUS_MASK(I) .EQ. "SEED POINT".AND.SEED_VALUE(I).LT.MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


999 ERRORSEXITS("SOLVE_PROBLEM_GEODESIC",ERR,ERROR)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_GEODESIC

  
  !
  !================================================================================================================================
  !

  !>Calculates status mask for the local nodes.
  SUBROUTINE GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    !subroutine variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)

    !Argument variables
    INTEGER(INTG) :: Err !<The error code
    !    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string

    !Local variables
    INTEGER(INTG) :: I
        
!    ENTERS("GENERATE_STATUS_MASK",Err,Error,*999)
    
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (SEED_VALUE(I) .LT. 100.0_DP) THEN
        STATUS_MASK(I) = "SEED POINT"
      ELSE
        STATUS_MASK(I) = "UNKNOWN"
      ENDIF

    ENDDO

!    EXITS("GENERATE_STATUS_MASK")
!    RETURN
!999 ERRORSEXITS("GENERATE_STATUS_MASK",ERR,ERROR)
!    RETURN 1

  END SUBROUTINE GENERATE_STATUS_MASK


  !
  !================================================================================================================================
  !

  !>Calculates minimum and maximum value at array A.
  SUBROUTINE FIND_MINIMAX(A,N,MIN_VALUE,MAX_VALUE,Err)

    !Argument variables
    REAL(DP), ALLOCATABLE :: A(:)

    REAL(DP), INTENT(OUT) :: MIN_VALUE
    REAL(DP), INTENT(OUT) :: MAX_VALUE
    INTEGER(INTG), INTENT(IN)  :: N
    INTEGER(INTG) :: ERR !<The error code
    !    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I
        
!    ENTERS("FIND_MINIMAX",Err,Error,*999)
    
    IF(SIZE(A,1).GT.2) THEN
      MIN_VALUE = A(1)
      MAX_VALUE = A(1)
      DO I=2,N
        IF (MIN_VALUE .GT. A(I)) THEN
          MIN_VALUE = A(I)
        ENDIF
        IF (MAX_VALUE .LT. A(I)) THEN
          MAX_VALUE = A(I)
        ENDIF
      ENDDO
!    ELSE
!      CALL FlagError("Invalid matrix sizes.",Err)
    ENDIF

!    EXITS("FIND_MINIMAX")
!    RETURN
!999 ERRORSEXITS("FIND_MINIMAX",ERR,ERROR)
!    RETURN 1

  END SUBROUTINE FIND_MINIMAX

  !
  !================================================================================================================================
  !


  !>Calculates and returns the VECTOR-VECTOR-prouct of the double precision VECTOR A*B in C.
  SUBROUTINE VECTOR_VECTOR_PRODUCT(A,B,C,Err)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(3) !<The A VECTOR
    REAL(DP), INTENT(IN) :: B(3) !<The B VECTOR
    REAL(DP), INTENT(OUT) :: C !<On exit, the product SCALAR C=A*B
    INTEGER(INTG) :: ERR !<The error code
    !    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string
    !Local variables
        
!    ENTERS("VECTOR_VECTOR_PRODUCT",Err,Error,*999)
    
    IF(SIZE(A,1)==SIZE(B,1)) THEN
      SELECT CASE(SIZE(A,1))
        CASE(1)
          C=A(1)*B(1)
        CASE(2)
          C=A(1)*B(1)+A(2)*B(2)
        CASE(3)
          C=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
!        CASE DEFAULT
!          CALL FlagError("Invalid matrix size.",Err)
      END SELECT
!    ELSE
!      CALL FlagError("Invalid matrix sizes.",Err)
    ENDIF

!    EXITS("VECTOR_VECTOR_PRODUCT")
!    RETURN
!999 ERRORSEXITS("VECTOR_VECTOR_PRODUCT",ERR,ERROR)
!    RETURN 1

  END SUBROUTINE VECTOR_VECTOR_PRODUCT

  !
  !================================================================================================================================
  !


  !>to EXPORT output.
  SUBROUTINE POST_PROCESS_DATA(MATERIAL_BEHAVIOUR,OUTPUT_FILE_NAME,OUTPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,NODE_LIST,&
&CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,OUTPUT_FILE_FIELD_TITLE,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

    !subroutine variables
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)

    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES_PER_ELEMENT
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    CHARACTER (LEN=300) :: OUTPUT_FILE_NAME
    CHARACTER (LEN=300) :: OUTPUT_FILE_FIELD_TITLE
    CHARACTER (LEN=10)  :: OUTPUT_FILE_FORMAT
    CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
    INTEGER(INTG) :: Err
    
    !Local variables
    CHARACTER (LEN=300) :: STRING
    INTEGER(INTG) :: TEXT_LENGTH
    INTEGER(INTG) :: I, J

!    ENTERS("GENERATE_STATUS_MASK",Err,Error,*999)

!    IF (OUTPUT_FILE_FORMAT .NE. "TABC") THEN
!      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,">> OUTPUT FILE FORMAT ERROR: identify a correct output format",Err)
!    ENDIF
    

! in the case the OUTPUT is in TABC format (2)
    IF (OUTPUT_FILE_FORMAT .EQ. "TABC") THEN 

!      EXPORT NODE TABLE list
      TEXT_LENGTH = INDEX(OUTPUT_FILE_NAME,' ') - 1
      STRING = OUTPUT_FILE_NAME(1:TEXT_LENGTH)//".tabc"
!  INQUIRE(FILE=STRING, EXIST=ex)
!      PRINT *, STRING
      OPEN (12,FILE=STRING)
!      OPEN (12,FILE=OUTPUT_FILE_NAME)

      WRITE(12,*)"VARIABLES=""NODE"",""X"",""Y"",""Z"",""U"",""V"",""W"",""Speed function along fibers"", &
 &                ""Speed function in transverse direction"",""Time"""
      WRITE(12,*)"zone i=",TOTAL_NUMBER_OF_NODES," , DATAPACKING=POINT"

      DO I=1,TOTAL_NUMBER_OF_NODES
        WRITE(12,*) I,NODE_LIST(I,1),NODE_LIST(I,2),NODE_LIST(I,3),CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),&
            &CONDUCTIVITY_TENSOR(I,3),SPEED_FUNCTION_TABLE(I,1),&
            &SPEED_FUNCTION_TABLE(I,2),SPEED_FUNCTION_TABLE(I,3),&
            &SEED_VALUE(I)
      ENDDO
!      EXPORT NODE CONNECTIVITY list
      WRITE(12,*) "Connectivity"
      DO I=1,TOTAL_NUMBER_OF_NODES
        WRITE(12,*) I,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
!        PRINT *,STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
      ENDDO
      CLOSE(12)

    ENDIF

! in the case the OUTPUT is in VTK TETrahedral format (2)
    IF (OUTPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

!      rename to VTK and OPEN the file
      TEXT_LENGTH = INDEX(OUTPUT_FILE_NAME,' ') - 1
      STRING = OUTPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"
      OPEN (12,FILE=STRING)

!      export HEADER list
      WRITE(12,'(A)')"# vtk DataFile Version 3.0"
      WRITE(12,'(A)')"vtk output"
      WRITE(12,'(A)')"ASCII"
      WRITE(12,'(A)')"DATASET UNSTRUCTURED_GRID"
      WRITE(12,'(A,I8,A6)')"POINTS",TOTAL_NUMBER_OF_NODES,"float"

!      export NODAL POSITION list
      DO I=1,TOTAL_NUMBER_OF_NODES
        WRITE(12,*) (NODE_LIST(I,J),J=1,3)
      ENDDO

!      export ELEMENT list
      WRITE(12,'(A,I8,I8)')"CELLS",TOTAL_NUMBER_OF_ELEMENTS,TOTAL_NUMBER_OF_ELEMENTS*(NUMBER_OF_NODES_PER_ELEMENT+1)
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        WRITE(12,*) NUMBER_OF_NODES_PER_ELEMENT,((ELEMENT_LIST(I,J)-1),J=1,NUMBER_OF_NODES_PER_ELEMENT)
      ENDDO

!      export ELEMENT TYPE list
      WRITE(12,'(A,I8)')"CELL_TYPES",TOTAL_NUMBER_OF_ELEMENTS
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        WRITE(12,'(A)') "10"
      ENDDO

!      export CELL and POINT DATA list
      WRITE(12,'(A,I8)')"CELL_DATA",TOTAL_NUMBER_OF_ELEMENTS
      WRITE(12,'(A,I8)')"POINT_DATA",TOTAL_NUMBER_OF_NODES

!      export FIELD information
      WRITE(12,'(A,A)')"FIELD number"," 1"
      WRITE(12,'(A,I3,I8,A6)')OUTPUT_FILE_FIELD_TITLE,1,TOTAL_NUMBER_OF_NODES,"float"
      DO I=1,TOTAL_NUMBER_OF_NODES
        WRITE(12,'(F15.10)') SEED_VALUE(I)
      ENDDO

!      export VECTORS information
      WRITE(12,'(A,A,A6)') "VECTORS ","fiber_vector","float"
      DO I=1,TOTAL_NUMBER_OF_NODES
        WRITE(12,'(3F8.5)') (CONDUCTIVITY_TENSOR(I,J),J=1,3)
      ENDDO

      CLOSE(12)

    ENDIF

!   FILE="cmgui"
!   METHOD="FORTRAN"

!   EXPORT_FIELD=.TRUE.
!   IF(EXPORT_FIELD) THEN
!     WRITE(*,*)'Now export fields...'
!     CALL FLUID_MECHANICS_IO_WRITE_CMGUI(REGION,FILE,Err)
!     CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, Err)  
!     CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, Err)
!     WRITE(*,*)'All fields exported...'
!   ENDIF

!    EXITS("GENERATE_STATUS_MASK")
!    RETURN
!999 ERRORSEXITS("GENERATE_STATUS_MASK",ERR,ERROR)
!    RETURN 1

  END SUBROUTINE POST_PROCESS_DATA

 
END MODULE HAMILTON_JACOBI_EQUATIONS_ROUTINES
