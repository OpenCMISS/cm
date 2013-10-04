!> \file
!> \author Chris Bradley
!> \brief This module handles all Laplace equations routines.
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

!>This module handles all Laplace equations routines.
MODULE LAPLACE_EQUATIONS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES  
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
  USE MATHS  
  USE MATRIX_VECTOR
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

!!MERGE: move

  PUBLIC LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE,LAPLACE_EQUATION_EQUATIONS_SET_SETUP,LAPLACE_EQUATION_ANALYTIC_CALCULATE, &
    & LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET,LaplaceEquation_EquationsSetSpecificationSet, &
    & LaplaceEquation_ProblemSpecificationSet,LAPLACE_EQUATION_PROBLEM_SETUP

CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE LAPLACE_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

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
    
    CALL ENTERS("LAPLACE_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR,*999)

    NULLIFY(GEOMETRIC_PARAMETERS)

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
                                CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)
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
                                        & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
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
                                CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2)
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
                                        & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
                                    CASE(GLOBAL_DERIV_S1)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
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
                                CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)
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
                                        & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S3)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
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
                                CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2)
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
                                        & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP !!TODO
                                    CASE(GLOBAL_DERIV_S1)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S3)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S3)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S2_S3)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(GLOBAL_DERIV_S1_S2_S3)
                                      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
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
    
    CALL EXITS("LAPLACE_EQUATION_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_ANALYTIC_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Laplace equation finite element equations set.
  SUBROUTINE LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,mi,ms,nh,nhs,ni,ns,nj,a,b,i,k,h
    REAL(DP) :: CONDUCTIVITY_MATERIAL(3,3),CONDUCTIVITY(3,3),CONDUCTIVITY_TEMP(3,3),RWG,SUM,SUM1,SUM2,PGMSI(3),PGNSI(3),K_VALUE(3)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS,INDEPENDENT_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD,MATERIALS_FIELD,FIBRE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERP_POINT_METRICS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_XI
    REAL(DP) :: DNUDXI(3,3),DXIDNU(3,3)
    
#ifdef TAUPROF
    CHARACTER(26) :: CVAR
    INTEGER :: GAUSS_POINT_LOOP_PHASE(2) = [ 0, 0 ]
    SAVE GAUSS_POINT_LOOP_PHASE
#endif

    CALL ENTERS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Laplace is put in.
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
              ENDDO !ms
            ENDDO !mh
          ENDIF
          
        CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
!!TODO: store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          FIBRE_FIELD=>EQUATIONS%INTERPOLATION%FIBRE_FIELD
          
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
                  
          NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
          NUMBER_OF_XI=DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY% &
            & ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_XI
                    
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            
          GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          GEOMETRIC_INTERP_POINT_METRICS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
          
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
#ifdef TAUPROF
              WRITE (CVAR,'(a17,i2)') 'Gauss Point Loop ',ng
              CALL TAU_PHASE_CREATE_DYNAMIC(GAUSS_POINT_LOOP_PHASE,CVAR)
              CALL TAU_PHASE_START(GAUSS_POINT_LOOP_PHASE)
#endif
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERP_POINT_METRICS, &
              & ERR,ERROR,*999)
              
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN*QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
              
            !conductivity in material coordinates 
            CONDUCTIVITY_MATERIAL=0.0_DP
            IF(NUMBER_OF_DIMENSIONS==2) THEN
              CONDUCTIVITY_MATERIAL(1,1)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1)
              CONDUCTIVITY_MATERIAL(2,2)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,1)
              CONDUCTIVITY_MATERIAL(1,2)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,1)
              CONDUCTIVITY_MATERIAL(2,1)=CONDUCTIVITY_MATERIAL(1,2)
            ELSE
              CONDUCTIVITY_MATERIAL(1,1)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1)
              CONDUCTIVITY_MATERIAL(2,2)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,1)
              CONDUCTIVITY_MATERIAL(3,3)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,1)
              CONDUCTIVITY_MATERIAL(1,2)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(4,1)
              CONDUCTIVITY_MATERIAL(2,1)=CONDUCTIVITY_MATERIAL(1,2)
              CONDUCTIVITY_MATERIAL(2,3)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(5,1)
              CONDUCTIVITY_MATERIAL(3,2)=CONDUCTIVITY_MATERIAL(2,3)
              CONDUCTIVITY_MATERIAL(1,3)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(6,1)
              CONDUCTIVITY_MATERIAL(3,1)=CONDUCTIVITY_MATERIAL(1,3)
            ENDIF
  
            !rotate the conductivity from material coordinates into xi-space to get the effective conductivity
            CALL CoordinateMaterialSystemCalculate(GEOMETRIC_INTERP_POINT_METRICS,FIBRE_INTERPOLATED_POINT, &
              & DNUDXI,DXIDNU,ERR,ERROR,*999)

            CALL MATRIX_PRODUCT(DNUDXI,CONDUCTIVITY_MATERIAL,CONDUCTIVITY_TEMP,ERR,ERROR,*999)
            CALL MATRIX_PRODUCT(CONDUCTIVITY_TEMP,DXIDNU,CONDUCTIVITY,ERR,ERROR,*999)
            
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                !Loop over field components
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    !Loop over element columns
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                      ENDDO !ni
                      
                      SUM=0.0_DP
                      DO i=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        DO k=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          DO h=1,DEPENDENT_BASIS%NUMBER_OF_XI
                            SUM=SUM+CONDUCTIVITY(i,k)*PGNSI(k)*PGMSI(h)*GEOMETRIC_INTERP_POINT_METRICS%GU(i,h)
                          ENDDO !h
                        ENDDO !k
                      ENDDO !i
                        
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
              ENDDO !ms
            ENDDO !mh
          ENDIF

        CASE(EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Laplace is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          INDEPENDENT_BASIS=>INDEPENDENT_FIELD%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & INDEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

            K_VALUE(1)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
            K_VALUE(2)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
            IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
              K_VALUE(3)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
            END IF 


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

                      IF(nh==mh) THEN 
                        SUM=0.0_DP
                        DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                            SUM=SUM+PGMSI(mi)*PGNSI(ni)*EQUATIONS%INTERPOLATION% &
                              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%GU(mi,ni)*K_VALUE(mh)
                          ENDDO !ni
                        ENDDO !mi
                        EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                      ENDIF

                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
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
              ENDDO !ms
            ENDDO !mh
          ENDIF

        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a Laplace equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the moving mesh Laplace equation.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_MOVING_MESH_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,GEOMETRIC_COMPONENT_NUMBER,MATERIAL_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,NUMBER_OF_DIMENSIONS,I,MATERIAL_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATION_SET_MOVING_MESH_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a moving mesh Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"Phi",ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
                & ERR,ERROR,*999)

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

              !calculate number of components with one component for each dimension and one for pressure
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)

              DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              END DO


              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)

              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                END DO
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
                !Other solutions not defined yet
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
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, & 
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)

              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)

              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
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
              & " is invalid for a moving mesh Laplace equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
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
              INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              !calculate number of components with one component for each dimension
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              DO I=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999) 
              END DO
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              !calculate number of components with one component for each dimension and one for pressure
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
            ENDIF    
          !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
            CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
              & FIELD_BOUNDARY_SET_TYPE,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=1!X
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS=1!Y

            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                !start field creation with name 'MATERIAL_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET% & 
                  & MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, & 
                  & ERR,ERROR,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                  & MATERIAL_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & MATERIAL_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
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
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
              ENDIF              
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1.0_DP,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a moving mesh Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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
              & " is invalid for a moving mesh Laplace equation."
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
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a moving mesh Laplace equation."
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
              & " is invalid for a moving mesh Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
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
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
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
              & " is invalid for a moving mesh Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a moving mesh Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a moving mesh Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_MOVING_MESH_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_MOVING_MESH_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_MOVING_MESH_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_MOVING_MESH_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Laplace equation type of a classical field equations set class.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Laplace equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
        CALL LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
        CALL LAPLACE_EQUATION_EQUATIONS_SET_MOVING_MESH_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
        CALL LAPLACE_EQUATION_EQUATIONS_SET_GENERALISED_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Laplace equation type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Laplace equation type of an classical field equations set class.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)        
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
      CASE(EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)        
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
      CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)        
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
          & " is not valid for a Laplace equation type of an classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Laplace equation type of a classical field equations set class.
  SUBROUTINE LaplaceEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    CALL Enters("LaplaceEquation_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE, &
          & EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE, &
          & EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a Laplace equation type of a classical field equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_LAPLACE_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    CALL Exits("LaplaceEquation_EquationsSetSpecificationSet")
    RETURN
999 CALL Errors("LaplaceEquation_EquationsSetSpecificationSet",err,error)
    CALL Exits("LaplaceEquation_EquationsSetSpecificationSet")
    RETURN 1
  END SUBROUTINE LaplaceEquation_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the standard Laplace equation.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATION_SET_STANDARD_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
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
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
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
              & " is invalid for a standard Laplace equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
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
              & " is invalid for a standard Laplace equation."
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
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a standard Laplace equation."
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
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
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
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
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
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a standard Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the generalised Laplace equation.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_GENERALISED_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS,NUMBER_OF_MATERIALS_COMPONENTS
    INTEGER(INTG) :: component_idx
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATIONS_SET_GENERALISED_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        !INITIAL
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a generalised Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !GEOMETRY
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        !DEPENDENT
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          !start action
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
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
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
          !finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a generalised Laplace equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !MATERIAL
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          !start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
              !Create the auto created materials field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%MATERIALS% &
                & MATERIALS_FIELD,ERR,ERROR,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,"Materials Field",ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                & ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"conductivity",ERR, &
                & ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              IF(NUMBER_OF_DIMENSIONS==1) THEN
                NUMBER_OF_MATERIALS_COMPONENTS=1
              ELSEIF(NUMBER_OF_DIMENSIONS==2) THEN
                NUMBER_OF_MATERIALS_COMPONENTS=3
              ELSEIF(NUMBER_OF_DIMENSIONS==3) THEN
                NUMBER_OF_MATERIALS_COMPONENTS=6
              ENDIF
              !Set the number of materials components
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)        
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              !Default the materials components to the first geometric component
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              DO component_idx=1,NUMBER_OF_MATERIALS_COMPONENTS                
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              ENDDO !components_idx
              
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_MATERIALS_COMPONENTS                
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
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
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              IF(NUMBER_OF_DIMENSIONS==1) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,ERR, &
                  & ERROR,*999)
              ELSEIF(NUMBER_OF_DIMENSIONS==2) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,3,ERR, &
                  & ERROR,*999)
              ELSEIF(NUMBER_OF_DIMENSIONS==3) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,6,ERR, &
                  & ERROR,*999)
              ENDIF              
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !do nothing
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
          !finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a generalised Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !SOURCE
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a generalised Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !ANALYTIC
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          !start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1
                  CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2)
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
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a generalised Laplace equation."
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
          !finish action
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
              & " is invalid for a generalised Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !EQUATIONS
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          !start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          !finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
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
              & " is invalid for a generalised Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a generalised Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a generalised Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_GENERALISED_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_GENERALISED_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_GENERALISED_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_GENERALISED_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Laplace problem.
  SUBROUTINE LAPLACE_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Laplace equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Laplace problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)
        CALL LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_LAPLACE_SUBTYPE)
        CALL LAPLACE_EQUATION_PROBLEM_GENERALISED_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Laplace equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Laplace equation type.
  SUBROUTINE LaplaceEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    CALL Enters("LaplaceEquation_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE, &
            & PROBLEM_GENERALISED_LAPLACE_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Laplace equation type of a classical field problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_LAPLACE_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Laplace problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    CALL Exits("LaplaceEquation_ProblemSpecificationSet")
    RETURN
999 CALL Errors("LaplaceEquation_ProblemSpecificationSet",err,error)
    CALL Exits("LaplaceEquation_ProblemSpecificationSet")
    RETURN 1
  END SUBROUTINE LaplaceEquation_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the standard Laplace equations problem.
  SUBROUTINE LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
    CALL ENTERS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Laplace problem.",err,error,*999)
      END IF
      IF(PROBLEM%SPECIFICATION(3)==PROBLEM_STANDARD_LAPLACE_SUBTYPE) THEN
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
              & " is invalid for a standard Laplace equation."
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
              & " is invalid for a standard Laplace equation."
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
                & " is invalid for a standard Laplace equation."
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
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a standard Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the generalised Laplace equations problem.
  SUBROUTINE LAPLACE_EQUATION_PROBLEM_GENERALISED_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
    CALL ENTERS("LAPLACE_EQUATION_PROBLEM_GENERALISED_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Laplace problem.",err,error,*999)
      END IF
      IF(PROBLEM%SPECIFICATION(3)==PROBLEM_GENERALISED_LAPLACE_SUBTYPE) THEN
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
              & " is invalid for a generalised Laplace equation."
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
              & " is invalid for a generalised Laplace equation."
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
                & " is invalid for a generalised Laplace equation."
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
              & " is invalid for a generalised Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a generalised Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a generalised Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_GENERALISED_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_PROBLEM_GENERALISED_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_GENERALISED_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_PROBLEM_GENERALISED_SETUP

  !
  !================================================================================================================================
  !
  
END MODULE LAPLACE_EQUATIONS_ROUTINES
