!> \file
!> \author Sebastian Krittian
!> \brief This module handles all fitting routines.
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

!>This module handles all Galerkin projection routines.
MODULE FITTING_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE DARCY_EQUATIONS_ROUTINES, ONLY: idebug1
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
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

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

!!MERGE: move

  PUBLIC FITTING_EQUATIONS_SET_SETUP
  PUBLIC FITTING_EQUATIONS_SET_VECTORDATA_SETUP
  PUBLIC FITTING_EQUATIONS_SET_SUBTYPE_SET
  PUBLIC FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC FITTING_PROBLEM_SETUP
  PUBLIC FITTING_PROBLEM_SUBTYPE_SET

  PUBLIC FITTING_FINITE_ELEMENT_CALCULATE
  PUBLIC FITTING_EQUATIONS_SET_CLASS_TYPE_SET
  PUBLIC FITTING_EQUATIONS_SET_CLASS_TYPE_GET
  PUBLIC FITTING_PROBLEM_CLASS_TYPE_SET
  PUBLIC FITTING_PROBLEM_CLASS_TYPE_GET

  PUBLIC FITTING_PRE_SOLVE
  PUBLIC FITTING_POST_SOLVE
  PUBLIC FITTING_PRE_SOLVE_UPDATE_INPUT_DATA

CONTAINS

  !
  !================================================================================================================================
  !


! ! !   !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
! ! !   SUBROUTINE FITTING_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*)
! ! ! 
! ! !     !Argument variables
! ! !     TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET 
! ! !     INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
! ! !     TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
! ! !     !Local Variables
! ! !     INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type
! ! !     REAL(DP) :: VALUE,X(3)
! ! !     REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
! ! !     TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
! ! !     TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
! ! !     TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
! ! !     TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
! ! !     TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
! ! !     TYPE(VARYING_STRING) :: LOCAL_ERROR    
! ! !     
#if DEBUG
! ! !     CALL ENTERS("FITTING_ANALYTIC_CALCULATE",ERR,ERROR,*999)
#endif
! ! !     
! ! !     IF(ASSOCIATED(EQUATIONS_SET)) THEN
! ! !       IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
! ! !         DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
! ! !         IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
! ! !           GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
! ! !           IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
! ! !             CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
! ! !             NULLIFY(GEOMETRIC_VARIABLE)
! ! !             CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
! ! !             CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
! ! !               & ERR,ERROR,*999)
! ! !             NULLIFY(BOUNDARY_CONDITIONS)
! ! !             CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
! ! !             DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
! ! !               variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
! ! !               FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
! ! !               IF(ASSOCIATED(FIELD_VARIABLE)) THEN
! ! !                 CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
! ! !                 DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
! ! !                   IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
! ! !                     DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
! ! !                     IF(ASSOCIATED(DOMAIN)) THEN
! ! !                       IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
! ! !                         DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
! ! !                         IF(ASSOCIATED(DOMAIN_NODES)) THEN
! ! !                           !Loop over the local nodes excluding the ghosts.
! ! !                           DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
! ! !                             !!TODO \todo We should interpolate the geometric field here and the node position.
! ! !                             DO dim_idx=1,NUMBER_OF_DIMENSIONS
! ! !                               local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
! ! !                               X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
! ! !                             ENDDO !dim_idx
! ! !                             !Loop over the derivatives
! ! !                             DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
! ! !                               SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
! ! !                               CASE(EQUATIONS_SET_FITTING_TWO_DIM_1)
! ! !                                 !u=x^2+2.x.y-y^2
! ! !                                 SELECT CASE(variable_type)
! ! !                                 CASE(FIELD_U_VARIABLE_TYPE)
! ! !                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=X(1)*X(1)-2.0_DP*X(1)*X(2)-X(2)*X(2)
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     VALUE=2.0_DP*X(1)+2.0_DP*X(2)
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     VALUE=2.0_DP*X(1)-2.0_DP*X(2)
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     VALUE=2.0_DP
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE(FIELD_DELUDELN_VARIABLE_TYPE)
! ! !                                  SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=0.0_DP !!TODO
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE DEFAULT
! ! !                                   LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
! ! !                                     & " is invalid."
! ! !                                   CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                 END SELECT                                
! ! !                               CASE(EQUATIONS_SET_FITTING_TWO_DIM_2)
! ! !                                 !u=cos(x).cosh(y)
! ! !                                 SELECT CASE(variable_type)
! ! !                                 CASE(FIELD_U_VARIABLE_TYPE)
! ! !                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=COS(X(1))*COSH(X(2))
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     VALUE=-SIN(X(1))*COSH(X(2))
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     VALUE=COS(X(1))*SINH(X(2))
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     VALUE=-SIN(X(1))*SINH(X(2))
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE(FIELD_DELUDELN_VARIABLE_TYPE)
! ! !                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=0.0_DP !!TODO
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE DEFAULT
! ! !                                   LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
! ! !                                     & " is invalid."
! ! !                                   CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                 END SELECT                                
! ! !                               CASE(EQUATIONS_SET_FITTING_THREE_DIM_1)
! ! !                                 !u=x^2+y^2-2.z^2
! ! !                                 SELECT CASE(variable_type)
! ! !                                 CASE(FIELD_U_VARIABLE_TYPE)
! ! !                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=X(1)*X(1)+X(2)*X(2)-2.0_DP*X(3)*X(3)
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     VALUE=2.0_DP*X(1)
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     VALUE=2.0_DP*X(2)
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     VALUE=0.0_DP
! ! !                                   CASE(GLOBAL_DERIV_S3)
! ! !                                     VALUE=-4.0_DP*X(3)
! ! !                                   CASE(GLOBAL_DERIV_S1_S3)
! ! !                                     VALUE=0.0_DP
! ! !                                   CASE(GLOBAL_DERIV_S2_S3)
! ! !                                     VALUE=0.0_DP
! ! !                                   CASE(GLOBAL_DERIV_S1_S2_S3)
! ! !                                     VALUE=0.0_DP
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE(FIELD_DELUDELN_VARIABLE_TYPE)
! ! !                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=0.0_DP !!TODO
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S3)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S1_S3)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S2_S3)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S1_S2_S3)
! ! !                                     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE DEFAULT
! ! !                                   LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
! ! !                                     & " is invalid."
! ! !                                   CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                 END SELECT                                
! ! !                               CASE(EQUATIONS_SET_FITTING_THREE_DIM_2)
! ! !                                 !u=cos(x).cosh(y).z
! ! !                                 SELECT CASE(variable_type)
! ! !                                 CASE(FIELD_U_VARIABLE_TYPE)
! ! !                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=COS(X(1))*COSH(X(2))*X(3)
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     VALUE=-SIN(X(1))*COSH(X(2))*X(3)
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     VALUE=COS(X(1))*SINH(X(2))*X(3)
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     VALUE=-SIN(X(1))*SINH(X(2))*X(3)
! ! !                                   CASE(GLOBAL_DERIV_S3)
! ! !                                     VALUE=COS(X(1))*COSH(X(2))
! ! !                                   CASE(GLOBAL_DERIV_S1_S3)
! ! !                                     VALUE=-SIN(X(1))*COSH(X(2))
! ! !                                   CASE(GLOBAL_DERIV_S2_S3)
! ! !                                     VALUE=COS(X(1))*SINH(X(2))
! ! !                                   CASE(GLOBAL_DERIV_S1_S2_S3)
! ! !                                     VALUE=-SIN(X(1))*SINH(X(2))
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE(FIELD_DELUDELN_VARIABLE_TYPE)
! ! !                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
! ! !                                   CASE(NO_GLOBAL_DERIV)
! ! !                                     VALUE=0.0_DP !!TODO
! ! !                                   CASE(GLOBAL_DERIV_S1)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S2)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
! ! !                                   CASE(GLOBAL_DERIV_S1_S2)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S3)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S1_S3)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S2_S3)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE(GLOBAL_DERIV_S1_S2_S3)
! ! !                                     !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
! ! !                                   CASE DEFAULT
! ! !                                     LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
! ! !                                       DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
! ! !                                       & " is invalid."
! ! !                                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                   END SELECT
! ! !                                 CASE DEFAULT
! ! !                                   LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
! ! !                                     & " is invalid."
! ! !                                   CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                                 END SELECT                                
! ! !                               CASE DEFAULT
! ! !                                 LOCAL_ERROR="The analytic function type of "// &
! ! !                                   & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                                   & " is invalid."
! ! !                                 CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                               END SELECT
! ! !                               local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
! ! !                                 & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
! ! !                               CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
! ! !                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
! ! !                               IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
! ! !                                 IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
! ! !                                   !If we are a boundary node then set the analytic value on the boundary
! ! !                                   CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
! ! !                                     & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
! ! !                                 ENDIF
! ! !                               ENDIF
! ! !                             ENDDO !deriv_idx
! ! !                           ENDDO !node_idx
! ! !                         ELSE
! ! !                           CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
! ! !                         ENDIF
! ! !                       ELSE
! ! !                         CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
! ! !                       ENDIF
! ! !                     ELSE
! ! !                       CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                   ELSE
! ! !                     CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
! ! !                   ENDIF
! ! !                 ENDDO !component_idx
! ! !                 CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
! ! !                   & ERR,ERROR,*999)
! ! !                 CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
! ! !                   & ERR,ERROR,*999)
! ! !               ELSE
! ! !                 CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
! ! !               ENDIF
! ! ! 
! ! !             ENDDO !variable_idx
! ! !             CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
! ! !             CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
! ! !               & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
! ! !           ELSE
! ! !             CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
! ! !           ENDIF            
! ! !         ELSE
! ! !           CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
! ! !         ENDIF
! ! !       ELSE
! ! !         CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
! ! !       ENDIF
! ! !     ELSE
! ! !       CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
! ! !     ENDIF
! ! !     
#if DEBUG
! ! !     CALL EXITS("FITTING_ANALYTIC_CALCULATE")
#endif
! ! !     RETURN
! ! ! 999 CALL ERRORS("FITTING_ANALYTIC_CALCULATE",ERR,ERROR)
#if DEBUG
! ! !     CALL EXITS("FITTING_ANALYTIC_CALCULATE")
#endif
! ! !     RETURN 1
! ! !   END SUBROUTINE FITTING_ANALYTIC_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Galerkin projection finite element equations set.
  SUBROUTINE FITTING_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: FIELD_VAR_TYPE,ng,mh,mhs,ms,nh,nhs,ns,mi,ni
    REAL(DP) :: RWG,SUM,jacobianGaussWeight
    REAL(DP) :: PGM,PGN,PGMSI(3),PGNSI(3)
    REAL(DP) :: U_VALUE(3)
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS,SOURCE_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD,DEPENDENT_FIELD,MATERIALS_FIELD,SOURCE_FIELD
    TYPE(FIELD_TYPE), POINTER :: independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: mappingVariable
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: MATERIALS_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: REFERENCE_GEOMETRIC_INTERPOLATED_POINT
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    TYPE(VARYING_STRING) :: localError

    REAL(DP), POINTER :: independentVectorParameters(:),independentWeightParameters(:)
    REAL(DP), ALLOCATABLE :: projectionXi(:)
    REAL(DP):: POROSITY_0, POROSITY, PERM_OVER_VIS_PARAM_0, PERM_OVER_VIS_PARAM,TAU_PARAM,KAPPA_PARAM
    REAL(DP):: tension,curvature
    REAL(DP):: MATERIAL_FACT
    REAL(DP):: DXDY(3,3), DXDXI(3,3), DYDXI(3,3), DXIDY(3,3), DXI_DX(3,3)
    REAL(DP):: Jxy, Jyxi
    REAL(DP):: dataPointWeight,dataPointVector(3)
    INTEGER(INTG) :: derivative_idx, component_idx, xi_idx, NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: dataPointIdx,dataPointUserNumber,dataPointLocalNumber,dataPointGlobalNumber
    INTEGER(INTG) :: numberOfXi
    INTEGER(INTG) :: componentIdx
    INTEGER(INTG) :: variableType,localDof

    INTEGER(INTG) NDOFS
    INTEGER(INTG) MESH_COMPONENT1,MESH_COMPONENT2


    
#if DEBUG
    CALL ENTERS("FITTING_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)
#endif

    NULLIFY(DEPENDENT_BASIS,GEOMETRIC_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(LINEAR_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(LINEAR_MATRICES)
    NULLIFY(RHS_VECTOR)
    NULLIFY(EQUATIONS_MATRIX)
    NULLIFY(DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD)
    NULLIFY(dataPoints)
    NULLIFY(dataProjection)
    NULLIFY(decompositionTopology)
    NULLIFY(independentField)
    NULLIFY(independentVectorParameters)
    NULLIFY(independentWeightParameters)
    NULLIFY(fieldVariable)
    NULLIFY(mappingVariable)
    NULLIFY(QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT)

    dataPointVector = 0.0_DP

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)

        CASE(EquationsSet_DataPointVectorStaticFittingSubtype, &
          &  EquationsSet_DataPointVectorQuasistaticFittingSubtype)
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          independentField=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD
          dataProjection=>independentField%dataProjection
          IF(.NOT.ASSOCIATED(dataProjection)) THEN
            localError="Data projection is not associated on independent field."
            CALL FLAG_ERROR(localError,err,error,*999)
          ENDIF
          decompositionTopology=>independentField%decomposition%topology
          IF(ASSOCIATED(decompositionTopology)) THEN
            dataPoints=>decompositionTopology%dataPoints
            IF(.NOT.ASSOCIATED(dataPoints)) THEN
              localError="Data points are not associated on the decomposition topology of the independent field."
              CALL FLAG_ERROR(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Decomposition topology is not associated on the independent field."
            CALL FLAG_ERROR(localError,err,error,*999)
          ENDIF
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          mappingVariable=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=mappingVariable%VARIABLE_TYPE
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
          numberOfXi = DEPENDENT_BASIS%NUMBER_OF_XI
          ALLOCATE(projectionXi(numberOfXi))
          projectionXi=0.0_DP
          ! Get data point vector parameters
          CALL FIELD_PARAMETER_SET_DATA_GET(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & independentVectorParameters,err,error,*999)
          ! Get data point weight parameters
          CALL FIELD_PARAMETER_SET_DATA_GET(independentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & independentWeightParameters,err,error,*999)

          !===========================================
          ! D a t a   P o i n t   V e c t o r    F i t
          !===========================================
          !Loop over data points
          DO dataPointIdx=1,dataPoints%elementDataPoint(ELEMENT_NUMBER)%numberOfProjectedData
            dataPointUserNumber = dataPoints%elementDataPoint(ELEMENT_NUMBER)%dataIndices(dataPointIdx)%userNumber
            dataPointLocalNumber = dataPoints%elementDataPoint(ELEMENT_NUMBER)%dataIndices(dataPointIdx)%localNumber
            dataPointGlobalNumber = dataPoints%elementDataPoint(ELEMENT_NUMBER)%dataIndices(dataPointIdx)%globalNumber
            ! Need to use global number to get the correct projection results
            projectionXi = dataProjection%data_projection_results(dataPointGlobalNumber)%xi
            CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,projectionXi,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,projectionXi,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

            ! Get data point vector value
            variableType=independentField%VARIABLES(1)%VARIABLE_TYPE
            fieldVariable=>independentField%VARIABLE_TYPE_MAP(variableType)%PTR
            DO componentIdx=1,numberOfXi
              localDof=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
               & DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointLocalNumber)
              dataPointVector(componentIdx)=independentVectorParameters(localDof)
            ENDDO

            variableType=independentField%VARIABLES(2)%VARIABLE_TYPE
            fieldVariable=>independentField%VARIABLE_TYPE_MAP(variableType)%PTR
            localDof=fieldVariable%COMPONENTS(1)%PARAM_TO_DOF_MAP% &
             & DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointLocalNumber)
            dataPointWeight=independentWeightParameters(localDof)

            mhs=0          
            !Loop over element rows
            DO mh=1,mappingVariable%NUMBER_OF_COMPONENTS
              MESH_COMPONENT1=mappingVariable%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                PGM=BASIS_EVALUATE_XI(DEPENDENT_BASIS1,ms,NO_PART_DERIV,projectionXi,err,error)
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,mappingVariable%NUMBER_OF_COMPONENTS
                    MESH_COMPONENT2=mappingVariable%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                    DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                      & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                    DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      PGN=BASIS_EVALUATE_XI(DEPENDENT_BASIS2,ns,NO_PART_DERIV,projectionXi,err,error)
                      SUM=0.0_DP
                      IF(mh==nh) THEN
                        SUM = SUM + PGM * PGN * dataPointWeight
                      ENDIF
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                SUM=0.0_DP
                IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                  SUM = SUM + PGM*dataPointVector(mh)*dataPointWeight
                  RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs) + SUM
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDDO !dataPointIdx

          !Restore data point vector parameters
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & independentVectorParameters,err,error,*999)
          !Restore data point weight parameters
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(independentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & independentWeightParameters,err,error,*999)
            
          !===========================================
          ! S o b e l o v   S m o o t h i n g 
          !===========================================
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            TAU_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
            KAPPA_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
            !Loop over field components
            jacobianGaussWeight=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

            mhs=0          
            DO mh=1,mappingVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
              MESH_COMPONENT1=mappingVariable%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                 !Loop over element columns
                  DO nh=1,mappingVariable%NUMBER_OF_COMPONENTS
                  MESH_COMPONENT2=mappingVariable%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP &
                    & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                    DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      SUM = 0.0_DP

                      !Calculate sobelov surface tension and curvature smoothing terms
                      tension = TAU_PARAM*2.0_DP* ( &
                        & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1,ng)* &
                        & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1,ng))
                      curvature = KAPPA_PARAM*2.0_DP* ( &
                        & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S1,ng)* &
                        & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,ng))

                      IF(mappingVariable%NUMBER_OF_COMPONENTS > 1) THEN
                        tension = tension + TAU_PARAM*2.0_DP* ( &
                          & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2,ng)* &
                          & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2,ng))
                        curvature = curvature + KAPPA_PARAM*2.0_DP* ( &
                          & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S2,ng)* &
                          & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,ng) + &
                          & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S2,ng)* &
                          & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,ng))

                        IF(mappingVariable%NUMBER_OF_COMPONENTS > 2) THEN
                          tension = tension + TAU_PARAM*2.0_DP* ( &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S3,ng))
                          curvature = curvature + KAPPA_PARAM*2.0_DP* ( &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S3_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,ng))
                        ENDIF ! 3D
                      ENDIF ! 2 or 3D

                      ! Add in smoothing terms to the element matrix
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) = &
                        & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) + (tension + curvature) * jacobianGaussWeight

                    ENDDO !ns
                  ENDDO !nh
                ENDIF ! update matrix
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            mhs=0          
            DO mh=1,mappingVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,mappingVariable%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                  RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                ENDIF 
              ENDDO !ms
            ENDDO !mh
          ENDIF

        CASE(EQUATIONS_SET_STANDARD_DATA_FITTING_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Galerkin projection is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          fieldVariable=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=fieldVariable%VARIABLE_TYPE
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
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
            DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns ???
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1

                      PGM=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                      PGN=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                      SUM = 0.0_DP
                      IF(mh==nh) THEN 
                        SUM = SUM + PGM * PGN
                      ENDIF
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) = &
                        & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) + SUM * RWG
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
            DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
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

        CASE(EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Galerkin projection is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD

          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          fieldVariable=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=fieldVariable%VARIABLE_TYPE

          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

          !--- Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS

            !--- Interpolation of Reference Geometry
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_INITIAL_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            REFERENCE_GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & REFERENCE_GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            !--- Retrieve local map DYDXI
            DO component_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
              DO xi_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7      
                DYDXI(component_idx,xi_idx)=REFERENCE_GEOMETRIC_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dy/dxi (y = referential)
              ENDDO
            ENDDO

            !--- Interpolation of (actual) Geometry and Metrics
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI, &
              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            !--- Retrieve local map DXDXI
            DO component_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
              DO xi_idx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7      
                DXDXI(component_idx,xi_idx)=GEOMETRIC_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dx/dxi
              ENDDO
            ENDDO

            !--- Compute deformation gradient tensor DXDY and its Jacobian Jxy
            CALL INVERT(DYDXI,DXIDY,Jyxi,ERR,ERROR,*999) !dy/dxi -> dxi/dy 
            CALL MATRIX_PRODUCT(DXDXI,DXIDY,DXDY,ERR,ERROR,*999) !dx/dxi * dxi/dy = dx/dy (deformation gradient tensor, F)
            Jxy=DETERMINANT(DXDY,ERR,ERROR)

            !--- Interpolation of Materials Field
            MATERIALS_INTERPOLATED_POINT => EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)

            !--- Retrieve reference material parameters:
            POROSITY_0            = MATERIALS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
            PERM_OVER_VIS_PARAM_0 = MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV)

            !--- Material dependence on structural deformation
            IF( ABS(Jxy) > 1.0E-10_DP ) THEN
              POROSITY = 1.0_DP - ( 1.0_DP - POROSITY_0 ) / Jxy
            ELSE
              localError="Jacobian Jxy is smaller than 1.0E-10_DP."
              CALL FLAG_ERROR(localError,ERR,ERROR,*999)
            END IF

            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE) THEN
              PERM_OVER_VIS_PARAM = PERM_OVER_VIS_PARAM_0
            ELSE
              MATERIAL_FACT = ( Jxy * POROSITY / POROSITY_0 )**2.0_DP
              PERM_OVER_VIS_PARAM = MATERIAL_FACT * PERM_OVER_VIS_PARAM_0
              !material modeling could use gradient information, or solve some PDE
            END IF

            IF(DIAGNOSTICS2) THEN
              IF(idebug1) THEN
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN = ", &
                  & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Jxy = ",Jxy,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"POROSITY = ",POROSITY,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"PERM_OVER_VIS_PARAM = ",PERM_OVER_VIS_PARAM,ERR,ERROR,*999)
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," ",ERR,ERROR,*999)
                idebug1 = .FALSE.
              ENDIF
            ENDIF

!!TODO: Think about symmetric problems. 
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

            !Loop over field components
            mhs=0          
            DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns ???
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN

                  !Loop over element columns
                  DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1

                      PGM=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                      PGN=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                      SUM = 0.0_DP
                      IF(mh==nh) THEN 
                        SUM = SUM + PGM * PGN
                      ENDIF
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) = &
                        & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) + SUM * RWG

                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                  PGM=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                  SUM = 0.0_DP
                  IF(mh==1) THEN 
                    SUM = SUM + PGM * POROSITY
                  ELSE IF(mh==2) THEN 
                    SUM = SUM + PGM * PERM_OVER_VIS_PARAM
                  END IF
                  RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs) = RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs) + SUM * RWG
                ENDIF

              ENDDO !ms
            ENDDO !mh
          ENDDO !ng


!-----------------------------------------------------------------------------------------------------------------------------------
! CHECK STIFFNESS MATRIX AND RHS VECTOR WITH CMHEART
          IF(DIAGNOSTICS5) THEN
            IF( ELEMENT_NUMBER == 1 ) THEN
              NDOFS = 0
              DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
                MESH_COMPONENT1 = fieldVariable%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                NDOFS = NDOFS + DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
              END DO

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element Matrix for element number 1 (Galerkin Projection):",ERR,ERROR,*999)
              DO mhs=1,NDOFS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"row number = ",mhs,ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS,NDOFS,NDOFS,&
                  & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,:), &
                  & '("",4(X,E13.6))','4(4(X,E13.6))',ERR,ERROR,*999)
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," ",ERR,ERROR,*999)
              END DO
            END IF
          END IF
!-----------------------------------------------------------------------------------------------------------------------------------
          
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            mhs=0          
            DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
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
       
        CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          fieldVariable=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=fieldVariable%VARIABLE_TYPE
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          SOURCE_BASIS=>SOURCE_FIELD%DECOMPOSITION%DOMAIN(SOURCE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
!             CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
            CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            TAU_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
            KAPPA_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
! WRITE(*,*)'TAU_PARAM ',TAU_PARAM          
            U_VALUE=0.0_DP
            IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, & 
                & EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              U_VALUE(1)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              U_VALUE(2)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
              IF(DEPENDENT_BASIS%NUMBER_OF_XI==3) THEN
                U_VALUE(3)=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
              ENDIF
            ENDIF
            !Calculate RWG.
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
              MESH_COMPONENT1=fieldVariable%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                 !Loop over element columns
                  DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
                  MESH_COMPONENT2=fieldVariable%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP &
                    & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                    DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      PGM=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                      PGN=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                      DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                        DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                          DXI_DX(mi,ni)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & DXI_DX(mi,ni)
                        END DO
                        PGMSI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        PGNSI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                      END DO !ni
                      SUM = 0.0_DP
                      !Calculate SUM 
                      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                        IF(mh==nh) THEN 
                          !This stiffness matrix contribution is without "integration" means ng=nd in fact = least square!
                          SUM = SUM + PGM * PGN
                        ENDIF
                        !
!                         IF(mh==nh) THEN 
!                           !This stiffness matrix happens with "integration" so the integral error is reduced!
!                           SUM = SUM + PGM * PGN * RWG
!                         ENDIF
!REDUCED SOBOLEV SMOOTHING
                          !This stiffness matrix contribution is with "integration" means ng=ng in fact!
                          SUM = SUM +    ( &
                            & TAU_PARAM*2.0_DP* ( &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S3,ng)) +&
                            & KAPPA_PARAM*2.0_DP* ( &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S1,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S2,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S3_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S2,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,ng))) !&
! no weighting either?
!                             & * RWG

                        EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) = &
                          & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) + SUM

                      ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                        IF(mh==nh.AND.mh<=NUMBER_OF_DIMENSIONS) THEN 
                          SUM = SUM + PGM * PGN
!REDUCED SOBOLEV SMOOTHING
                          !This stiffness matrix contribution is with "integration" means ng=ng in fact!
                        ENDIF
!REDUCED SOBOLEV SMOOTHING
                          !This stiffness matrix contribution is with "integration" means ng=ng in fact!
                          SUM = SUM +    ( &
                            & TAU_PARAM*2.0_DP* ( &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S3,ng)) +&
                            & KAPPA_PARAM*2.0_DP* ( &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S1,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S2,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S3_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S2,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,ng)+ &
                            & QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S3,ng)* &
                            & QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,ng))) !& 
! no weighting either?
!                             & * RWG

                          EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) = &
                            & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) + SUM
                        IF(nh==fieldVariable%NUMBER_OF_COMPONENTS.AND.mh<=NUMBER_OF_DIMENSIONS) THEN 
                          SUM=0.0_DP
                          !Calculate SUM 
                          DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            SUM=SUM+PGN*PGMSI(ni)*DXI_DX(ni,mh)
                          ENDDO !ni
                          EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) = &
                            & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) + SUM * RWG
                          EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(nhs,mhs) = &
                            & EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(nhs,mhs) + SUM * RWG
                        ENDIF
                      ENDIF 
                    ENDDO !ns
                  ENDDO !nh

                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
                IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                  IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                    SUM=0.0_DP
                    PGM=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                    SUM=U_VALUE(mh)*PGM
!                     SUM=42.0_DP*PGM
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+SUM
                  ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                    IF(mh<=NUMBER_OF_DIMENSIONS) THEN 
                      SUM=0.0_DP
                      PGM=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                      SUM=U_VALUE(mh)*PGM
!                       SUM=42.0_DP*PGM
                      SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+SUM
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            mhs=0          
            DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                  SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                ENDIF 
                IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                  RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                ENDIF 
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Galerkin projection type of a data fitting equations set class."
          CALL FLAG_ERROR(localError,ERR,ERROR,*999)
        END SELECT
        
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_FINITE_ELEMENT_CALCULATE")
#endif
    RETURN
999 CALL ERRORS("FITTING_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_FINITE_ELEMENT_CALCULATE")
#endif
    RETURN 1
  END SUBROUTINE FITTING_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the update-materials Galerkin projection.
  SUBROUTINE FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

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
!     TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP",ERR,ERROR,*999)
#endif

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE.OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)

        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing

        !-----------------------------------------------------------------
        ! d e p e n d e n t   f i e l d
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
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
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

              !component 1: dependent porosity variable, component 2: dependent permeability variable
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
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

              !component 1: dependent porosity variable, component 2: dependent permeability variable
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
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
              & " is invalid for an update-materials Galerkin projection"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! I N d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
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
                & (/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
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
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              !calculate number of components with one component for each dimension
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
              & " is invalid for an update-materials Galerkin projection"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=1
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS=2

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
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"Fitting Materials", &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & MATERIAL_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                DO I = 1, MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & I,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                END DO
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
              & " is invalid for an update-materials Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   s o u r c e   t y p e  
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
              & " is invalid for an update-materials Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
! ! !         CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
! ! !           SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
! ! !           CASE(EQUATIONS_SET_SETUP_START_ACTION)
! ! !             IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
! ! !               IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
! ! !                 GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
! ! !                 IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
! ! !                   CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
! ! !                   SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_1)
! ! !                     !Check that we are in 2D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=2) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_2)
! ! !                     !Check that we are in 2D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=2) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_2
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_1)
! ! !                     !Check that we are in 3D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=3) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_2)
! ! !                     !Check that we are in 3D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=3) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_2
! ! !                   CASE DEFAULT
! ! !                     LOCAL_ERROR="The specified analytic function type of "// &
! ! !                       & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                       & " is invalid for a moving mesh Galerkin projection."
! ! !                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                   END SELECT
! ! !                 ELSE
! ! !                   CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
! ! !                 ENDIF
! ! !              ELSE
! ! !                 CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
! ! !             IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
! ! !               ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
! ! !               IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
! ! !                 IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
! ! !                   CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
! ! !                 ENDIF
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
! ! !             IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
! ! !                 IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
! ! !                   CALL FITTING_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
! ! !                 ELSE
! ! !                   CALL FLAG_ERROR("Equations set analtyic has not been finished.",ERR,ERROR,*999)
! ! !                 ENDIF
! ! !               ELSE
! ! !                 CALL FLAG_ERROR("Equations set analtyic is not associated.",ERR,ERROR,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set dependent has not been finished.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE DEFAULT
! ! !             LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
! ! !               & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
! ! !               & " is invalid for an update-materials Galerkin projection."
! ! !             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !           END SELECT

        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e   
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
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
              & " is invalid for an update-materials Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   c a s e   d e f a u l t
        !-----------------------------------------------------------------
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for an update-materials Galerkin projection."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal an update-materials Galerkin projection subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Galerkin projection type of a data fitting equations set class.
  SUBROUTINE FITTING_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Galerkin projection on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATIONS_SET_SETUP",ERR,ERROR,*999)
#endif

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_DATA_FITTING_SUBTYPE)
        CALL FITTING_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE)
        CALL FITTING_EQUATIONS_SET_VECTORDATA_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
        CALL FITTING_EQUATIONS_SET_VECTORDATA_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EquationsSet_DataPointVectorStaticFittingSubtype, &
        &  EquationsSet_DataPointVectorQuasistaticFittingSubtype)
        CALL FITTING_EQUATIONS_SET_VECTORDATA_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
        & EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE)
        CALL FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Galerkin projection type of a data fitting equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_SETUP")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_SETUP",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_SETUP")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Galerkin projection type of an data fitting equations set class.
  SUBROUTINE FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
#endif
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_DATA_FITTING_SUBTYPE)        
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
      CASE(EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
        & EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE)        
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
      CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)        
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
      CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
        & EquationsSet_DataPointVectorStaticFittingSubtype,EquationsSet_DataPointVectorQuasistaticFittingSubtype)
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
          & " is not valid for a Galerkin projection type of an data fitting equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a data fitting equation set class.
  SUBROUTINE FITTING_EQUATIONS_SET_CLASS_TYPE_GET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_TYPE !<On return, the equation type
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SUBTYPE !<On return, the equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR,*999)
#endif

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_FITTING_CLASS) THEN
        EQUATIONS_TYPE=EQUATIONS_SET%TYPE
        EQUATIONS_SUBTYPE=EQUATIONS_SET%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Equations set is not the data fitting type",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_CLASS_TYPE_GET")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_CLASS_TYPE_GET")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a data fitting equation set class.
  SUBROUTINE FITTING_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_TYPE !<The equation type
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SUBTYPE !<The equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATIONS_SET_CLASS_SET",ERR,ERROR,*999)
#endif

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_TYPE)
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        CALL FITTING_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SUBTYPE,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a data fitting equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_CLASS_TYPE_SET")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_CLASS_TYPE_SET",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_CLASS_TYPE_SET")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a Galerkin projection type of a data fitting equations set class.
  SUBROUTINE FITTING_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
#endif
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_DATA_FITTING_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_STANDARD_DATA_FITTING_SUBTYPE
      CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE
      CASE(EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE
      CASE(EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE
      CASE(EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE
      CASE(EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE
      CASE(EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE
      CASE(EquationsSet_DataPointVectorStaticFittingSubtype)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EquationsSet_DataPointVectorStaticFittingSubtype
      CASE(EquationsSet_DataPointVectorQuasistaticFittingSubtype)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FITTING_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EquationsSet_DataPointVectorQuasistaticFittingSubtype
      CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Galerkin projection type of a data fitting equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_SUBTYPE_SET")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_SUBTYPE_SET")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the standard Galerkin projection.
  SUBROUTINE FITTING_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE !,NUMBER_OF_DIMENSIONS
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
!     TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATION_SET_STANDARD_SETUP",ERR,ERROR,*999)
#endif

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_DATA_FITTING_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)

        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing

        !-----------------------------------------------------------------
        ! d e p e n d e n t   f i e l d
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
              & " is invalid for a standard Galerkin projection"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   s o u r c e   t y p e 
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
              & " is invalid for a standard Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
! ! !         CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
! ! !           SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
! ! !           CASE(EQUATIONS_SET_SETUP_START_ACTION)
! ! !             IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
! ! !               IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
! ! !                 GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
! ! !                 IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
! ! !                   CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
! ! !                   SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_1)
! ! !                     !Check that we are in 2D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=2) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_2)
! ! !                     !Check that we are in 2D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=2) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_2
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_1)
! ! !                     !Check that we are in 3D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=3) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_2)
! ! !                     !Check that we are in 3D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=3) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_2
! ! !                   CASE DEFAULT
! ! !                     LOCAL_ERROR="The specified analytic function type of "// &
! ! !                       & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                       & " is invalid for a standard Galerkin projection."
! ! !                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                   END SELECT
! ! !                 ELSE
! ! !                   CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
! ! !                 ENDIF
! ! !              ELSE
! ! !                 CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
! ! !             IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
! ! !               ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
! ! !               IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
! ! !                 IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
! ! !                   CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
! ! !                 ENDIF
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
! ! !             IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
! ! !                 IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
! ! !                   CALL FITTING_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
! ! !                 ELSE
! ! !                   CALL FLAG_ERROR("Equations set analtyic has not been finished.",ERR,ERROR,*999)
! ! !                 ENDIF
! ! !               ELSE
! ! !                 CALL FLAG_ERROR("Equations set analtyic is not associated.",ERR,ERROR,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set dependent has not been finished.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE DEFAULT
! ! !             LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
! ! !               & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
! ! !               & " is invalid for a standard Galerkin projection."
! ! !             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !           END SELECT

        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e   
        !-----------------------------------------------------------------
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
              & " is invalid for a standard Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   c a s e   d e f a u l t
        !-----------------------------------------------------------------
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Galerkin projection."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard Galerkin projection subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_STANDARD_SETUP")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_STANDARD_SETUP")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_STANDARD_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the vector data Galerkin projection.
  SUBROUTINE FITTING_EQUATIONS_SET_VECTORDATA_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,GEOMETRIC_COMPONENT_NUMBER
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,I !,MATERIAL_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: dependentFieldNumberOfVariables
    INTEGER(INTG) :: dimensionIdx
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
!     TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_EQUATION_SET_VECTORDATA_SETUP",ERR,ERROR,*999)
#endif

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE, &
        & EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing

        !-----------------------------------------------------------------
        ! d e p e n d e n t   f i e l d
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
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
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
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
! !                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
! !                   & 1,ERR,ERROR,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
! !                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
! !                   & 1,ERR,ERROR,*999)

!                 DO I=1,1
                DO I=1,NUMBER_OF_DIMENSIONS
                  !Default to the geometric interpolation setup
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                END DO
              ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & NUMBER_OF_DIMENSIONS+1,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & NUMBER_OF_DIMENSIONS+1,ERR,ERROR,*999)
                DO I=1,NUMBER_OF_DIMENSIONS+1
                  !Default to the geometric interpolation setup
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                END DO
              ENDIF
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE) THEN
!                   DO I=1,NUMBER_OF_DIMENSIONS
                  DO I=1,1
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE) THEN
                  DO I=1,NUMBER_OF_DIMENSIONS+1
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                ENDIF
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
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE) THEN
                  DO I=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE) THEN
                  DO I=1,NUMBER_OF_DIMENSIONS+1
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                ENDIF
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
                !Other solutions not defined yet
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
              & " is invalid for an update-materials Galerkin projection"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! I N d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
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
                & (/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
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
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              !calculate number of components with one component for each dimension
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
              & " is invalid for an update-materials Galerkin projection"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
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
                  & 1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
                  & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,0.0_DP,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   s o u r c e   t y p e  
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
          CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE, &
            & EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
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
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                  & 1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                  & (/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
              IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                !These 2 parameter sets will contain the fitted hermite/lagrange velocity field
!                 CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)
!                 CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)

!                 CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_INPUT_DATA3_SET_TYPE,ERR,ERROR,*999)
!                 CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_BOUNDARY_SET_TYPE,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard PEE problem"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for a PPE equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
! ! !         CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
! ! !           SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
! ! !           CASE(EQUATIONS_SET_SETUP_START_ACTION)
! ! !             IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
! ! !               IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
! ! !                 GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
! ! !                 IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
! ! !                   CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
! ! !                   SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_1)
! ! !                     !Check that we are in 2D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=2) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_2)
! ! !                     !Check that we are in 2D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=2) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_2
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_1)
! ! !                     !Check that we are in 3D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=3) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_2)
! ! !                     !Check that we are in 3D
! ! !                     IF(NUMBER_OF_DIMENSIONS/=3) THEN
! ! !                       LOCAL_ERROR="The number of geometric dimensions of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_2
! ! !                   CASE DEFAULT
! ! !                     LOCAL_ERROR="The specified analytic function type of "// &
! ! !                       & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
! ! !                       & " is invalid for a standard Galerkin projection."
! ! !                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !                   END SELECT
! ! !                 ELSE
! ! !                   CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
! ! !                 ENDIF
! ! !              ELSE
! ! !                 CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
! ! !             IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
! ! !               ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
! ! !               IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
! ! !                 IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
! ! !                   CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
! ! !                 ENDIF
! ! !               ENDIFstandard
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
! ! !             IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
! ! !                 IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
! ! !                   CALL FITTING_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
! ! !                 ELSE
! ! !                   CALL FLAG_ERROR("Equations set analtyic has not been finished.",ERR,ERROR,*999)
! ! !                 ENDIF
! ! !               ELSE
! ! !                 CALL FLAG_ERROR("Equations set analtyic is not associated.",ERR,ERROR,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FLAG_ERROR("Equations set dependent has not been finished.",ERR,ERROR,*999)
! ! !             ENDIF
! ! !           CASE DEFAULT
! ! !             LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
! ! !               & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
! ! !               & " is invalid for a standard Galerkin projection."
! ! !             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
! ! !           END SELECT

        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e   
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
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
              & " is invalid for a vector data Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a vector data Galerkin projection."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      CASE(EquationsSet_DataPointVectorStaticFittingSubtype, &
        &  EquationsSet_DataPointVectorQuasistaticFittingSubtype)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL FITTING_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          ! Do nothing
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        ! (this field will hold the mesh fitted data from the data points field)
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              !start field creation with name 'DEPENDENT_FIELD'
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              !start creation of a new field
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                !define new created field to be dependent
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              !look for decomposition rule already defined
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              !apply decomposition rule found on new created field
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              !point new field to geometric field
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET% & 
                & GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
              !set number of variables to 2 (U, delUdelN)
              dependentFieldNumberOfVariables=2
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & dependentFieldNumberOfVariables,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                & [FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
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
              !calculate number of components with one component for each dimension
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              DO I=1,NUMBER_OF_DIMENSIONS
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              END DO
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, & 
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
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
              & " is invalid for an update-materials Galerkin projection"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
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
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                ! Sobelov smoothing material parameters- tau and kappa
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
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
                  & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,0.0_DP,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   t y p e  
        ! (this field holds the data point based field of vectors to map to the dependent field)
        !-----------------------------------------------------------------
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
              !label the field
              CALL FIELD_LABEL_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",ERR,ERROR, & 
                & *999)
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
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                & 2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
              ! U Variable: data point vectors
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              ! V Variable: data point weights
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                & FIELD_V_VARIABLE_TYPE,1,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_V_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO dimensionIdx = 1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,dimensionIdx,FIELD_DATA_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,1,FIELD_DATA_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
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
              ! U (vector) variable
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              !calculate number of components with one component for each dimension
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              ! V (weight) variable
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_DATA_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_DATA_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
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
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard PEE problem"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e   
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              IF (EQUATIONS_SET%SUBTYPE==EquationsSet_DataPointVectorStaticFittingSubtype) THEN
                CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
              ELSE IF (EQUATIONS_SET%SUBTYPE==EquationsSet_DataPointVectorQuasistaticFittingSubtype) THEN
                CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_QUASISTATIC,ERR,ERROR,*999)                
              ENDIF
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
              & " is invalid for a vector data Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a vector data Galerkin projection."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a vector data Galerkin projection subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_VECTORDATA_SETUP")
#endif
    RETURN
999 CALL ERRORS("FITTING_EQUATIONS_SET_VECTORDATA_SETUP",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_EQUATIONS_SET_VECTORDATA_SETUP")
#endif
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_VECTORDATA_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the Galerkin projection problem.
  SUBROUTINE FITTING_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Galerkin projection on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_PROBLEM_SETUP",ERR,ERROR,*999)
#endif

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(Problem_DataPointVectorStaticFittingSubtype)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(Problem_DataPointVectorQuasistaticFittingSubtype)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Galerkin projection type of a data fitting problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_SETUP")
#endif
    RETURN
999 CALL ERRORS("FITTING_PROBLEM_SETUP",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_SETUP")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a Galerkin projection type .
  SUBROUTINE FITTING_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
#endif
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_FITTING_CLASS
        PROBLEM%TYPE=PROBLEM_DATA_FITTING_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STANDARD_DATA_FITTING_SUBTYPE     
      CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_FITTING_CLASS
        PROBLEM%TYPE=PROBLEM_DATA_FITTING_TYPE
        PROBLEM%SUBTYPE=PROBLEM_VECTOR_DATA_FITTING_SUBTYPE     
      CASE(PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_FITTING_CLASS
        PROBLEM%TYPE=PROBLEM_DATA_FITTING_TYPE
        PROBLEM%SUBTYPE=PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE     
      CASE(Problem_DataPointVectorStaticFittingSubtype)
        PROBLEM%CLASS=PROBLEM_FITTING_CLASS
        PROBLEM%TYPE=PROBLEM_DATA_FITTING_TYPE
        PROBLEM%SUBTYPE=Problem_DataPointVectorStaticFittingSubtype     
      CASE(Problem_DataPointVectorQuasistaticFittingSubtype)
        PROBLEM%CLASS=PROBLEM_FITTING_CLASS
        PROBLEM%TYPE=PROBLEM_DATA_FITTING_TYPE
        PROBLEM%SUBTYPE=Problem_DataPointVectorQuasistaticFittingSubtype     
      CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Galerkin projection type of a data fitting problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_SUBTYPE_SET")
#endif
    RETURN
999 CALL ERRORS("FITTING_PROBLEM_SUBTYPE_SET",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_SUBTYPE_SET")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the standard Galerkin projections problem.
  SUBROUTINE FITTING_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
#if DEBUG
    CALL ENTERS("FITTING_PROBLEM_STANDARD_SETUP",ERR,ERROR,*999)
#endif

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_STANDARD_DATA_FITTING_SUBTYPE) THEN
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
              & " is invalid for a standard Galerkin projection."
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
              & " is invalid for a standard Galerkin projection."
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
                & " is invalid for a standard Galerkin projection."
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
              & " is invalid for a standard Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Galerkin projection."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard Galerkin projection subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_STANDARD_SETUP")
#endif
    RETURN
999 CALL ERRORS("FITTING_PROBLEM_STANDARD_SETUP",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_STANDARD_SETUP")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_STANDARD_SETUP

 !
  !================================================================================================================================
  !

  !>Sets up the vector data Galerkin projections problem.
  SUBROUTINE FITTING_PROBLEM_VECTORDATA_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

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
    
#if DEBUG
    CALL ENTERS("FITTING_PROBLEM_VECTORDATA_SETUP",ERR,ERROR,*999)
#endif

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_VECTOR_DATA_FITTING_SUBTYPE.OR. &
        & PROBLEM%SUBTYPE==Problem_DataPointVectorStaticFittingSubtype .OR. &
        & PROBLEM%SUBTYPE==Problem_DataPointVectorQuasistaticFittingSubtype .OR. &
        & PROBLEM%SUBTYPE==PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE) THEN
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
              & " is invalid for a vector data Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            IF(PROBLEM%SUBTYPE==Problem_DataPointVectorStaticFittingSubtype) THEN
              CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
            ELSE
              CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            ENDIF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a vector data Galerkin projection."
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
                & " is invalid for a vector data Galerkin projection."
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
            IF(PROBLEM%SUBTYPE==Problem_DataPointVectorStaticFittingSubtype) THEN
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            ELSE
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
            ENDIF
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
              & " is invalid for a vector data Galerkin projection."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a vector data Galerkin projection."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a vector data Galerkin projection subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_VECTORDATA_SETUP")
#endif
    RETURN
999 CALL ERRORS("FITTING_PROBLEM_VECTORDATA_SETUP",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_VECTORDATA_SETUP")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_VECTORDATA_SETUP

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a data fitting problem class.
  SUBROUTINE FITTING_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<On return, the problem type
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<On return, the proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
#if DEBUG
    CALL ENTERS("FITTING_PROBLEM_CLASS_TYPE_GET",ERR,ERROR,*999)
#endif

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%CLASS==PROBLEM_FITTING_CLASS) THEN
        PROBLEM_EQUATION_TYPE=PROBLEM%TYPE
        PROBLEM_SUBTYPE=PROBLEM%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Problem is not data fitting class",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_CLASS_TYPE_GET")
#endif
    RETURN
999 CALL ERRORS("FITTING_PROBLEM_CLASS_TYPE_GET",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_CLASS_TYPE_GET")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a data fitting problem class.
  SUBROUTINE FITTING_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE !<The problem type
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_PROBLEM_CLASS_SET",ERR,ERROR,*999)
#endif
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_EQUATION_TYPE)
       CASE(PROBLEM_DATA_FITTING_TYPE)
        CALL FITTING_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem equation type "//TRIM(NUMBER_TO_VSTRING(PROBLEM_EQUATION_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a data fitting problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_CLASS_TYPE_SET")
#endif
    RETURN
999 CALL ERRORS("FITTING_PROBLEM_CLASS_TYPE_SET",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PROBLEM_CLASS_TYPE_SET")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !   

  !>Evaluates the deformation gradient tensor at a given Gauss point
  SUBROUTINE FITTING_GAUSS_DEFORMATION_GRADIENT_TENSOR(REFERENCE_GEOMETRIC_INTERPOLATED_POINT, &
    & GEOMETRIC_INTERPOLATED_POINT, DXDY, Jxy, ERR, ERROR, *)    

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: REFERENCE_GEOMETRIC_INTERPOLATED_POINT, GEOMETRIC_INTERPOLATED_POINT
    REAL(DP) :: DXDY(3,3)  !DXDY - Deformation Gradient Tensor  
    REAL(DP) :: Jxy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string   
    !Local Variables
    INTEGER(INTG) :: derivative_idx,component_idx,xi_idx 
    REAL(DP) :: DXDXI(3,3),DYDXI(3,3),DXIDY(3,3)
    REAL(DP) :: Jyxi

#if DEBUG
    CALL ENTERS("FITTING_GAUSS_DEFORMATION_GRADIENT_TENSOR",ERR,ERROR,*999)
#endif

    !--- ToDo: Needs to be generalized such that it also works for 2D
    DO component_idx=1,3 !Always 3 components - 3D
      DO xi_idx=1,3 !Thus 3 element coordinates
        derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7      
        DXDXI(component_idx,xi_idx)=GEOMETRIC_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dx/dxi
        DYDXI(component_idx,xi_idx)=REFERENCE_GEOMETRIC_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dy/dxi (y = referential)
      ENDDO
    ENDDO


    CALL INVERT(DYDXI,DXIDY,Jyxi,ERR,ERROR,*999) !dy/dxi -> dxi/dy 

    CALL MATRIX_PRODUCT(DXDXI,DXIDY,DXDY,ERR,ERROR,*999) !dx/dxi * dxi/dy = dx/dy (deformation gradient tensor, F)

    Jxy=DETERMINANT(DXDY,ERR,ERROR)


#if DEBUG
    CALL EXITS("FITTING_GAUSS_DEFORMATION_GRADIENT_TENSOR")
#endif
    RETURN
999 CALL ERRORS("FITTING_GAUSS_DEFORMATION_GRADIENT_TENSOR",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_GAUSS_DEFORMATION_GRADIENT_TENSOR")
#endif
    RETURN 1
  END SUBROUTINE FITTING_GAUSS_DEFORMATION_GRADIENT_TENSOR

  !
  !================================================================================================================================
  !


  !>Sets up the output type for a data fitting problem class.
  SUBROUTINE FITTING_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#if DEBUG
    CALL ENTERS("FITTING_PRE_SOLVE",ERR,ERROR,*999)
#endif

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE)
!               do nothing
            CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
!               do nothing
            CASE(PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE)
!               do nothing
            CASE(Problem_DataPointVectorStaticFittingSubtype)
!               do nothing
            CASE(Problem_DataPointVectorQuasistaticFittingSubtype)
!               do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
! !               IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1)THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
                !Update indpendent data fields
                CALL FITTING_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
! !                 CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"While loop... ",ERR,ERROR,*999)
! !               ELSE
! !                 CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"While loop... ",ERR,ERROR,*999)
! !               ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a data fitting problem class."
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
       
#if DEBUG
    CALL EXITS("FITTING_PRE_SOLVE")
#endif
    RETURN
999 CALL ERRORS("FITTING_PRE_SOLVE",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PRE_SOLVE")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a data fitting problem class.
  SUBROUTINE FITTING_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

#if DEBUG
    CALL ENTERS("FITTING_POST_SOLVE",ERR,ERROR,*999)
#endif
    NULLIFY(SOLVER2)
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE,PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE, &
              & PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE)
!               do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE, &
              & Problem_DataPointVectorStaticFittingSubtype)
              CALL FITTING_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(Problem_DataPointVectorQuasistaticFittingSubtype)
              ! do nothing
            CASE(PROBLEM_VECTOR_DATA_PRE_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
!               do nothing
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a fitting type of a classical field problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF   
#if DEBUG
    CALL EXITS("FITTING_POST_SOLVE")
#endif
    RETURN
999 CALL ERRORS("FITTING_POST_SOLVE",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_POST_SOLVE")
#endif
    RETURN 1
  END SUBROUTINE FITTING_POST_SOLVE


  !
  !================================================================================================================================
  !


  !>Output data post solve
  SUBROUTINE FITTING_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

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
    CHARACTER(7) :: FILE
    CHARACTER(7) :: OUTPUT_FILE

#if DEBUG
    CALL ENTERS("FITTING_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)
#endif

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!       write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE,PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE, &
              & PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
              & Problem_DataPointVectorStaticFittingSubtype)
!               do nothing
            CASE(PROBLEM_VECTOR_DATA_PRE_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
!               do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
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
                         WRITE(OUTPUT_FILE,'("DATA_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("DATA_",I0)') CURRENT_LOOP_ITERATION
                        END IF
                        FILE=OUTPUT_FILE
!          FILE="TRANSIENT_OUTPUT"
                        METHOD="FORTRAN"
                        EXPORT_FIELD=.TRUE.
                        IF(EXPORT_FIELD) THEN          
                          IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",ERR,ERROR,*999)
                            CALL FLUID_MECHANICS_IO_WRITE_FITTED_FIELD(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER, &
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
                & " is not valid for a fitting equation of a classical field problem class."
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
#if DEBUG
    CALL EXITS("FITTING_POST_SOLVE_OUTPUT_DATA")
#endif
    RETURN
999 CALL ERRORS("FITTING_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_POST_SOLVE_OUTPUT_DATA")
#endif
    RETURN 1
  END SUBROUTINE FITTING_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Update input data conditions for field fitting
  SUBROUTINE FITTING_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

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

!     REAL(DP) :: CURRENT_TIME,TIME_INCREMENT

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,CURRENT_LOOP_ITERATION
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_VEL_NEW_DATA(:)!,INPUT_VEL_OLD_DATA(:)
!     REAL(DP), POINTER :: INPUT_VEL_LABEL_DATA(:) !,INPUT_VEL_U_DATA(:),INPUT_VEL_V_DATA(:),INPUT_VEL_W_DATA(:)
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control loop to solve.
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    LOGICAL :: BOUNDARY_UPDATE

    BOUNDARY_UPDATE=.FALSE.

#if DEBUG
    CALL ENTERS("FITTING_PRE_SOLVE_UPDATE_INPUT_DATA",ERR,ERROR,*999)
#endif

    NULLIFY(INPUT_VEL_NEW_DATA)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE)
!               do nothing
            CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
!               do nothing
            CASE(PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE)
!               do nothing
            CASE(Problem_DataPointVectorStaticFittingSubtype)
!               do nothing
            CASE(Problem_DataPointVectorQuasistaticFittingSubtype)
!               do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
!               do nothing
                CONTROL_TIME_LOOP=>CONTROL_LOOP
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
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
                      CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                      !this is the current time step
!\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
                      INPUT_TYPE=1
                      INPUT_OPTION=1
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
! ! !                         & FIELD_INPUT_DATA1_SET_TYPE,INPUT_VEL_NEW_DATA,ERR,ERROR,*999)
                        & FIELD_VALUES_SET_TYPE,INPUT_VEL_NEW_DATA,ERR,ERROR,*999)
! ! !                       CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_NEW_DATA, & 
! ! !                         & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_NEW_DATA, & 
                        & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CURRENT_LOOP_ITERATION,1.0_DP)
                      !this is the previous time step
! ! ! !\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
! ! !                       INPUT_TYPE=1
! ! !                       INPUT_OPTION=2
! ! !                       CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
! ! !                         & FIELD_INPUT_DATA2_SET_TYPE,INPUT_VEL_OLD_DATA,ERR,ERROR,*999)
! ! !                       CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_OLD_DATA, & 
! ! !                         & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
! ! !                       !this is the interior flag
! ! !                       INPUT_TYPE=1
! ! !                       INPUT_OPTION=3
! ! !                       CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
! ! !                         & FIELD_INPUT_LABEL_SET_TYPE,INPUT_VEL_LABEL_DATA,ERR,ERROR,*999)
! ! !                       CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_LABEL_DATA, & 
! ! !                         & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
! ! !                       !this is the reference U velocity
! ! !                       INPUT_TYPE=1
! ! !                       INPUT_OPTION=4
! ! !                       CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
! ! !                         & FIELD_INPUT_VEL1_SET_TYPE,INPUT_VEL_U_DATA,ERR,ERROR,*999)
! ! !                       CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_U_DATA, & 
! ! !                         & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
! ! !                       !this is the reference V velocity
! ! !                       INPUT_TYPE=1
! ! !                       INPUT_OPTION=5
! ! !                       CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
! ! !                         & FIELD_INPUT_VEL2_SET_TYPE,INPUT_VEL_V_DATA,ERR,ERROR,*999)
! ! !                       CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_V_DATA, & 
! ! !                         & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
! ! !                       !this is the reference W velocity
! ! !                       INPUT_TYPE=1
! ! !                       INPUT_OPTION=6
! ! !                       CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
! ! !                         & FIELD_INPUT_VEL3_SET_TYPE,INPUT_VEL_W_DATA,ERR,ERROR,*999)
! ! !                       CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_VEL_W_DATA, & 
! ! !                         & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
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
                & " is not valid for a vector data type of a fitting field problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
             & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, & 
             & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
#if DEBUG
    CALL EXITS("FITTING_PRE_SOLVE_UPDATE_INPUT_DATA")
#endif
    RETURN
999 CALL ERRORS("FITTING_PRE_SOLVE_UPDATE_INPUT_DATA",ERR,ERROR)
#if DEBUG
    CALL EXITS("FITTING_PRE_SOLVE_UPDATE_INPUT_DATA")
#endif
    RETURN 1
  END SUBROUTINE FITTING_PRE_SOLVE_UPDATE_INPUT_DATA

  !
  !================================================================================================================================
  !

 
END MODULE FITTING_ROUTINES
