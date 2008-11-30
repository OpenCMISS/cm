!> \file
!> $Id: boundary_condition_rountines.f90 204 2008-10-31 05:02:38Z tingy $
!> \author Ting Yu
!> \brief This module set the boudary conditions for the given equation set
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

!>This module handles all analytic analysis routines.
MODULE BOUNDARY_CONDITION_ROUTINES

  USE BASE_ROUTINES
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE NODE_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE, BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE
    
CONTAINS  

  !
  !================================================================================================================================
  !  
  
  !>Sets a fixed condition for the equation set on the specified node.
  SUBROUTINE BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE(EQUATIONS_SET,SET_TYPE,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER, &
    & VARIABLE_NUMBER,CONDITION,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the fixed condition for
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The field parameter set type
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The node_number to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<The variable number to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The fixed condition to set
    REAL(DP), INTENT(IN) :: VALUE !<The value of the fixed condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny, global_ny,GLOBAL_NODE_NUMBER, USER_NODE_NUMBER
    LOGICAL :: NODE_EXISTS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
          CALL FLAG_ERROR("Fixed conditions have been finished for this equations set.",ERR,ERROR,*999)
        ELSE
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              USER_NODE_NUMBER=DEPENDENT_FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES% &
                & NODES(NODE_NUMBER)%USER_NUMBER
              CALL NODE_CHECK_EXISTS(USER_NODE_NUMBER,DEPENDENT_FIELD%REGION,NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
              IF(.NOT.NODE_EXISTS) THEN
                CALL FLAG_ERROR("Invalid node number.",ERR,ERROR,*999)
              ELSE
                local_ny=DEPENDENT_FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                  & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,NODE_NUMBER,1)
                global_ny=DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                SELECT CASE(CONDITION)
                CASE(EQUATIONS_SET_NOT_FIXED)
                  FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED
                CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                  FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,local_ny, &
                    & EQUATIONS_SET_FIXED_BOUNDARY_CONDITION,ERR,ERROR,*999)
                  !Should not update start and finish here as these are collective operations.
                  !CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(DEPENDENT_FIELD,SET_TYPE,DERIVATIVE_NUMBER, &
                    & NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER,VALUE,ERR,ERROR,*999)
                  !CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,SET_TYPE,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The specified condition of "//TRIM(NUMBER_TO_VSTRING(CONDITION,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
            ELSE
              CALL FLAG_ERROR("The dependent field is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem dependent has not been finished",ERR,ERROR,*999)              
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set fixed conditions are not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE
  
  !
  !================================================================================================================================
  !
  
  !>Update the boundary condition from analytic value.
  SUBROUTINE BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:)
    TYPE(FIELD_TYPE), POINTER :: FIELD 
    INTEGER(INTG) :: var_idx,comp_idx,node_idx,node_number,dev_idx
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    REAL(DP) :: VALUE
    
    
    CALL ENTERS("BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE",ERR,ERROR,*999)
    
    FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET(FIELD,FIELD_ANALYTIC_SET_TYPE,ANALYTIC_PARAMETERS,ERR,ERROR,*999)
      DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
        DO comp_idx=1,FIELD%VARIABLES(var_idx)%NUMBER_OF_COMPONENTS
          DOMAIN_NODES=>FIELD%VARIABLES(1)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
          IF(ASSOCIATED(DOMAIN_NODES)) THEN
            NODES_MAPPING=>FIELD%VARIABLES(1)%COMPONENTS(comp_idx)%DOMAIN%MAPPINGS%NODES
            DO node_idx=1,NODES_MAPPING%NUMBER_OF_INTERNAL
              node_number=NODES_MAPPING%INTERNAL_LIST(node_idx)
              DO dev_idx=1,DOMAIN_NODES%NODES(node_number)%NUMBER_OF_DERIVATIVES
                ! Set the boundary condition for dependent field \todo for rectangular system only
                IF(var_idx==1.AND.DOMAIN_NODES%NODES(node_number)%NUMBER_OF_SURROUNDING_ELEMENTS &
                  & < 2**FIELD%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN   
                  VALUE = ANALYTIC_PARAMETERS(FIELD%VARIABLES(var_idx)%COMPONENTS(comp_idx)%PARAM_TO_DOF_MAP% &
                    & NODE_PARAM2DOF_MAP(dev_idx,node_number,var_idx))
                  CALL BOUNDARY_CONDITION_FIXED_CONDITIONS_SET_NODE(EQUATIONS_SET,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,dev_idx, & 
                    & node_number,comp_idx,var_idx,EQUATIONS_SET_FIXED_BOUNDARY_CONDITION,VALUE,ERR,ERROR,*999)
                ENDIF       
              ENDDO ! dev_idx
            ENDDO ! node_idx
          ELSE
            CALL FLAG_ERROR("Domain nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ENDDO ! comp_idx
      ENDDO ! var_idx
    ELSE
      CALL FLAG_ERROR("The field is not associated.",ERR,ERROR,*999) 
    ENDIF
    
    CALL EXITS("BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE


  !
  !================================================================================================================================
  !  
  
 
END MODULE BOUNDARY_CONDITION_ROUTINES
