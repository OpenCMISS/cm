!> \file
!> $Id: boundary_conditions_rountines.f90 204 2008-10-31 05:02:38Z tingy $
!> \author Ting Yu
!> \brief This module set the boundary conditions for the given equation set
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

!>This module handles all boundary conditions routines.
MODULE BOUNDARY_CONDITIONS_ROUTINES

  USE BASE_ROUTINES
  USE CMISS_MPI
  USE CONSTANTS
  USE COMP_ENVIRONMENT
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MPI
  USE NODE_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions BOUNDARY_CONDITIONS_ROUTINES::BoundaryConditions
  !> \brief Boundary conditions type parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NOT_FIXED=0 !<The dof is not fixed. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED=1 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MIXED=2 !<The dof is set as a mixed boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE BOUNDARY_CONDITIONS_SET_LOCAL_DOF
    MODULE PROCEDURE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1
    MODULE PROCEDURE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS
  END INTERFACE !BOUNDARY_CONDITIONS_SET_LOCAL_DOF

  PUBLIC BOUNDARY_CONDITION_NOT_FIXED,BOUNDARY_CONDITION_FIXED,BOUNDARY_CONDITION_MIXED

  PUBLIC BOUNDARY_CONDITIONS_CREATE_FINISH,BOUNDARY_CONDITIONS_CREATE_START,BOUNDARY_CONDITIONS_DESTROY
  
  PUBLIC BOUNDARY_CONDITIONS_SET_CONSTANT,BOUNDARY_CONDITIONS_SET_LOCAL_DOF,BOUNDARY_CONDITIONS_SET_ELEMENT, &
    & BOUNDARY_CONDITIONS_SET_NODE
  
  PUBLIC EQUATIONS_SET_BOUNDARY_CONDITIONS_GET
  
  PUBLIC BOUNDARY_CONDITION_PARAM_SET_UPDATE_FROM_ANAL_VALUE
    
CONTAINS  

  !
  !================================================================================================================================
  !

  !>Finish the creation of boundary conditions.
  SUBROUTINE BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditiosn to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: computational_node_idx,DOF_START,MPI_IERROR,SEND_COUNT,variable_type_idx
    INTEGER(INTG), ALLOCATABLE :: DISPLACEMENTS(:),RECEIVE_COUNTS(:)    
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITION_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: VARIABLE_DOMAIN_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have already been finished.",ERR,ERROR,*999)        
      ELSE
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP)) THEN
          IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>1) THEN
            !Transfer all the boundary conditions to all the computational nodes.
            ALLOCATE(RECEIVE_COUNTS(0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate recieve counts.",ERR,ERROR,*999)
            ALLOCATE(DISPLACEMENTS(0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate displacements.",ERR,ERROR,*999)
!!TODO \todo Look at this.
            DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
              BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type_idx)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
                FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(VARIABLE_DOMAIN_MAPPING)) THEN
                    SEND_COUNT=0
                    DOF_START=0
                    DO computational_node_idx=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
                      DISPLACEMENTS(computational_node_idx)=DOF_START
                      RECEIVE_COUNTS(computational_node_idx)=VARIABLE_DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(computational_node_idx)
                      IF(RECEIVE_COUNTS(computational_node_idx)>SEND_COUNT) SEND_COUNT=RECEIVE_COUNTS(computational_node_idx)
                      DOF_START=DOF_START+VARIABLE_DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(computational_node_idx)
                    ENDDO !computational_node_idx
                    CALL MPI_ALLGATHERV(MPI_IN_PLACE,SEND_COUNT,MPI_INTEGER,BOUNDARY_CONDITIONS% &
                      & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type_idx)%PTR%GLOBAL_BOUNDARY_CONDITIONS, &
                      & RECEIVE_COUNTS,DISPLACEMENTS,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="Field variable domain mapping is not associated for variable type "// &
                      & TRIM(NUMBER_TO_VSTRING(variable_type_idx,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Field variable is not associated for variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(variable_type_idx,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Boundary condition variable is not associated for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(variable_type_idx,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO ! variable_type_idx
            IF(ALLOCATED(RECEIVE_COUNTS)) DEALLOCATE(RECEIVE_COUNTS)
            IF(ALLOCATED(DISPLACEMENTS)) DEALLOCATE(DISPLACEMENTS)
          ENDIF
          !Set the finished flag
          BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.TRUE.
        ELSE
          CALL FLAG_ERROR("Boundary conditions variable type map is not allocated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Boundary conditions:",ERR,ERROR,*999)
      DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type_idx)%PTR
        FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
        VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Variable type = ",variable_type_idx,ERR,ERROR,*999)
        IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global dofs = ",VARIABLE_DOMAIN_MAPPING% &
            & NUMBER_OF_GLOBAL,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,8,8, &
            & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS,'("    Global BCs:",8(X,I8))','(15X,8(X,I8))', &
            & ERR,ERROR,*999)
        ELSE
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Not mapped",ERR,ERROR,*999)
        ENDIF
      ENDDO !variable_type_idx
    ENDIF
    
    CALL EXITS("BOUNDARY_CONDITIONS_CREATE_FINISH")
    RETURN
999 IF(ALLOCATED(RECEIVE_COUNTS)) DEALLOCATE(RECEIVE_COUNTS)
    IF(ALLOCATED(DISPLACEMENTS)) DEALLOCATE(DISPLACEMENTS)
    CALL ERRORS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of boundary conditions for the equation set.
  SUBROUTINE BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<On exit, a pointer to the created boundary conditions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%BOUNDARY_CONDITIONS)) THEN
        CALL FLAG_ERROR("Boundary conditions are already associated for the equations set.",ERR,ERROR,*999)        
      ELSE
        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
          CALL FLAG_ERROR("Boundary conditions is already associated.",ERR,ERROR,*999)
        ELSE
          !Initialise the boundary conditions
          CALL BOUNDARY_CONDITIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
          !Return the pointer
          BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_CREATE_START")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_CREATE_START")
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys boundary conditions
  SUBROUTINE BOUNDARY_CONDITIONS_DESTROY(BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      CALL BOUNDARY_CONDITIONS_FINALISE(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_DESTROY")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_DESTROY")
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions and deallocate all memory.
  SUBROUTINE BOUNDARY_CONDITIONS_FINALISE(BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_type_idx
    
    CALL ENTERS("BOUNDARY_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP)) THEN
        DO variable_type_idx=1,SIZE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP,1)
          CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
            & variable_type_idx)%PTR,ERR,ERROR,*999)
        ENDDO !varaible_type_idx
        DEALLOCATE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP)
      ENDIF
      DEALLOCATE(BOUNDARY_CONDITIONS)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_FINALISE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_FINALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_FINALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the boundary conditions for an equations set.
  SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,linear_variable_idx,linear_variable_type,variable_type_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING    
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%BOUNDARY_CONDITIONS)) THEN
        CALL FLAG_ERROR("Boundary conditions is already associated for this equations set.",ERR,ERROR,*998)
      ELSE
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          IF(EQUATIONS%EQUATIONS_FINISHED) THEN
            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
            IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
              IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
                ALLOCATE(EQUATIONS_SET%BOUNDARY_CONDITIONS,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate boundary conditions.",ERR,ERROR,*999)
                EQUATIONS_SET%BOUNDARY_CONDITIONS%EQUATIONS_SET=>EQUATIONS_SET
                EQUATIONS_SET%BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.FALSE.
                ALLOCATE(EQUATIONS_SET%BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(FIELD_NUMBER_OF_VARIABLE_TYPES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate boundary conditions variable type map.",ERR,ERROR,*999)
                DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                  NULLIFY(EQUATIONS_SET%BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type_idx)%PTR)
                ENDDO !variable_type_idx
                SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
                CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                  SELECT CASE(EQUATIONS%LINEARITY)
                  CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                    LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                      DO linear_variable_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                        linear_variable_type=LINEAR_MAPPING%LINEAR_MATRIX_VARIABLE_TYPES(linear_variable_idx)
                        IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(linear_variable_type)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,LINEAR_MAPPING% &
                            & VAR_TO_EQUATIONS_MATRICES_MAPS(linear_variable_type)%VARIABLE,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !linear_variable_idx
                    ELSE
                      CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,RHS_MAPPING% &
                        & RHS_VARIABLE,ERR,ERROR,*999)
                    ENDIF
                  CASE(EQUATIONS_NONLINEAR)
                    NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,NONLINEAR_MAPPING% &
                        & RESIDUAL_VARIABLE,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,RHS_MAPPING% &
                        & RHS_VARIABLE,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                  SELECT CASE(EQUATIONS%LINEARITY)
                  CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                    DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                    IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,DYNAMIC_MAPPING% &
                        & DYNAMIC_VARIABLE,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,RHS_MAPPING% &
                        & RHS_VARIABLE,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE(EQUATIONS_NONLINEAR)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                CASE DEFAULT
                  LOCAL_ERROR="The equations time dependence type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Equations mapping has not been finished.",ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations has not been finished.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_FINALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_INITIALISE")
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE
  
  !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified constant.
  SUBROUTINE BOUNDARY_CONDITIONS_SET_CONSTANT(BOUNDARY_CONDITIONS,VARIABLE_TYPE,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_CONSTANT",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD           
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              CALL FIELD_COMPONENT_DOF_GET_CONSTANT(DEPENDENT_FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,local_ny,global_ny, &
                & ERR,ERROR,*999)
              BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                SELECT CASE(CONDITION)
                CASE(EQUATIONS_SET_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED
                CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                    & VALUE,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The specified boundary condition type of "//TRIM(NUMBER_TO_VSTRING(CONDITION,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " has not been created."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The equations set dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The equations set dependent field has not been finished.",ERR,ERROR,*999)              
          ENDIF
        ELSE
          CALL FLAG_ERROR("The boundary conditions equations set is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITION_SET_CONSTANT")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_SET_CONSTANT",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_SET_CONSTANT")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_CONSTANT
  
  !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified DOF.
  SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1(BOUNDARY_CONDITIONS,VARIABLE_TYPE,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOFS(BOUNDARY_CONDITIONS,VARIABLE_TYPE,(/DOF_INDEX/),(/CONDITION/),(/VALUE/),ERR,ERROR,*999)
           
    CALL EXITS("BOUNDARY_CONDITION_SET_LOCAL_DOF1")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_SET_LOCAL_DOF1",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_SET_LOCAL_DOF1")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1
  
  !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified DOFs.
  SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS(BOUNDARY_CONDITIONS,VARIABLE_TYPE,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDICES(:) !<DOF_INDICES(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITIONS(:) !<CONDITIONS(:). The boundary condition type to set for the i'th dof \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(:). The value of the boundary condition for the i'th dof to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,global_ny,local_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              NULLIFY(DEPENDENT_VARIABLE)
              CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,VARIABLE_TYPE,DEPENDENT_VARIABLE,ERR,ERROR,*999)
              DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
              IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
                BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
                IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                  IF(SIZE(DOF_INDICES,1)==SIZE(CONDITIONS,1)) THEN
                    IF(SIZE(DOF_INDICES,1)==SIZE(VALUES,1)) THEN
                      DO i=1,SIZE(DOF_INDICES,1)
                        local_ny=DOF_INDICES(i)
                        IF(local_ny>=1.AND.local_ny<=DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                          global_ny=DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                          SELECT CASE(CONDITIONS(i))
                          CASE(EQUATIONS_SET_NOT_FIXED)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED
                          CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)= &
                              & EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            LOCAL_ERROR="The specified boundary condition type for dof index "// &
                              & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" of "// &
                              & TRIM(NUMBER_TO_VSTRING(CONDITIONS(i),"*",ERR,ERROR))//" is invalid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        ELSE 
                          LOCAL_ERROR="The local dof of  "//&
                            & TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))//" at dof index "// &
                            & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))// &
                            & " is invalid. The dof should be between 1 and "// &
                            & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)                            
                        ENDIF
                      ENDDO !i
                    ELSE
                      LOCAL_ERROR="The size of the dof indices array ("// &
                        & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                        & ") does not match the size of the values array ("// &
                        & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The size of the dof indices array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                      & ") does not match the size of the fixed conditions array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(CONDITIONS,1),"*",ERR,ERROR))//")."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Boundary conditions variable is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The dependent field variable domain mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The equations set dependent field is not associated..",ERR,ERROR,*999)              
            ENDIF
          ELSE
            CALL FLAG_ERROR("The equations set dependent field has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The boundary conditions equations set is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITION_SET_LOCAL_DOFS")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_SET_LOCAL_DOFS",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_SET_LOCAL_DOFS")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS
  
  !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified user element.
  SUBROUTINE BOUNDARY_CONDITIONS_SET_ELEMENT(BOUNDARY_CONDITIONS,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_ELEMENT",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD           
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              CALL FIELD_COMPONENT_DOF_GET_USER_ELEMENT(DEPENDENT_FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
                & local_ny,global_ny,ERR,ERROR,*999)
              BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                SELECT CASE(CONDITION)
                CASE(EQUATIONS_SET_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED
                CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                    & VALUE,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The specified boundary condition type of "//TRIM(NUMBER_TO_VSTRING(CONDITION,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " has not been created."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The equations set dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The equations set dependent field has not been finished.",ERR,ERROR,*999)              
          ENDIF
        ELSE
          CALL FLAG_ERROR("The boundary conditions equations set is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITION_SET_ELEMENT")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_SET_ELEMENT",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_SET_ELEMENT")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_ELEMENT
  
  !
  !================================================================================================================================
  !

 
  !>Sets a boundary condition on the specified user node.
  SUBROUTINE BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,VARIABLE_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER, &
    & CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD           
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              CALL FIELD_COMPONENT_DOF_GET_USER_NODE(DEPENDENT_FIELD,VARIABLE_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
                & COMPONENT_NUMBER,local_ny,global_ny,ERR,ERROR,*999)
              BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                SELECT CASE(CONDITION)
                CASE(EQUATIONS_SET_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED
                CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The specified boundary condition type of "//TRIM(NUMBER_TO_VSTRING(CONDITION,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " has not been created."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The equations set dependent field is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The equations set dependent has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Boundary conditions equations set is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITION_SET_NODE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_SET_NODE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_SET_NODE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_NODE
  
  !
  !================================================================================================================================
  !
  
  !>Update the boundary condition from analytic value.
  SUBROUTINE BOUNDARY_CONDITION_PARAM_SET_UPDATE_FROM_ANAL_VALUE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:)
    TYPE(FIELD_TYPE), POINTER :: FIELD 
    INTEGER(INTG) :: var_idx,var_type,comp_idx,node_idx,node_number,dev_idx
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    REAL(DP) :: VALUE
    
    
    CALL ENTERS("BOUNDARY_CONDITION_PARAM_SET_UPDATE_FROM_ANAL_VALUE",ERR,ERROR,*999)
    
    FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
    IF(ASSOCIATED(FIELD)) THEN
      DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
        var_type=FIELD%VARIABLES(var_idx)%VARIABLE_TYPE
        CALL FIELD_PARAMETER_SET_DATA_GET(FIELD,var_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ANALYTIC_PARAMETERS,ERR,ERROR,*999)
        DO comp_idx=1,FIELD%VARIABLES(var_idx)%NUMBER_OF_COMPONENTS
          DOMAIN_NODES=>FIELD%VARIABLES(1)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
          IF(ASSOCIATED(DOMAIN_NODES)) THEN
            NODES_MAPPING=>FIELD%VARIABLES(1)%COMPONENTS(comp_idx)%DOMAIN%MAPPINGS%NODES
            DO node_idx=NODES_MAPPING%INTERNAL_START,NODES_MAPPING%INTERNAL_FINISH
              node_number=NODES_MAPPING%DOMAIN_LIST(node_idx)
              DO dev_idx=1,DOMAIN_NODES%NODES(node_number)%NUMBER_OF_DERIVATIVES
                ! Set the boundary condition for dependent field \todo for rectangular system only
                IF(var_idx==1.AND.DOMAIN_NODES%NODES(node_number)%NUMBER_OF_SURROUNDING_ELEMENTS &
                  & < 2**FIELD%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN   
                  VALUE = ANALYTIC_PARAMETERS(FIELD%VARIABLES(var_idx)%COMPONENTS(comp_idx)%PARAM_TO_DOF_MAP% &
                    & NODE_PARAM2DOF_MAP(dev_idx,node_number))
                  CALL BOUNDARY_CONDITIONS_SET_NODE(EQUATIONS_SET%BOUNDARY_CONDITIONS,var_idx, &
                    & dev_idx,node_number,comp_idx,BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                ENDIF       
              ENDDO ! dev_idx
            ENDDO ! node_idx
          ELSE
            CALL FLAG_ERROR("Domain nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ENDDO ! comp_idx
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(FIELD,var_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ANALYTIC_PARAMETERS,ERR,ERROR,*999)
      ENDDO ! var_idx
    ELSE
      CALL FLAG_ERROR("The field is not associated.",ERR,ERROR,*999) 
    ENDIF
    
    CALL EXITS("BOUNDARY_CONDITION_PARAMSET_UPDATE_FROM_ANAL_VALUE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_PARAM_SET_UPDATE_FROM_ANAL_VALUE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_PARAM_SET_UPDATE_FROM_ANAL_VALUE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITION_PARAM_SET_UPDATE_FROM_ANAL_VALUE


  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions variable and deallocate all memory.
  SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE !<A pointer to the boundary conditions variable to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ALLOCATED(BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS))  &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS)
      !IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%BOUNDARY_CONDITIONS_VALUES)) &
      !  & CALL DISTRIBUTED_VECTOR_DESTROY(BOUNDARY_CONDITIONS_VARIABLE%BOUNDARY_CONDITIONS_VALUES,ERR,ERROR,*999)
      DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the boundary conditions variable for a variable type.
  SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(BOUNDARY_CONDITIONS,FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to initialise a variable t ype for
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable to initialise the boundary conditions variable for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_type
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: VARIABLE_DOMAIN_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP)) THEN
          VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
          IF(ASSOCIATED(VARIABLE_DOMAIN_MAPPING)) THEN
            variable_type=FIELD_VARIABLE%VARIABLE_TYPE
            IF(ASSOCIATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR)) THEN
              LOCAL_ERROR="The boundary conditions variable is already associated for variable type "// &
                & TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
            ELSE
              ALLOCATE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate boundary condition variable.",ERR,ERROR,*999)
              BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%BOUNDARY_CONDITIONS=> &
                & BOUNDARY_CONDITIONS
              BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%VARIABLE_TYPE=variable_type
              BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%VARIABLE=>FIELD_VARIABLE            
              ALLOCATE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR% &
                & GLOBAL_BOUNDARY_CONDITIONS(VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global boundary conditions.",ERR,ERROR,*999)
              BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%GLOBAL_BOUNDARY_CONDITIONS= &
                & BOUNDARY_CONDITION_NOT_FIXED
            ENDIF
          ELSE
            CALL FLAG_ERROR("Field variable domain mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Boundary conditions variable type map is not allocated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
      & variable_type)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for an equations set.
  SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to get the boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<On exit, a pointer to the boundary conditions in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_BOUNDARY_CONDITIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
          CALL FLAG_ERROR("Boundary conditions is already associated.",ERR,ERROR,*999)
        ELSE
          BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
          IF(.NOT.ASSOCIATED(BOUNDARY_CONDITIONS)) CALL FLAG_ERROR("Equations set boundary conditions is not associated.", &
            & ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_GET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_BOUNDARY_CONDITIONS_GET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_GET")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_GET

  !
  !================================================================================================================================
  !  
 
END MODULE BOUNDARY_CONDITIONS_ROUTINES
