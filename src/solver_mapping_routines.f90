!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all solver mapping routines.
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

!> This module handles all solver mapping routines.
MODULE SOLVER_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE COMP_ENVIRONMENT
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup SOLVER_MAPPING_EquationsMatrixTypes SOLVER_MAPPING::EquationsMatrixTypes
  !> \brief Equations matrix types
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX=1 !<The equations matrix in the solver mapping is a dynamic equations matrix \see SOLVER_MAPPING_EquationsMatrixTypes,SOLVER_MAPPING
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX=2 !<The equations matrix in the solver mapping is a linear equations matrix \see SOLVER_MAPPING_EquationsMatrixTypes,SOLVER_MAPPING
 !>@}
 
  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX,SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX

  PUBLIC SOLVER_MAPPING_CREATE_FINISH,SOLVER_MAPPING_CREATE_START,SOLVER_MAPPING_DESTROY, &
    & SOLVER_MAPPING_EQUATIONS_SET_ADD,SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET, &
    & SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the solver mappings
  SUBROUTINE SOLVER_MAPPING_CALCULATE(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,dof_idx,DEPENDENT_VARIABLE_TYPE,DYNAMIC_EQUATIONS_MATRIX_OFFSET,equations_column, &
      & equations_matrix,equations_matrix_idx,equations_set_idx,global_dof,global_row,jacobian_column,local_dof, &
      & LINEAR_EQUATIONS_MATRIX_OFFSET,local_row,matrix_number,myrank,myrank_local_dof,NUMBER_OF_COLUMNS, &
      & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,NUMBER_OF_LINEAR_EQUATIONS_MATRICES,NUMBER_OF_GLOBAL_SOLVER_COLS, &
      & LOCAL_SOLVER_DOF_OFFSET,NUMBER_OF_GHOST_SOLVER_DOFS,NUMBER_OF_GLOBAL_SOLVER_ROWS,NUMBER_OF_LOCAL_SOLVER_COLS, &
      & NUMBER_OF_LOCAL_SOLVER_DOFS,NUMBER_OF_LOCAL_SOLVER_ROWS,NUMBER_OF_VARIABLES,rank,rank_idx,row_idx,SOLVER_DOF, &
      & solver_matrix_idx,TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS,variable_idx,variable_type,equations_row_number
    LOGICAL :: INCLUDE_ROW,MYRANK_DOF,RANK_DOF
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COL_DOMAIN_MAPPING,COL_DOFS_MAPPING,ROW_DOMAIN_MAPPING,ROW_DOFS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MAPPING_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
        SOLVER_EQUATIONS=>SOLVER_MAPPING%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN          
          !
          !--- Row mappings ---
          !
          !Calculate the row mappings.
          !We do not have any couplings defined at the moment there is only a 1-1 mapping.
          myrank=COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER
          NUMBER_OF_GLOBAL_SOLVER_ROWS=0
          NUMBER_OF_LOCAL_SOLVER_ROWS=0                                 
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING                
                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                  DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                  LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                  NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                  RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                  SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
                  ROW_DOFS_MAPPING=>EQUATIONS_MAPPING%ROW_DOFS_MAPPING
                  IF(ASSOCIATED(ROW_DOFS_MAPPING)) THEN                    
                    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                      BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN                      
!!TODO: see how slow this is. At the moment we go through number of ranks*number of global rows. We could presort the global rows into a list for each rank. This would take additional memory.
                        DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
                          DO global_row=1,EQUATIONS_MAPPING%NUMBER_OF_GLOBAL_ROWS
                            RANK_DOF=.FALSE.
                            local_row=0
                            DO rank_idx=1,ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%NUMBER_OF_DOMAINS
                              IF(ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%DOMAIN_NUMBER(rank_idx)==rank &
                                & .AND.ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%LOCAL_TYPE(rank_idx)/= &
                                & DOMAIN_LOCAL_GHOST) THEN
                                RANK_DOF=.TRUE.
                                EXIT
                              ENDIF
                            ENDDO !rank_idx
                            IF(RANK_DOF) THEN
                              INCLUDE_ROW=.FALSE.
                              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                                DEPENDENT_VARIABLE_TYPE=DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE
                                BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                                  & DEPENDENT_VARIABLE_TYPE)%PTR
                                IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                  !This is wrong as we only have the mappings for the local rank not the global ranks.
                                  !For now assume 1-1 mapping between rows and dofs.
                                  global_dof=global_row                                  
                                  INCLUDE_ROW=BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)== &
                                    & BOUNDARY_CONDITION_NOT_FIXED
                                ELSE
                                  CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                                  DEPENDENT_VARIABLE_TYPE=NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                                  BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                                    & DEPENDENT_VARIABLE_TYPE)%PTR
                                  IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                    !This is wrong as we only have the mappings for the local rank not the global ranks.
                                    !For now assume 1-1 mapping between rows and dofs.
                                    global_dof=global_row                                  
                                    INCLUDE_ROW=BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)== &
                                      & BOUNDARY_CONDITION_NOT_FIXED
                                  ELSE
                                    CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                  !Loop over the variables in the equations set. Don't include the row in the solver matrices if
                                  !all the variable dofs associated with this equations row are fixed.
                                  DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                                    DEPENDENT_VARIABLE_TYPE=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)% &
                                      & VARIABLE_TYPE
                                    BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                                      & DEPENDENT_VARIABLE_TYPE)%PTR
                                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                      !This is wrong as we only have the mappings for the local rank not the global ranks.
                                      !For now assume 1-1 mapping between rows and dofs.
                                      !
                                      !local_dof=EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(local_row)% &
                                      !  & ROW_TO_DOFS_MAP(equations_matrix_idx)
                                      !variable_type=EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES(equations_matrix_idx)`
                                      !global_dof=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_dof)
                                      global_dof=global_row                                    
                                      INCLUDE_ROW=INCLUDE_ROW.OR.BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS( &
                                        & global_dof)==BOUNDARY_CONDITION_NOT_FIXED
                                    ELSE
                                      CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ENDDO !matrix_idx
                                ENDIF
                              ENDIF
                              IF(INCLUDE_ROW) THEN
                                NUMBER_OF_GLOBAL_SOLVER_ROWS=NUMBER_OF_GLOBAL_SOLVER_ROWS+1
                                IF(rank==myrank) NUMBER_OF_LOCAL_SOLVER_ROWS=NUMBER_OF_LOCAL_SOLVER_ROWS+1 !1-1 mapping
                              ENDIF !include row
                            ENDIF !rank dof
                          ENDDO !global_row
                        ENDDO !rank
                      ELSE
                        CALL FLAG_ERROR("Equations set boundary conditions is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set row degree of freedom mappings is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations equations mapping is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !equations set idx
          IF(NUMBER_OF_LOCAL_SOLVER_ROWS==0) &
            & CALL FLAG_ERROR("Invalid problem setup. The number of local solver rows is zero.",ERR,ERROR,*999)
          IF(NUMBER_OF_GLOBAL_SOLVER_ROWS==0) &
            & CALL FLAG_ERROR("Invalid problem setup. The number of global solver rows is zero.",ERR,ERROR,*999)
          !Allocate memory for the rows mapping
          !Allocate equations set to solver map
          ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping equations set to solver map.",ERR,ERROR,*999)      
          !Allocate the solver rows to equations set maps
          ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping solver row to equations set map.",ERR,ERROR,*999)
          !Set the number of rows
          SOLVER_MAPPING%NUMBER_OF_ROWS=NUMBER_OF_LOCAL_SOLVER_ROWS
          SOLVER_MAPPING%NUMBER_OF_GLOBAL_ROWS=NUMBER_OF_GLOBAL_SOLVER_ROWS
          !Allocate the solver rows domain mapping
          ALLOCATE(SOLVER_MAPPING%ROW_DOFS_MAPPING,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping row dofs mapping.",ERR,ERROR,*999)
!!TODO: what is the real number of domains for a solver???
          CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(SOLVER_MAPPING%ROW_DOFS_MAPPING,COMPUTATIONAL_ENVIRONMENT% &
            & NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)
          ROW_DOMAIN_MAPPING=>SOLVER_MAPPING%ROW_DOFS_MAPPING
          ALLOCATE(ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row dofs mapping global to local map.",ERR,ERROR,*999)
          ROW_DOMAIN_MAPPING%NUMBER_OF_GLOBAL=NUMBER_OF_GLOBAL_SOLVER_ROWS              
          !Calculate the row mappings
          NUMBER_OF_GLOBAL_SOLVER_ROWS=0
          !Loop over the ranks to  ensure that the lowest ranks have the lowest numbered solver variables
          DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            NUMBER_OF_LOCAL_SOLVER_ROWS=0
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              !Note that pointers have been checked for association above
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
              ROW_DOFS_MAPPING=>EQUATIONS_MAPPING%ROW_DOFS_MAPPING
              DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
              IF(rank==myrank) THEN
                !Initialise the equations set to solver map
                CALL SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                  & equations_set_idx),ERR,ERROR,*999)
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_SET_INDEX=equations_set_idx
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%SOLVER_MAPPING=>SOLVER_MAPPING
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS=>EQUATIONS
                !Allocate the equations set to solver maps for solver matrix (sm) indexing
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                  & SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set to solver map equations to solver matrix maps sm.", &
                  & ERR,ERROR,*999)
                DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                  CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE(SOLVER_MAPPING% &
                    & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx),ERR,ERROR,*999)
                ENDDO !solver_matrix_idx
                IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                  !Allocate the equations set to solver maps for equations matrix (em) indexing
                  ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                    & DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
                  IF(ERR/=0) &
                    & CALL FLAG_ERROR("Could not allocate equations set to solver map equations to solver matrix maps em.", &
                    & ERR,ERROR,*999)
                  DO equations_matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                    CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE(SOLVER_MAPPING% &
                      & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                      & equations_matrix_idx),ERR,ERROR,*999)
                  ENDDO !equations_matrix_idx
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
                    CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE(SOLVER_MAPPING% &
                      & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx),ERR,ERROR,*999)
                  ENDIF
                ELSE
                  IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                    !Allocate the equations set to solver maps for equations matrix (em) indexing
                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                      & LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                    IF(ERR/=0) &
                      & CALL FLAG_ERROR("Could not allocate equations set to solver map equations to solver matrix maps em.", &
                      & ERR,ERROR,*999)
                    DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                      CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE(SOLVER_MAPPING% &
                        & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                        & equations_matrix_idx),ERR,ERROR,*999)
                    ENDDO !equations_matrix_idx
                  ENDIF
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
                    CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE(SOLVER_MAPPING% &
                      & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx),ERR,ERROR,*999)
                  ENDIF
                ENDIF
                !Allocate the equations row to solver rows maps
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                  & EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set to solver map equations row to solver rows maps.", &
                  & ERR,ERROR,*999)
                DO equations_row_number=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
                  !Initialise
                  CALL SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE(SOLVER_MAPPING% &
                    & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number), &
                    & ERR,ERROR,*999)
                ENDDO
              ENDIF !rank==my rank
!!TODO: see how slow this is. At the moment we go through number of ranks*number of global rows. We could presort the global rows into a list for each rank. This would take additional memory.
              DO global_row=1,ROW_DOFS_MAPPING%NUMBER_OF_GLOBAL
                RANK_DOF=.FALSE.
                local_row=0
                DO rank_idx=1,ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%NUMBER_OF_DOMAINS
                  IF(ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%DOMAIN_NUMBER(rank_idx)==rank &
                    & .AND.ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%LOCAL_TYPE(rank_idx)/= &
                    & DOMAIN_LOCAL_GHOST) THEN
                    RANK_DOF=.TRUE.
                    local_row=ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%LOCAL_NUMBER(rank_idx)
                    EXIT
                  ENDIF
                ENDDO !rank_idx                            
                IF(RANK_DOF) THEN
                  INCLUDE_ROW=.FALSE.
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    DEPENDENT_VARIABLE_TYPE=DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE
                    BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                      & DEPENDENT_VARIABLE_TYPE)%PTR
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                      !This is wrong as we only have the mappings for the local rank not the global ranks.
                      !For now assume 1-1 mapping between rows and dofs.
                      global_dof=global_row                                  
                      INCLUDE_ROW=BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)==BOUNDARY_CONDITION_NOT_FIXED
                    ELSE
                      CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      DEPENDENT_VARIABLE_TYPE=NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                      BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                        & DEPENDENT_VARIABLE_TYPE)%PTR
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                        !This is wrong as we only have the mappings for the local rank not the global ranks.
                        !For now assume 1-1 mapping between rows and dofs.
                        global_dof=global_row                                  
                        INCLUDE_ROW=BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)== &
                          & BOUNDARY_CONDITION_NOT_FIXED
                      ELSE
                        CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                      DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                        DEPENDENT_VARIABLE_TYPE=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)% &
                          & VARIABLE_TYPE
                        BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                          & DEPENDENT_VARIABLE_TYPE)%PTR
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          !This is wrong as we only have the mappings for the local rank not the global ranks.
                          !For now assume 1-1 mapping between rows and dofs.
                          !
                          !local_dof=EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(local_row)% &
                          !  & ROW_TO_DOFS_MAP(equations_matrix_idx)
                          !variable_type=EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES(equations_matrix_idx)`
                          !global_dof=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_dof)
                          global_dof=global_row                                    
                          INCLUDE_ROW=INCLUDE_ROW.OR.BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)== &
                            & BOUNDARY_CONDITION_NOT_FIXED
                        ELSE
                          CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO !matrix_idx
                    ENDIF
                  ENDIF
                  IF(INCLUDE_ROW) THEN
                    NUMBER_OF_GLOBAL_SOLVER_ROWS=NUMBER_OF_GLOBAL_SOLVER_ROWS+1
                    NUMBER_OF_LOCAL_SOLVER_ROWS=NUMBER_OF_LOCAL_SOLVER_ROWS+1
                    !Set up the row domain mappings.
                    !There are no ghosted rows for the solver matrices so there is only one domain for the global to local map.
                    !Initialise
                    CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                      & NUMBER_OF_GLOBAL_SOLVER_ROWS),ERR,ERROR,*999)
                    !Allocate the global to local map arrays
                    ALLOCATE(ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS)%LOCAL_NUMBER(1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row global to local map local number.",ERR,ERROR,*999)
                    ALLOCATE(ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS)%DOMAIN_NUMBER(1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row global to local map domain number.",ERR,ERROR,*999)
                    ALLOCATE(ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS)%LOCAL_TYPE(1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row global to local map local type.",ERR,ERROR,*999)
                    !Set the global to local mappings
                    ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS)%NUMBER_OF_DOMAINS=1
                    ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS)%LOCAL_NUMBER(1)= &
                      & NUMBER_OF_LOCAL_SOLVER_ROWS
                    ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS)%DOMAIN_NUMBER(1)=rank
                    ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                    !If this is my rank then set up the solver->equations and equations->solver row mappings
                    IF(rank==myrank) THEN
                      !Set up the solver row -> equations row mappings. 1-1 mapping as no coupling at the moment
                      !Initialise
                      CALL SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE(SOLVER_MAPPING% &
                        & SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS),ERR,ERROR,*999)
                      !Allocate the solver row to equations row mapping arrays
                      ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_SET(1), &
                        & STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row to equations set map equations set.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)% &
                        & EQUATIONS_ROW_NUMBER(1),STAT=ERR)
                      IF(ERR/=0) &
                        & CALL FLAG_ERROR("Could not allocate row to equations set map equations row number.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)% &
                        & COUPLING_COEFFICIENTS(1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row to equations set map coupling coefficients.", &
                        & ERR,ERROR,*999)
                      !Set the mappings
                      SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%NUMBER_OF_ROWS=1
                      SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_SET(1)= &
                        & equations_set_idx
                      SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_ROW_NUMBER(1)= &
                        & local_row
                      SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%COUPLING_COEFFICIENTS(1)= &
                        & 1.0_DP
                      !Set up the equations row -> solver row mappings
                      !Allocate the equations row to solver row mappings arrays
                      ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%SOLVER_ROWS(1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to solver rows maps solver rows.", &
                        & ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to solver rows maps solver rows.", &
                        & ERR,ERROR,*999)
                      !Set the mappings
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%NUMBER_OF_SOLVER_ROWS=1
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%SOLVER_ROWS(1)=NUMBER_OF_LOCAL_SOLVER_ROWS
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%COUPLING_COEFFICIENTS(1)=1.0_DP
                    ENDIF !rank==my rank
                  ELSE
                    IF(rank==myrank) THEN
                      !Set up the equations row -> solver row mappings
                      !Set the mappings
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%NUMBER_OF_SOLVER_ROWS=0
                    ENDIF !rank==my rank
                  ENDIF !include row
                ENDIF !rank dof
              ENDDO !global_row
            ENDDO !equations_set_idx
          ENDDO !rank
          CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(ROW_DOMAIN_MAPPING,ERR,ERROR,*999)
          !
          !--- Column mappings ---
          !
          !Allocate solver column to equations sets mapping array
          ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping solver column to equations sets map.",ERR,ERROR,*999)
          !Calculate the column mappings for each solver matrix
          DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
            !Calculate the number of solver variables
            NUMBER_OF_GLOBAL_SOLVER_COLS=0
            NUMBER_OF_LOCAL_SOLVER_COLS=0
            TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS=0
            !Initialise solver column to equations sets mapping array
            CALL SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE(SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx),ERR,ERROR,*999)
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_MATRIX_NUMBER=solver_matrix_idx
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_MAPPING=>SOLVER_MAPPING
            !Allocate the solver col to equations set maps array
            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver col to equations sets map solver col to equation set maps.", &
              & ERR,ERROR,*999)
            !Loop over the ranks
            DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
              !Loop over the equations sets
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                !The pointers below have been checked for association above.
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
                IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                  NUMBER_OF_VARIABLES=1
                ELSE
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    NUMBER_OF_VARIABLES=1
                  ELSE
                    NUMBER_OF_VARIABLES=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,equations_set_idx, &
                      & solver_matrix_idx)
                  ENDIF
                ENDIF
                IF(rank==0) THEN
                  !Initialise equations set to solver map (sm)
                  CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE(SOLVER_MAPPING% &
                    & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx), &
                    & ERR,ERROR,*999)
                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                  !Allocate the equations set to solver map variables arrays
                  ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%VARIABLE_TYPES(NUMBER_OF_VARIABLES),STAT=ERR)
                  IF(ERR/=0)  &
                    & CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm variable types.",ERR,ERROR,*999)
                  ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%VARIABLES(NUMBER_OF_VARIABLES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm variables.",ERR,ERROR,*999)
                  ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(NUMBER_OF_VARIABLES),STAT=ERR)
                  IF(ERR/=0)  &
                    & CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm variables to solver col maps.", &
                    & ERR,ERROR,*999)                  
                  !Setup
                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%NUMBER_OF_VARIABLES=NUMBER_OF_VARIABLES
                  NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
                  NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                ENDIF !rank==0
!!TODO: see how slow this is. At the moment we go through number of ranks*number of globals. We could presort the global nys into a list for each rank. This would take additional memory.
                DO variable_idx=1,NUMBER_OF_VARIABLES
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(equations_set_idx)
                  ELSE
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(equations_set_idx)
                    ELSE
                      variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,equations_set_idx, &
                        & solver_matrix_idx)
                    ENDIF
                  ENDIF
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(COL_DOFS_MAPPING)) THEN
                    BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                      IF(rank==0) THEN
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%VARIABLE_TYPES(variable_idx)=variable_type
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%VARIABLES(variable_idx)%PTR=>DEPENDENT_VARIABLE
                        !Allocate the variable to solver col maps arrays
                        ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)% &
                          & COLUMN_NUMBERS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps column numbers.", &
                          & ERR,ERROR,*999)
                        ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)% &
                          & COUPLING_COEFFICIENTS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps coupling coefficients.", &
                          & ERR,ERROR,*999)
                        ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)% &
                          & ADDITIVE_CONSTANTS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps additive constants.", &
                          & ERR,ERROR,*999)
                        !Setup
                        IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                          NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+DYNAMIC_MAPPING% &
                            & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES
                        ELSE
                          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                            NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES+LINEAR_MAPPING% &
                              & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES
                          ENDIF
                        ENDIF
                      ENDIF !rank==0
                      DO global_dof=1,DEPENDENT_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                        RANK_DOF=.FALSE.
                        DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%NUMBER_OF_DOMAINS
                          IF(COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%DOMAIN_NUMBER(rank_idx)==rank &
                            & .AND.COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_TYPE(rank_idx)/= &
                            & DOMAIN_LOCAL_GHOST) THEN
                            RANK_DOF=.TRUE.
                            local_dof=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_NUMBER(rank_idx)
                            EXIT
                          ENDIF
                        ENDDO !rank_idx                                          
                        IF(RANK_DOF) THEN
                          MYRANK_DOF=.FALSE.
                          DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%NUMBER_OF_DOMAINS
                            IF(COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%DOMAIN_NUMBER(rank_idx)== &
                              & myrank) THEN
                              MYRANK_DOF=.TRUE.
                              EXIT
                            ENDIF
                          ENDDO
                          IF(BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)==BOUNDARY_CONDITION_NOT_FIXED) THEN
                            NUMBER_OF_GLOBAL_SOLVER_COLS=NUMBER_OF_GLOBAL_SOLVER_COLS+1
                            IF(MYRANK_DOF) TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS=TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS+1
                            IF(rank==myrank) NUMBER_OF_LOCAL_SOLVER_COLS=NUMBER_OF_LOCAL_SOLVER_COLS+1
                          ENDIF !global dof not fixed
                        ENDIF !rank dof
                      ENDDO !global_dof
                    ELSE
                      CALL FLAG_ERROR("Boundary condition variable not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations matrix columns degree of freedom mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !variable_idx
              ENDDO !equations_set_idx
            ENDDO !rank
            IF(NUMBER_OF_LOCAL_SOLVER_COLS==0) THEN
              LOCAL_ERROR="Invalid problem setup. The number of local solver columns for solver matrix "// &
                & TRIM(NUMBER_TO_VSTRING(solver_matrix_idx,"*",ERR,ERROR))//" is zero."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
            IF(NUMBER_OF_GLOBAL_SOLVER_COLS==0) THEN
              LOCAL_ERROR="Invalid problem setup. The number of global solver columns for solver matrix "// &
                & TRIM(NUMBER_TO_VSTRING(solver_matrix_idx,"*",ERR,ERROR))//" is zero."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
            !Allocate memory for this solver matrix
            !Allocate solver columns to equations sets maps
            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
              & TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps.",ERR,ERROR,*999)
            !Set the number of columns
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS=NUMBER_OF_GLOBAL_SOLVER_COLS
            !Set the number of variables
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_DOFS=NUMBER_OF_LOCAL_SOLVER_COLS
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS= &
              & TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_GLOBAL_DOFS=NUMBER_OF_GLOBAL_SOLVER_COLS
            !Allocate the columns domain mapping
            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%COLUMN_DOFS_MAPPING,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver col to equations sets map column dofs mapping.",ERR,ERROR,*999)
!!TODO: what is the real number of domains for a solver???
            CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
              & COLUMN_DOFS_MAPPING,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)            
            COL_DOMAIN_MAPPING=>SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%COLUMN_DOFS_MAPPING
            ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column dofs mapping global to local.",ERR,ERROR,*999)
            COL_DOMAIN_MAPPING%NUMBER_OF_GLOBAL=NUMBER_OF_GLOBAL_SOLVER_COLS
            !Calculate the column mappings
            NUMBER_OF_COLUMNS=NUMBER_OF_GLOBAL_SOLVER_COLS
            NUMBER_OF_GLOBAL_SOLVER_COLS=0
            LOCAL_SOLVER_DOF_OFFSET=NUMBER_OF_LOCAL_SOLVER_COLS
            !Loop over the ranks to ensure that the lowest ranks have the lowest numbered solver variables
!!TODO: see how slow this is. At the moment we go through number of ranks*number of globals. We could presort the global nys into a list for each rank. This would take additional memory.
            NUMBER_OF_LOCAL_SOLVER_DOFS=0
            NUMBER_OF_GHOST_SOLVER_DOFS=0
            DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
              NUMBER_OF_LOCAL_SOLVER_COLS=0
              TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS=0
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                !The pointers below have been checked for association above.
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
                IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                  NUMBER_OF_VARIABLES=1
                ELSE
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    NUMBER_OF_VARIABLES=1
                  ELSE
                    NUMBER_OF_VARIABLES=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,equations_set_idx, &
                      & solver_matrix_idx)
                  ENDIF
                ENDIF
                IF(rank==0) THEN
                  !Allocate memory
                  !Initialise solver columns to equations set map
                  CALL SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE(SOLVER_MAPPING% &
                    & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                    & equations_set_idx),ERR,ERROR,*999)
                  !Allocate the solver columns to equations set map arrays
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                      & equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(NUMBER_OF_COLUMNS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver columns to dynamic equations map.",ERR,ERROR,*999)
                  ELSE
                    ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                      & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(NUMBER_OF_COLUMNS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver columns to static equations map.",ERR,ERROR,*999)
                  ENDIF
                  !Set the solver column to equations set map
                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                    & equations_set_idx)%EQUATIONS=>EQUATIONS
                  !Allocate the equations to solver matrix maps sm equations to solver maps
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm dynamic equations "// &
                      & "to solver matrix maps.",ERR,ERROR,*999)
                    !Set up dynamic arrays
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                    DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                      NULLIFY(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR)
                      ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR,STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps.",ERR,ERROR,*999)
                      CALL SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                        & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                        & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR,ERR,ERROR,*999)
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                        & EQUATIONS_MATRIX_TYPE=SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                        & equations_matrix_idx)%EQUATIONS_MATRIX_NUMBER=equations_matrix_idx
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                        & equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                        & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES+1
                    ENDDO !equations_matrix_idx                    
                    !Set up nonlinear arrays
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP,STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Jacobian to solver matrix maps.",ERR,ERROR,*999)
                      CALL SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                        & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP, &
                        & ERR,ERROR,*999)
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM% &
                        & JACOBIAN_TO_SOLVER_MATRIX_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
                    ENDIF
                  ELSE
                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                    IF(ERR/=0) &
                      & CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm equations to solver matrix maps.", &
                      & ERR,ERROR,*999)
                    !Set up linear arrays
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                      NULLIFY(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR)
                      ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR,STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps.",ERR,ERROR,*999)
                      CALL SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                        & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                        & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR,ERR,ERROR,*999)
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                        & EQUATIONS_MATRIX_TYPE=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                        & equations_matrix_idx)%EQUATIONS_MATRIX_NUMBER=equations_matrix_idx
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                        & equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                        & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES+1
                    ENDDO !equations_matrix_idx
                    !Set up nonlinear arrays
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP,STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Jacobian to solver matrix maps.",ERR,ERROR,*999)
                      CALL SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                        & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP, &
                        & ERR,ERROR,*999)
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM% &
                        & JACOBIAN_TO_SOLVER_MATRIX_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
                    ENDIF
                  ENDIF
                ENDIF !rank==0
                !Loop over the variables
                DYNAMIC_EQUATIONS_MATRIX_OFFSET=0
                LINEAR_EQUATIONS_MATRIX_OFFSET=0
                DO variable_idx=1,NUMBER_OF_VARIABLES
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(equations_set_idx)
                  ELSE
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(equations_set_idx)
                    ELSE
                      variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,equations_set_idx, &
                        & solver_matrix_idx)
                    ENDIF
                  ENDIF
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                      & NUMBER_OF_EQUATIONS_MATRICES
                  ELSE
                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                      NUMBER_OF_LINEAR_EQUATIONS_MATRICES=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                        & NUMBER_OF_EQUATIONS_MATRICES
                    ENDIF
                  ENDIF
                  IF(rank==0) THEN
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%VARIABLES(variable_idx)%PTR=>DEPENDENT_VARIABLE
                    IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                      !Allocate dynamic equations to solver matrix maps equations column to solver columns maps
                      DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                        MATRIX_NUMBER=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                          & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(DYNAMIC_EQUATIONS_MATRIX_OFFSET+ &
                          & equations_matrix_idx)%PTR%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(DYNAMIC_EQUATIONS_MATRIX_OFFSET+ &
                          & equations_matrix_idx)%PTR%EQUATIONS_MATRIX_NUMBER=MATRIX_NUMBER
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(DYNAMIC_EQUATIONS_MATRIX_OFFSET+ &
                          & equations_matrix_idx)%PTR%EQUATIONS_MATRIX=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                          & MATRIX_NUMBER)%EQUATIONS_MATRIX
                        NUMBER_OF_COLUMNS=DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                        ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                          & DYNAMIC_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                          & NUMBER_OF_COLUMNS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic equations column to solver columns map.", &
                          & ERR,ERROR,*999)
                      ENDDO !equations_matrix_idx
                      IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_MATRIX=>NONLINEAR_MAPPING% &
                          & JACOBIAN_TO_VAR_MAP%JACOBIAN
                        NUMBER_OF_COLUMNS=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
                        ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                          & JACOBIAN_COL_SOLVER_COLS_MAP(NUMBER_OF_COLUMNS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Jacobian column to solver columns map.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      !Allocate linear equations to solver matrix maps equations column to solver columns maps
                      IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                        DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                            & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                            & equations_matrix_idx)%PTR%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                            & equations_matrix_idx)%PTR%EQUATIONS_MATRIX_NUMBER=MATRIX_NUMBER
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                            & equations_matrix_idx)%PTR%EQUATIONS_MATRIX=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                            & MATRIX_NUMBER)%EQUATIONS_MATRIX
                          NUMBER_OF_COLUMNS=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                          ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & LINEAR_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                            & NUMBER_OF_COLUMNS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear equations column to solver columns map.", &
                            & ERR,ERROR,*999)
                        ENDDO !equations_matrix_idx
                      ENDIF
                      IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                        SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_MATRIX=>NONLINEAR_MAPPING% &
                          & JACOBIAN_TO_VAR_MAP%JACOBIAN
                        NUMBER_OF_COLUMNS=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
                        ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                          & JACOBIAN_COL_SOLVER_COLS_MAP(NUMBER_OF_COLUMNS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Jacobian column to solver columns map.",ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ENDIF !rank==0
                  DO global_dof=1,DEPENDENT_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                    local_dof=0
                    RANK_DOF=.FALSE.
                    DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%NUMBER_OF_DOMAINS
                      IF(COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%DOMAIN_NUMBER(rank_idx)==rank &
                        & .AND.COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_TYPE(rank_idx)/=DOMAIN_LOCAL_GHOST) THEN
                        RANK_DOF=.TRUE.
                        local_dof=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_NUMBER(rank_idx)
                        EXIT
                      ENDIF
                    ENDDO !rank_idx 
                    IF(RANK_DOF) THEN
                      MYRANK_DOF=.FALSE.
                      DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%NUMBER_OF_DOMAINS
                        IF(COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%DOMAIN_NUMBER(rank_idx)==myrank) THEN
                          MYRANK_DOF=.TRUE.
                          myrank_local_dof=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_NUMBER(rank_idx)
                          EXIT
                        ENDIF
                      ENDDO                                          
                      IF(BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)==BOUNDARY_CONDITION_NOT_FIXED) THEN
                        !DOF is not fixed so map the variable/equation dof to a new solver dof
                        NUMBER_OF_GLOBAL_SOLVER_COLS=NUMBER_OF_GLOBAL_SOLVER_COLS+1
                        !Initialise_sm
                        CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                          & NUMBER_OF_GLOBAL_SOLVER_COLS),ERR,ERROR,*999)
                        IF(MYRANK_DOF) THEN
                          TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS=TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS+1
                          IF(rank==myrank) THEN
                            NUMBER_OF_LOCAL_SOLVER_COLS=NUMBER_OF_LOCAL_SOLVER_COLS+1
                            NUMBER_OF_LOCAL_SOLVER_DOFS=NUMBER_OF_LOCAL_SOLVER_DOFS+1
                            SOLVER_DOF=NUMBER_OF_LOCAL_SOLVER_DOFS
                          ELSE
                            NUMBER_OF_GHOST_SOLVER_DOFS=NUMBER_OF_GHOST_SOLVER_DOFS+1
                            SOLVER_DOF=LOCAL_SOLVER_DOF_OFFSET+NUMBER_OF_GHOST_SOLVER_DOFS
                          ENDIF
                          !Set up the column domain mappings.
                          !There are no ghosted coLs for the solver matrices so there is only 1 domain for the global to local map.
                          !Allocate the global to local map arrays
                          ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS)%LOCAL_NUMBER(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column domain global to local map local number", &
                            & ERR,ERROR,*999)
                          ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS)%DOMAIN_NUMBER(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column domain global to local map domain number", &
                            & ERR,ERROR,*999)
                          ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS)%LOCAL_TYPE(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column domain global to local map domain number", &
                            & ERR,ERROR,*999)
                          !Set up the global to local mappings
                          COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS)%NUMBER_OF_DOMAINS=1
                          COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS)%LOCAL_NUMBER(1)= &
                            & TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS
                          COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS)%DOMAIN_NUMBER(1)=rank
                          COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_COLS)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                          !Set up the solver column -> equations column mappings. 1-1 as no coupling yet
                          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                            !Initialise
                            CALL SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_INITIALISE(SOLVER_MAPPING% &
                              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                              & equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(NUMBER_OF_GLOBAL_SOLVER_COLS), &
                              & ERR,ERROR,*999)
                            !Allocate the solver column to equations column mapping arrays
                            !No coupling yet so the number of columns the solver column is mapped to is just the number of matrices
                            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_MATRIX_NUMBERS(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES), &
                              & STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic equations matrix numbers.",ERR,ERROR,*999)
                            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_COL_NUMBERS(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic equations column numbers.",ERR,ERROR,*999)
                            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%COUPLING_COEFFICIENTS(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic equations coupling coefficients.",ERR,ERROR,*999)
                            !Set the mappings
                            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES= &
                              & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                            !Loop over the dynamic equations matrices associated with the variable and set the column maps
                            DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                              !Set the column map
                              MATRIX_NUMBER=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                              equations_column=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                              SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)=MATRIX_NUMBER
                              SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_COL_NUMBERS(equations_matrix_idx)=equations_column
                              SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%COUPLING_COEFFICIENTS(equations_matrix_idx)= &
                                & DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%MATRIX_COEFFICIENT
                            ENDDO !equations_matrix_idx
                          ELSE
                            !Initialise
                            CALL SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_INITIALISE(SOLVER_MAPPING% &
                              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                              & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(NUMBER_OF_GLOBAL_SOLVER_COLS),ERR,ERROR,*999)
                            !Allocate the solver column to equations column mapping arrays
                            !No coupling yet so the number of columns the solver column is mapped to is just the number of matrices
                            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_MATRIX_NUMBERS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES), &
                              & STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear equations matrix numbers.",ERR,ERROR,*999)
                            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_COL_NUMBERS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear equations column numbers.",ERR,ERROR,*999)
                            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%COUPLING_COEFFICIENTS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear equations coupling coefficients.",ERR,ERROR,*999)
                            !Set the mappings
                            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES= &
                              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                              !Loop over the linear equations matrices associated with the variable and set the column maps
                              DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                !Set the column map
                                MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                  & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                  & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                  & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                                  & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)=MATRIX_NUMBER
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                  & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                                  & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_COL_NUMBERS(equations_matrix_idx)=equations_column
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                  & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                                  & NUMBER_OF_GLOBAL_SOLVER_COLS)%COUPLING_COEFFICIENTS(equations_matrix_idx)= &
                                  & LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%MATRIX_COEFFICIENT
                              ENDDO !equations_matrix_idx
                            ENDIF
                            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                              jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(myrank_local_dof)
                              SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%JACOBIAN_COL_NUMBER=jacobian_column
                              SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%JACOBIAN_COUPLING_COEFFICIENT= &
                                & NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT
                            ENDIF
                          ENDIF
                          !Set up the solver dofs -> variable dofs map
                          !Initialise
                          CALL SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE(SOLVER_MAPPING% &
                            & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF),ERR,ERROR,*999)
                          !Allocate the solver dofs to variable dofs arrays
                          !No coupling so there is only one equations set at the moment
                          ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%EQUATIONS_SET_INDICES(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps equations set indices.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%VARIABLE(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable type.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%VARIABLE_DOF(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable dof.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%VARIABLE_COEFFICIENT(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable coefficient.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%ADDITIVE_CONSTANT(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps additive constant.", &
                            & ERR,ERROR,*999)
                          !Setup
                          SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%NUMBER_OF_EQUATIONS_SETS=1
                          SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%EQUATIONS_SET_INDICES(1)=equations_set_idx
                          SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%VARIABLE(1)%PTR=>DEPENDENT_VARIABLE
                          SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%VARIABLE_DOF(1)=myrank_local_dof
                          SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%VARIABLE_COEFFICIENT(1)=1.0_DP
                          SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%ADDITIVE_CONSTANT(1)=0.0_DP
                          !Set up the equations variables -> solver columns mapping
                          !No coupling yet so the mapping is 1-1
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COLUMN_NUMBERS(myrank_local_dof)= &
                            & NUMBER_OF_GLOBAL_SOLVER_COLS
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COUPLING_COEFFICIENTS( &
                            & myrank_local_dof)=1.0_DP
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%ADDITIVE_CONSTANTS( &
                            & myrank_local_dof)=0.0_DP
                          !Set up the equations columns -> solver columns mapping
                          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                            DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                              MATRIX_NUMBER=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                              equations_column=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                              !Allocate the equation to solver map column items.
                              !No coupling yet so the mapping is 1-1
                              ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                                & DYNAMIC_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                & equations_column)%SOLVER_COLS(1),STAT=ERR)
                              IF(ERR/=0) CALL  FLAG_ERROR("Could not allocate dynamic equations column to solver columns map "// &
                                & "solver colums.",ERR,ERROR,*999)
                              ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                                & DYNAMIC_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                & equations_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic equations column to solver columns map "// &
                                & "coupling coefficients.",ERR,ERROR,*999)
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(DYNAMIC_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=1
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(DYNAMIC_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)%SOLVER_COLS(1)= &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(DYNAMIC_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)% &
                                & COUPLING_COEFFICIENTS(1)=DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)% &
                                & MATRIX_COEFFICIENT                              
                            ENDDO !equations_matrix_idx
                          ELSE
                            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                              DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                  & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                  & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                                !Allocate the equation to solver map column items.
                                !No coupling yet so the mapping is 1-1
                                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                                  & LINEAR_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                  & equations_column)%SOLVER_COLS(1),STAT=ERR)
                                IF(ERR/=0) CALL  FLAG_ERROR("Could not allocate linear equations column to solver columns map "// &
                                  & "solver colums.",ERR,ERROR,*999)
                                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                                  & LINEAR_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                  & equations_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear equations column to solver columns map "// &
                                  & "coupling coefficients.",ERR,ERROR,*999)
                                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                  & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                  & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)% &
                                  & NUMBER_OF_SOLVER_COLS=1
                                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                  & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                  & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)%SOLVER_COLS(1)= &
                                  & NUMBER_OF_GLOBAL_SOLVER_COLS
                                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                  & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                  & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)% &
                                  & COUPLING_COEFFICIENTS(1)=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)% &
                                  & MATRIX_COEFFICIENT                              
                              ENDDO !equations_matrix_idx
                            ENDIF
                            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                              jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(myrank_local_dof)
                              !Allocate the Jacobian to solver map column items.
                              !No coupling yet so the mapping is 1-1
                              ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)%SOLVER_COLS(1),STAT=ERR)
                              IF(ERR/=0) CALL  &
                                & FLAG_ERROR("Could not allocate Jacobian column to solver columns map solver colums.", &
                                & ERR,ERROR,*999)
                              ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                              IF(ERR/=0) CALL &
                                & FLAG_ERROR("Could not allocate Jacobain column to solver columns map coupling coefficients.",&
                                & ERR,ERROR,*999)
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)% &
                                & NUMBER_OF_SOLVER_COLS=1
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)% &
                                & SOLVER_COLS(1)=NUMBER_OF_GLOBAL_SOLVER_COLS
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)% &
                                & COUPLING_COEFFICIENTS(1)=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT
                            ENDIF
                          ENDIF
                        ENDIF
                      ELSE
                        IF(MYRANK_DOF) THEN
                          !Set up the equations variables -> solver columns mapping
                          !No coupling yet so the mapping is 1-1
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COLUMN_NUMBERS(myrank_local_dof)=0
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COUPLING_COEFFICIENTS( &
                            & myrank_local_dof)=0.0_DP
                          SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%ADDITIVE_CONSTANTS( &
                            & myrank_local_dof)=0.0_DP
                          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                            !Set up the equations columns -> solver columns mapping
                            DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                              MATRIX_NUMBER=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                              equations_column=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                              !No coupling yet so the mapping is 1-1
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(DYNAMIC_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=0
                            ENDDO !equations_matrix_idx
                          ELSE
                            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                              !Set up the equations columns -> solver columns mapping
                              DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                  & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                  & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                                !No coupling yet so the mapping is 1-1
                                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                  & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                  & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)% &
                                  & NUMBER_OF_SOLVER_COLS=0
                              ENDDO !equations_matrix_idx
                            ENDIF
                            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                              jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(myrank_local_dof)
                              !No coupling yet so the mapping is 1-1
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                & jacobian_column)%NUMBER_OF_SOLVER_COLS=0
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF !field dof is fixed
                    ENDIF !rank dof
                  ENDDO !global_dof
                  DYNAMIC_EQUATIONS_MATRIX_OFFSET=DYNAMIC_EQUATIONS_MATRIX_OFFSET+NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  LINEAR_EQUATIONS_MATRIX_OFFSET=LINEAR_EQUATIONS_MATRIX_OFFSET+NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                ENDDO !variable_idx
              ENDDO !equations_set_idx
            ENDDO !rank            
            CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(COL_DOMAIN_MAPPING,ERR,ERROR,*999)
          ENDDO !solver_matrix_idx
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            EQUATIONS=>EQUATIONS_SET%EQUATIONS
            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
              DO equations_matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                  & equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                  & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps.",ERR,ERROR,*999)
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                  & equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=0
              ENDDO !equations_matrix_idx
              DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                    & SOLVER_MATRIX_NUMBER==solver_matrix_idx) THEN
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                      & equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                      & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES+1
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                      & equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                      & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES)% &
                      & PTR=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR
                  ENDIF
                ENDDO !equations_matrix_idx
              ENDDO !solver_matrix_idx             
            ELSE IF(ASSOCIATED(LINEAR_MAPPING)) THEN
              DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                  & equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                  & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps.",ERR,ERROR,*999)
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                  & equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=0
              ENDDO !equations_matrix_idx
              DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                    & SOLVER_MATRIX_NUMBER==solver_matrix_idx) THEN
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                      & equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                      & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES+1
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                      & equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                      & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES)% &
                      & PTR=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR
                  ENDIF
                ENDDO !equations_matrix_idx
              ENDDO !solver_matrix_idx             
            ENDIF
          ENDDO !equations_set_idx
        ELSE
          CALL FLAG_ERROR("The solver mapping solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver mapping create values cache is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Solver mappings:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of equations sets = ",SOLVER_MAPPING% &
        & NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equations_set_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Region user number = ",EQUATIONS_SET%REGION%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set user number = ",EQUATIONS_SET%USER_NUMBER, &
          & ERR,ERROR,*999)                
      ENDDO !equations_set_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Row mappings:",ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",SOLVER_MAPPING%NUMBER_OF_ROWS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global rows = ",SOLVER_MAPPING%NUMBER_OF_GLOBAL_ROWS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver rows to equations sets rows mappings:",ERR,ERROR,*999)      
      DO row_idx=1,SOLVER_MAPPING%NUMBER_OF_ROWS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver row : ",row_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of rows mapped to = ",SOLVER_MAPPING% &
          & SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%NUMBER_OF_ROWS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)% &
          & NUMBER_OF_ROWS,5,5,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%EQUATIONS_SET, &
          & '("      Equations sets indices  :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)% &
          & NUMBER_OF_ROWS,5,5,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%EQUATIONS_ROW_NUMBER, &
          & '("      Equations row numbers   :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)% &
          & NUMBER_OF_ROWS,5,5,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%COUPLING_COEFFICIENTS, &
          & '("      Coupling coefficients   :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999) 
      ENDDO !row_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equations sets rows to solver rows mappings:",ERR,ERROR,*999)
      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equations_set_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations set rows = ",EQUATIONS_MAPPING% &
          & TOTAL_NUMBER_OF_ROWS,ERR,ERROR,*999)
        DO row_idx=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set row : ",row_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows mapped to = ",SOLVER_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS, &
            & ERR,ERROR,*999)
          IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)% &
            & NUMBER_OF_SOLVER_ROWS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS,5,5,SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%SOLVER_ROWS, &
              & '("        Solver row numbers    :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS,5,5,SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%COUPLING_COEFFICIENTS, &
              & '("        Coupling coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !row_idx
      ENDDO !equations_set_idx            
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Column mappings:",ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of solver matrices = ",SOLVER_MAPPING% &
        & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
      DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver columns to equations sets columns mappings:",ERR,ERROR,*999)        
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",SOLVER_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of DOFs = ",SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
          & solver_matrix_idx)%NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of DOFs = ",SOLVER_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global DOFs = ",SOLVER_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_GLOBAL_DOFS,ERR,ERROR,*999)
        DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Equations set index : ",equations_set_idx,ERR,ERROR,*999)          
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Column mappings:",ERR,ERROR,*999)        
          DO column_idx=1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS           
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver column : ",column_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of dynamic equations matrices mapped to = ", &
              & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES, &
              & ERR,ERROR,*999)
            IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES>0) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                & column_idx)%EQUATIONS_MATRIX_NUMBERS,'("      Equation matrices numbers  :",5(X,I13))','(34X,5(X,I13))', &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                & column_idx)%EQUATIONS_COL_NUMBERS,'("      Equation column numbers    :",5(X,I13))','(34X,5(X,I13))', &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                & column_idx)%COUPLING_COEFFICIENTS,'("      Coupling coefficients      :",5(X,E13.6))','(34X,5(X,E13.6))', &
                & ERR,ERROR,*999)
            ENDIF
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of linear equations matrices mapped to = ", &
              & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES, &
              & ERR,ERROR,*999)
            IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                & column_idx)%EQUATIONS_MATRIX_NUMBERS,'("      Equation matrices numbers  :",5(X,I13))','(34X,5(X,I13))', &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                & column_idx)%EQUATIONS_COL_NUMBERS,'("      Equation column numbers    :",5(X,I13))','(34X,5(X,I13))', &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                & column_idx)%COUPLING_COEFFICIENTS,'("      Coupling coefficients      :",5(X,E13.6))','(34X,5(X,E13.6))', &
                & ERR,ERROR,*999)
            ENDIF
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian column number     = ", &
              & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%JACOBIAN_COL_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian coupling coeff    = ", &
              & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%JACOBIAN_COUPLING_COEFFICIENT,ERR,ERROR,*999)
          ENDDO !column_idx
        ENDDO !equations_set_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Variable mappings:",ERR,ERROR,*999)        
        DO dof_idx=1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS     
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver dof : ",dof_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations sets mapped to = ",SOLVER_MAPPING% &
            & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
            & NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
          IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
            & NUMBER_OF_EQUATIONS_SETS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%EQUATIONS_SET_INDICES, &
              & '("      Equations set indices  :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
!!TODO: write out the variables somehow.
            !CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
            !  & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(column_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING% &
            !  & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%VARIABLE_TYPE, &
            !  & '("      Variable types         :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%VARIABLE_DOF, &
              & '("      Variable dofs          :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
              & VARIABLE_COEFFICIENT,'("      Variable coefficients  :",5(X,E13.6))','(28X,5(X,E13.6))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
              & ADDITIVE_CONSTANT,'("      Additive constants     :",5(X,E13.6))','(28X,5(X,E13.6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !dof_idx
      ENDDO !solver_matrix_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equations sets columns to solver columns mappings:",ERR,ERROR,*999)
      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
        LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Equations set index : ",equations_set_idx,ERR,ERROR,*999)          
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix indexing:",ERR,ERROR,*999)
        DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Dynamic equations matrix columns to solver matrix columns:", &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of dynamic equations matrices = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
              & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,ERR,ERROR,*999)
            DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                & equations_matrix_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations matrix index : ",equations_matrix_idx,ERR,ERROR,*999)
              equations_matrix=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX_NUMBER
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix number = ",equations_matrix,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",EQUATIONS_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
              DO column_idx=1,DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix)%NUMBER_OF_COLUMNS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix column : ",column_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver columns mapped to = ", &
                  & EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS,ERR,ERROR,*999)
                IF(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS>0) THEN
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)% &
                    & SOLVER_COLS,'("          Solver columns         :",5(X,I13))','(33X,5(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)% &
                    & COUPLING_COEFFICIENTS,'("          Coupling coefficients  :",5(X,E13.6))','(33X,5(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDDO !column_idx
            ENDDO !equations_matrix_idx
          ELSE
            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Linear equations matrix columns to solver matrix columns:", &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of liner equations matrices = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*999)
              DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                  & equations_matrix_idx)%PTR
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations matrix index : ",equations_matrix_idx, &
                  & ERR,ERROR,*999)
                equations_matrix=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX_NUMBER
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix number = ",equations_matrix,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",EQUATIONS_TO_SOLVER_MAP% &
                  & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
                DO column_idx=1,LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix)%NUMBER_OF_COLUMNS
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix column : ",column_idx,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver columns mapped to = ", &
                    & EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS,ERR,ERROR,*999)
                  IF(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS>0) THEN
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                      & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)% &
                      & SOLVER_COLS,'("          Solver columns         :",5(X,I13))','(33X,5(X,I13))',ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                      & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP(column_idx)% &
                      & COUPLING_COEFFICIENTS,'("          Coupling coefficients  :",5(X,E13.6))','(33X,5(X,E13.6))',ERR,ERROR,*999)
                  ENDIF
                ENDDO !column_idx
              ENDDO !equations_matrix_idx
            ENDIF
            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian equations matrix columns to solver matrix columns:", &
                & ERR,ERROR,*999)
              JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",JACOBIAN_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
              DO column_idx=1,NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix column : ",column_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver columns mapped to = ", &
                  & JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS,ERR,ERROR,*999)
                IF(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS>0) THEN
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(column_idx)% &
                    & SOLVER_COLS,'("          Solver columns         :",5(X,I13))','(33X,5(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(column_idx)% &
                    & COUPLING_COEFFICIENTS,'("          Coupling coefficients  :",5(X,E13.6))','(33X,5(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDDO !column_idx
            ENDIF
          ENDIF
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Variable dofs to solver matrix dofs:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",SOLVER_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
            & NUMBER_OF_VARIABLES,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_VARIABLES,5,5,SOLVER_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
            & VARIABLE_TYPES,'("      Variable types :",5(X,I13))','(21X,5(X,I13))',ERR,ERROR,*999)
          DO variable_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_VARIABLES
            DEPENDENT_VARIABLE=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLES(variable_idx)%PTR
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable index : ",variable_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of variable dofs = ",DEPENDENT_VARIABLE% &
              & NUMBER_OF_DOFS,ERR,ERROR,*999)
            DO local_dof=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Variable dof : ",local_dof,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Solver column number = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COLUMN_NUMBERS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Coupling coefficient = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COUPLING_COEFFICIENTS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Additive constant    = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%ADDITIVE_CONSTANTS(local_dof),ERR,ERROR,*999)              
            ENDDO !local_dof
          ENDDO !variable_idx
        ENDDO !solver_matrix_idx
        IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic equations matrix indexing:",ERR,ERROR,*999)
          DO equations_matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations matrix : ",equations_matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of solver matrices = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)% &
              & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
            DO solver_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES
              EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                & solver_matrix_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix index : ",solver_matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix number = ",EQUATIONS_TO_SOLVER_MAP% &
                & EQUATIONS_MATRIX_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",EQUATIONS_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)            
            ENDDO !solver_matrix_idx
          ENDDO !equations_matrix_idx
        ELSE
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Linear equations matrix indexing:",ERR,ERROR,*999)
            DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations matrix : ",equations_matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of solver matrices = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)% &
                & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
              DO solver_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES
                EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                  & solver_matrix_idx)%PTR
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix index : ",solver_matrix_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix number = ",EQUATIONS_TO_SOLVER_MAP% &
                  & EQUATIONS_MATRIX_NUMBER,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",EQUATIONS_TO_SOLVER_MAP% &
                  & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)            
              ENDDO !solver_matrix_idx
            ENDDO !equations_matrix_idx
          ENDIF
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Jacobian matrix indexing:",ERR,ERROR,*999)
            JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM%JACOBIAN_TO_SOLVER_MATRIX_MAP
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",JACOBIAN_TO_SOLVER_MAP% &
              & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDDO ! equations_set_idx       
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_CALCULATE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_CALCULATE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_CALCULATE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver mapping
  SUBROUTINE SOLVER_MAPPING_CREATE_FINISH(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MAPPING_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solver mapping has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
          CALL SOLVER_MAPPING_CALCULATE(SOLVER_MAPPING,ERR,ERROR,*999)
          CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLVER_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
          SOLVER_MAPPING%SOLVER_MAPPING_FINISHED=.TRUE.            
        ELSE
          CALL FLAG_ERROR("Solver mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLVER_MAPPING_CREATE_FINISH")
    RETURN
999 CALL SOLVER_MAPPING_FINALISE(SOLVER_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MAPPING_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver mapping for a problem solver
  SUBROUTINE SOLVER_MAPPING_CREATE_START(SOLVER_EQUATIONS,SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to create the solver mapping on.
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<On return, a pointer to the solver mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solver equations has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          CALL FLAG_ERROR("Solver mapping is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(SOLVER_MAPPING)
          CALL SOLVER_MAPPING_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLVER_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_CREATE_START",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_CREATE_START")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises a solver mapping create values cache and deallocates all memory
  SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solver mapping create values cache 
  SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Solver mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache.",ERR,ERROR,*999)
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache dynamic variable type.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
          & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache matrix variable types.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache residual variable type.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache RHS variable type.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache source variable type.", &
          & ERR,ERROR,*999)
        SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE=0
        SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=0
        SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE=0
        SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=0
        SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLVER_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Destroy a solver mapping.
  SUBROUTINE SOLVER_MAPPING_DESTROY(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      CALL SOLVER_MAPPING_FINALISE(SOLVER_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_DESTROY")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_DESTROY",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a equations column to solver columns map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(EQUATIONS_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE) :: EQUATIONS_COL_SOLVER_COLS_MAP !<The equations col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_COL_SOLVER_COLS_MAP%SOLVER_COLS)) &
      & DEALLOCATE(EQUATIONS_COL_SOLVER_COLS_MAP%SOLVER_COLS)
    IF(ALLOCATED(EQUATIONS_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(EQUATIONS_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations column to solver columns map
  SUBROUTINE SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE(EQUATIONS_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE) :: EQUATIONS_COL_SOLVER_COLS_MAP !<The equations column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_COL_SOLVER_COLS_MAP%NUMBER_OF_SOLVER_COLS=0
    
    CALL EXITS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the mapping of global variables to a solver matrix for the solver mapping
  SUBROUTINE SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET(SOLVER_MAPPING,SOLVER_MATRIX,EQUATIONS_SET_INDEX, &
    & VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: SOLVER_MATRIX !<The solver matrix number to set the equations variables for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_INDEX !<The equations set index in the solver mapping to specify the variable types for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPES(:) !<The variable types to map to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solver mappings has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(SOLVER_MATRIX>=1.AND.SOLVER_MATRIX<=SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
            IF(EQUATIONS_SET_INDEX>=1.AND.EQUATIONS_SET_INDEX<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS) THEN
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(EQUATIONS_SET_INDEX)%PTR
              IF(ASSOCIATED(EQUATIONS_SET)) THEN                
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                  IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                    LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                      IF(SIZE(VARIABLE_TYPES,1)>=1.AND.SIZE(VARIABLE_TYPES,1)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                        DO variable_idx=1,SIZE(VARIABLE_TYPES,1)
!!TODO: CHECK THAT THE VARIABLE TYPE IS NOT REPEATED
                          IF(VARIABLE_TYPES(variable_idx)<1.OR. &
                            & VARIABLE_TYPES(variable_idx)>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                            LOCAL_ERROR="The variable type of "// &
                              & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))// &
                              & " at position "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                              & " in the array is invalid. The number must be >=1 and <= "// &
                              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                          IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(VARIABLE_TYPES(variable_idx))% &
                            & NUMBER_OF_EQUATIONS_MATRICES==0) THEN
                            LOCAL_ERROR="The variable type of "// &
                              & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))// &
                              & " at position "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                              & " in the array is invalid. That variable type is not mapped to any equations matrices"
                          ENDIF
                        ENDDO !variable_idx
                        SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,EQUATIONS_SET_INDEX,SOLVER_MATRIX)= &
                          & SIZE(VARIABLE_TYPES,1)
                        SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:SIZE(VARIABLE_TYPES,1),EQUATIONS_SET_INDEX, &
                          & SOLVER_MATRIX)=VARIABLE_TYPES
                      ELSE
                        LOCAL_ERROR="The supplied size of variable types array of "// &
                          & TRIM(NUMBER_TO_VSTRING(SIZE(VARIABLE_TYPES,1),"*",ERR,ERROR))// &
                          & " is invalid. The size must be between 1 and "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The equations set index of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_INDEX,"*",ERR,ERROR))// &
                & " is invalid. The number must be >= 1 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The solver matrix number of "//TRIM(NUMBER_TO_VSTRING(SOLVER_MATRIX,"*",ERR,ERROR))// &
              & " is invalid. The number must be >= 1 and <= "// &
              & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises a equations row to solver rows map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE(EQUATIONS_ROW_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_ROW_TO_SOLVER_ROWS_MAP_TYPE) :: EQUATIONS_ROW_SOLVER_ROWS_MAP !<The equations row to solver rows map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_ROW_SOLVER_ROWS_MAP%SOLVER_ROWS)) &
      & DEALLOCATE(EQUATIONS_ROW_SOLVER_ROWS_MAP%SOLVER_ROWS)
    IF(ALLOCATED(EQUATIONS_ROW_SOLVER_ROWS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(EQUATIONS_ROW_SOLVER_ROWS_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations row to solver rows map
  SUBROUTINE SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE(EQUATIONS_ROW_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_ROW_TO_SOLVER_ROWS_MAP_TYPE) :: EQUATIONS_ROW_SOLVER_ROWS_MAP !<The equations row to solver rows map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_ROW_SOLVER_ROWS_MAP%NUMBER_OF_SOLVER_ROWS=0
    
    CALL EXITS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds an equations set to a solver mapping
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_ADD(SOLVER_MAPPING,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to add the equations set to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_INDEX !<On exit, the index of the equations set in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,matrix_idx,variable_idx,variable_type
    INTEGER(INTG), ALLOCATABLE :: OLD_DYNAMIC_VARIABLE_TYPE(:),OLD_MATRIX_VARIABLE_TYPES(:,:,:),OLD_RHS_VARIABLE_TYPE(:), &
      & OLD_RESIDUAL_VARIABLE_TYPE(:),OLD_SOURCE_VARIABLE_TYPE(:)
    LOGICAL :: MATRIX_DONE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_SET_PTR_TYPE), ALLOCATABLE :: OLD_EQUATIONS_SETS(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_SET_ADD",ERR,ERROR,*999)

    EQUATIONS_SET_INDEX=0
    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solver mapping has been finished.",ERR,ERROR,*999)
      ELSE
        SOLVER_EQUATIONS=>SOLVER_MAPPING%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          IF(ASSOCIATED(EQUATIONS_SET)) THEN
            IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN                
                  IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
                    IF(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS>0) THEN
                      ALLOCATE(OLD_DYNAMIC_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old dynamic variable type.",ERR,ERROR,*999)
                      ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES,SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types.",ERR,ERROR,*999)
                      ALLOCATE(OLD_RESIDUAL_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old residual variable type.",ERR,ERROR,*999)
                      ALLOCATE(OLD_RHS_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old RHS variable type.",ERR,ERROR,*999)
                      ALLOCATE(OLD_SOURCE_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old source variable type.",ERR,ERROR,*999)
                      OLD_DYNAMIC_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE
                      OLD_MATRIX_VARIABLE_TYPES=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
                      OLD_RESIDUAL_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE
                      OLD_RHS_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE
                      OLD_SOURCE_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE
                      ALLOCATE(OLD_EQUATIONS_SETS(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old equations sets.",ERR,ERROR,*999)
                      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        OLD_EQUATIONS_SETS(equations_set_idx)%PTR=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      ENDDO !equations_set_idx
                      DEALLOCATE(SOLVER_MAPPING%EQUATIONS_SETS)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate residual variable type.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                        & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate residual variable type.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate RHS variable type.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate source variable type.",ERR,ERROR,*999)
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=OLD_DYNAMIC_VARIABLE_TYPE
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS,:)=OLD_MATRIX_VARIABLE_TYPES
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=OLD_RESIDUAL_VARIABLE_TYPE
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=OLD_RHS_VARIABLE_TYPE
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=OLD_SOURCE_VARIABLE_TYPE
                    ELSE IF(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS==0) THEN
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)
                      IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic variable type.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                        & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate residual variable type.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate RHS variable type.",ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate source variable type.",ERR,ERROR,*999)
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS,:)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(1:SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS)=0
                    ELSE
                      CALL FLAG_ERROR("The number of equations sets is < 0.",ERR,ERROR,*999)
                    ENDIF
                    SELECT CASE(SOLVER_EQUATIONS%TIME_DEPENDENCE)
                    CASE(SOLVER_EQUATIONS_STATIC,SOLVER_EQUATIONS_QUASISTATIC)
                      SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
                      CASE(SOLVER_EQUATIONS_LINEAR)
                        IF(ASSOCIATED(EQUATIONS_MAPPING%LINEAR_MAPPING)) THEN
                          !Linear matrices to map. 
                          !Map the first matrix variable found in the equations set to the first solver matrix, the second
                          !variable found to the second, etc.
                          variable_type=1
                          DO matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                            MATRIX_DONE=.FALSE.
                            DO WHILE(variable_type<=FIELD_NUMBER_OF_VARIABLE_TYPES.AND..NOT.MATRIX_DONE)
                              IF(EQUATIONS_MAPPING%LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & NUMBER_OF_EQUATIONS_MATRICES>0) THEN                  
                                SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0, &
                                  & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,matrix_idx)=1
                                SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1, &
                                  & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,matrix_idx)=variable_type
                                MATRIX_DONE=.TRUE.
                              ELSE
                                variable_type=variable_type+1
                              ENDIF
                            ENDDO
                            IF(.NOT.MATRIX_DONE) THEN
                              !Error - could not find any more variables to map to this solver matrix
                              LOCAL_ERROR="Could not find any unmapped variables for solver matrix "// &
                                & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//"."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !matrix_idx
                          !Check if there are still unmapped matrix variables.
                          DO variable_idx=variable_type+1,FIELD_NUMBER_OF_VARIABLE_TYPES
                            IF(EQUATIONS_MAPPING%LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_idx)% &
                              & NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                              LOCAL_ERROR="Variable type "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                                & " is mapped to a linear matrix but has not been mapped to any solver matrices."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !variable_idx
                        ELSE
                          CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
                        ENDIF
                      CASE(SOLVER_EQUATIONS_NONLINEAR)
                        IF(ASSOCIATED(EQUATIONS_MAPPING%NONLINEAR_MAPPING)) THEN
                          !Map the residual variable
                          SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                            & EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                          IF(ASSOCIATED(EQUATIONS_MAPPING%LINEAR_MAPPING)) THEN
                            !If there are linear matrices operating on the residual variable then map them to the
                            !solver matrix (Jacobian)
                            variable_type=EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                            IF(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES==1) THEN
                              IF(EQUATIONS_MAPPING%LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & NUMBER_OF_EQUATIONS_MATRICES>0) THEN                  
                                SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0, &
                                  & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,1)=1
                                SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1, &
                                  & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,1)=variable_type
                              ENDIF
                            ELSE
                              LOCAL_ERROR="Invalid number of solve matrices. For nonlinear solver equations there should "// &
                                & "be 1 solver matrix and there are "// &
                                & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))// &
                                & " solver matrices."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                        ENDIF
                      CASE DEFAULT
                        LOCAL_ERROR="The solver equations linearity type of "// &
                          & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
                     SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
                      CASE(SOLVER_EQUATIONS_LINEAR)
                        IF(ASSOCIATED(EQUATIONS_MAPPING%DYNAMIC_MAPPING)) THEN
                          SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                            & EQUATIONS_MAPPING%DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE
                        ELSE
                          CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                        ENDIF
                     CASE(SOLVER_EQUATIONS_NONLINEAR)
! SEBK 16/09/09 NOT SURE ABOUT SOLVER MAPPING HERE
!|
                        IF(ASSOCIATED(EQUATIONS_MAPPING%DYNAMIC_MAPPING)) THEN
                          SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                            & EQUATIONS_MAPPING%DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE
! new.... need to double check
                          SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                            & EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                        ELSE
                          CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                        ENDIF
!|
! SEBK 16/09/09 NOT SURE ABOUT SOLVER MAPPING HERE
                      CASE DEFAULT
                        LOCAL_ERROR="The solver equations linearity type of "// &
                          & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    CASE DEFAULT
                      LOCAL_ERROR="The solver equations time dependence type of "// &
                        & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                    IF(ASSOCIATED(EQUATIONS_MAPPING%RHS_MAPPING)) THEN
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                        & EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE
                    ENDIF
                    IF(ASSOCIATED(EQUATIONS_MAPPING%SOURCE_MAPPING)) THEN
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                        & EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_VARIABLE_TYPE
                    ENDIF
                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SETS(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations sets.",ERR,ERROR,*999)
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR=>OLD_EQUATIONS_SETS(equations_set_idx)%PTR
                    ENDDO !equations_set_idx
                    SOLVER_MAPPING%EQUATIONS_SETS(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)%PTR=>EQUATIONS_SET
                    SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1
                    EQUATIONS_SET_INDEX=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    IF(ALLOCATED(OLD_DYNAMIC_VARIABLE_TYPE)) DEALLOCATE(OLD_DYNAMIC_VARIABLE_TYPE)
                    IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
                    IF(ALLOCATED(OLD_RESIDUAL_VARIABLE_TYPE)) DEALLOCATE(OLD_RESIDUAL_VARIABLE_TYPE)
                    IF(ALLOCATED(OLD_RHS_VARIABLE_TYPE)) DEALLOCATE(OLD_RHS_VARIABLE_TYPE)
                    IF(ALLOCATED(OLD_SOURCE_VARIABLE_TYPE)) DEALLOCATE(OLD_SOURCE_VARIABLE_TYPE)
                    IF(ALLOCATED(OLD_EQUATIONS_SETS)) DEALLOCATE(OLD_EQUATIONS_SETS)
                  ELSE
                    CALL FLAG_ERROR("Solvers mapping create values cache is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver mapping solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_SET_ADD")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
    IF(ALLOCATED(OLD_RESIDUAL_VARIABLE_TYPE)) DEALLOCATE(OLD_RESIDUAL_VARIABLE_TYPE)
    IF(ALLOCATED(OLD_RHS_VARIABLE_TYPE)) DEALLOCATE(OLD_RHS_VARIABLE_TYPE)
    IF(ALLOCATED(OLD_SOURCE_VARIABLE_TYPE)) DEALLOCATE(OLD_SOURCE_VARIABLE_TYPE)
    IF(ALLOCATED(OLD_EQUATIONS_SETS)) DEALLOCATE(OLD_EQUATIONS_SETS)
    CALL ERRORS("SOLVER_MAPPING_EQUATIONS_SET_ADD",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_SET_ADD")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_ADD

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_matrix_idx,row_idx,solver_matrix_idx
    
    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM)) THEN
      DO solver_matrix_idx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM,1)
        CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx),ERR,ERROR,*999)
      ENDDO !solver_matrix_idx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM)
    ENDIF
    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM)) THEN
      DO equations_matrix_idx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM,1)
        CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx),ERR,ERROR,*999)
      ENDDO !equations_matrix_idx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM)
    ENDIF
    CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,ERR,ERROR,*999)
    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS)) THEN
      DO row_idx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS,1)
        CALL SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx),ERR,ERROR,*999)
      ENDDO !row_idx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver map.
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_SET_INDEX=0
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS)
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)
        
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE(EQUATIONS_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx
    
    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
      IF(ALLOCATED(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP)) THEN
        DO column_idx=1,SIZE(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP,1)
          CALL SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(EQUATIONS_TO_SOLVER_MAP% &
            & EQUATIONS_COL_SOLVER_COLS_MAP(column_idx),ERR,ERROR,*999)
        ENDDO !column_idx
        DEALLOCATE(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP)
      ENDIF
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver maps
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE(EQUATIONS_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP !<The equations to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
      EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX_TYPE=0
      EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX_NUMBER=0
      EQUATIONS_TO_SOLVER_MAP%SOLVER_MATRIX_NUMBER=0
      NULLIFY(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX)
      NULLIFY(EQUATIONS_TO_SOLVER_MAP%SOLVER_MATRIX)
    ELSE
      CALL FLAG_ERROR("Equations to solver map is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map em and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM !<The equations set to solver matrix maps em to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_TO_SOLVER_MATRIX_MAPS)) THEN
      DO matrix_idx=1,SIZE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_TO_SOLVER_MATRIX_MAPS,1)
        CALL SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS(matrix_idx)%PTR,ERR,ERROR,*999)        
      ENDDO !variable_idx
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_TO_SOLVER_MATRIX_MAPS)
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps em.
  SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM !<The equations to solver matrix maps em to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_MATRIX_NUMBER=0
    EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%NUMBER_OF_SOLVER_MATRICES=0
        
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map jm and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM !<The equations set to solver matrix maps jm to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)) THEN
      CALL SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM%JACOBIAN_TO_SOLVER_MATRIX_MAP, &
        & ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps jm.
  SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to initialise the euqations_to_solver_matrix_map_jm for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)) THEN
      CALL FLAG_ERROR("Equations to solver matrix maps jm is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps jm.",ERR,ERROR,*999)
      NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM%JACOBIAN_TO_SOLVER_MATRIX_MAP)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE")
    RETURN
999 CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map sm and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM !<The equations set to solver matrix maps sm to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,variable_idx
    
    CALL ENTERS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TYPES)) DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TYPES)
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLES)) DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLES)
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TO_SOLVER_COL_MAPS)) THEN
      DO variable_idx=1,SIZE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TO_SOLVER_COL_MAPS,1)
        CALL SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM% &
          & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx),ERR,ERROR,*999)        
      ENDDO !variable_idx
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TO_SOLVER_COL_MAPS)
    ENDIF
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS)) THEN
      DO matrix_idx=1,SIZE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS,1)
        CALL SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM% &
          & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(matrix_idx)%PTR,ERR,ERROR,*999)        
      ENDDO !variable_idx
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS)
    ENDIF
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS)) THEN
      DO matrix_idx=1,SIZE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS,1)
        CALL SOLVER_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM% &
          & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(matrix_idx)%PTR,ERR,ERROR,*999)        
      ENDDO !variable_idx
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS)
    ENDIF
    CALL SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%JACOBIAN_TO_SOLVER_MATRIX_MAP, &
      & ERR,ERROR,*999)
    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps sm.
  SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM !<The equations to solver matrix maps sm to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%SOLVER_MATRIX_NUMBER=0
    EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_VARIABLES=0
    EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
    EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
    NULLIFY(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%JACOBIAN_TO_SOLVER_MATRIX_MAP)
        
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_FINALISE(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,row_idx,solver_matrix_idx

    CALL ENTERS("SOLVER_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ALLOCATED(SOLVER_MAPPING%EQUATIONS_SETS)) DEALLOCATE(SOLVER_MAPPING%EQUATIONS_SETS)        
      IF(ALLOCATED(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)) THEN
        DO equations_set_idx=1,SIZE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP,1)
          CALL SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
            & equations_set_idx),ERR,ERROR,*999)
        ENDDO !equations_set_idx
        DEALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP)) THEN
        DO solver_matrix_idx=1,SIZE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP,1)
          CALL SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
            & solver_matrix_idx),ERR,ERROR,*999)
        ENDDO !solver_matrix_idx
        DEALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP)
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS)) THEN
        DO row_idx=1,SIZE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS,1)
          CALL SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS( &
            & row_idx),ERR,ERROR,*999)
        ENDDO !row_idx
        DEALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS)
      ENDIF
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(SOLVER_MAPPING%ROW_DOFS_MAPPING,ERR,ERROR,*999)
      CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLVER_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(SOLVER_MAPPING)
    ENDIF
       
    CALL EXITS("SOLVER_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to initialise the solver mapping on.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
        CALL FLAG_ERROR("Solver equations solver mapping is already associated",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER_EQUATIONS%SOLVER_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver equations solver mapping",ERR,ERROR,*999)
        SOLVER_EQUATIONS%SOLVER_MAPPING%SOLVER_EQUATIONS=>SOLVER_EQUATIONS
        SOLVER_EQUATIONS%SOLVER_MAPPING%SOLVER_MAPPING_FINISHED=.FALSE.
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES=1
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_ROWS=0
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_GLOBAL_ROWS=0
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS=0
        NULLIFY(SOLVER_EQUATIONS%SOLVER_MAPPING%ROW_DOFS_MAPPING)
        NULLIFY(SOLVER_EQUATIONS%SOLVER_MAPPING%CREATE_VALUES_CACHE)
        CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE(SOLVER_EQUATIONS%SOLVER_MAPPING,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_INITIALISE")
    RETURN
999 CALL SOLVER_MAPPING_FINALISE(SOLVER_EQUATIONS%SOLVER_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a Jacobian column to solver columns map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE(JACOBIAN_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE) :: JACOBIAN_COL_SOLVER_COLS_MAP !<The Jacobian col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(JACOBIAN_COL_SOLVER_COLS_MAP%SOLVER_COLS)) &
      & DEALLOCATE(JACOBIAN_COL_SOLVER_COLS_MAP%SOLVER_COLS)
    IF(ALLOCATED(JACOBIAN_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(JACOBIAN_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an Jacobian column to solver columns map
  SUBROUTINE SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE(JACOBIAN_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE) :: JACOBIAN_COL_SOLVER_COLS_MAP !<The Jacobian column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR,*999)

    JACOBIAN_COL_SOLVER_COLS_MAP%NUMBER_OF_SOLVER_COLS=0
    
    CALL EXITS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(JACOBIAN_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP !<The jacobian to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx
    
    CALL ENTERS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
      IF(ALLOCATED(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP)) THEN
        DO column_idx=1,SIZE(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP,1)
          CALL SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE(JACOBIAN_TO_SOLVER_MAP% &
            & JACOBIAN_COL_SOLVER_COLS_MAP(column_idx),ERR,ERROR,*999)
        ENDDO !column_idx
        DEALLOCATE(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP)
      ENDIF
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a Jacobian to solver maps
  SUBROUTINE SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE(JACOBIAN_TO_SOLVER_MATRIX_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MATRIX_MAP !<The Jacobian to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MATRIX_MAP)) THEN
      JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX_NUMBER=0
      NULLIFY(JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_MATRIX)
      NULLIFY(JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX)
    ELSE
      CALL FLAG_ERROR("Jacobian to solver matrix map is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE")
    RETURN
999 CALL SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(JACOBIAN_TO_SOLVER_MATRIX_MAP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to dynamic equations map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_FINALISE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP !<The solver column to dynamic equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP%EQUATIONS_MATRIX_NUMBERS)) &
      & DEALLOCATE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP%EQUATIONS_MATRIX_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP%EQUATIONS_COL_NUMBERS)) &
      & DEALLOCATE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP%EQUATIONS_COL_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP%COUPLING_COEFFICIENTS)
       
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to dynamic equations mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_INITIALISE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP !<The solver column to dynamic equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
     
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to static equations map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_STATIC_EQUATIONS_MAP !<The solver column to static equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_STATIC_EQUATIONS_MAP%EQUATIONS_MATRIX_NUMBERS)) &
      & DEALLOCATE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP%EQUATIONS_MATRIX_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_STATIC_EQUATIONS_MAP%EQUATIONS_COL_NUMBERS)) &
      & DEALLOCATE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP%EQUATIONS_COL_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_STATIC_EQUATIONS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP%COUPLING_COEFFICIENTS)
       
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to static equations mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_INITIALISE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_STATIC_EQUATIONS_MAP !<The solver column to static equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_STATIC_EQUATIONS_MAP%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
    SOLVER_COL_TO_STATIC_EQUATIONS_MAP%JACOBIAN_COL_NUMBER=0    
    SOLVER_COL_TO_STATIC_EQUATIONS_MAP%JACOBIAN_COUPLING_COEFFICIENT=0.0_DP
    
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations set map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SET_MAP !<The solver column to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: col_idx
    
    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS)) THEN
      DO col_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS,1)
        CALL SOLVER_MAPPING_SOLVER_COL_TO_D_EQUATIONS_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SET_MAP% &
          & SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(col_idx),ERR,ERROR,*999)
      ENDDO !col_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS)
    ENDIF
    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS)) THEN
      DO col_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS,1)
        CALL SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SET_MAP% &
          & SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(col_idx),ERR,ERROR,*999)
      ENDDO !col_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS)
    ENDIF
       
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations set mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE(SOLVER_COL_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SET_MAP !<The solver column to equations set map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR,*999)

    NULLIFY(SOLVER_COL_TO_EQUATIONS_SET_MAP%EQUATIONS)
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations sets map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SETS_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SETS_MAP !<The solver column to equations sets map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: col_idx,equations_set_idx
    
    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_S_EQUATS_SETS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_COL_TO_EQUATIONS_SET_MAPS)) THEN
      DO equations_set_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_COL_TO_EQUATIONS_SET_MAPS,1)
        CALL SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP% &
          SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx),ERR,ERROR,*999)
      ENDDO !equations_set_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_COL_TO_EQUATIONS_SET_MAPS)
    ENDIF
    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS)) THEN
      DO col_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS,1)
        CALL SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS( &
          & col_idx),ERR,ERROR,*999)
      ENDDO !col_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS)
    ENDIF
    CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%COLUMN_DOFS_MAPPING,ERR,ERROR,*999)
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_S_EQUATS_SETS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_S_EQUATS_SETS_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations sets mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SETS_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SETS_MAP !<The solver column to equations sets map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_MATRIX_NUMBER=0
    NULLIFY(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_MATRIX)
    NULLIFY(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_MAPPING)
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%NUMBER_OF_COLUMNS=0
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%NUMBER_OF_DOFS=0
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%TOTAL_NUMBER_OF_DOFS=0
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%NUMBER_OF_GLOBAL_DOFS=0
    NULLIFY(SOLVER_COL_TO_EQUATIONS_SETS_MAP%COLUMN_DOFS_MAPPING)
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver dof to variable mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE(SOLVER_DOF_TO_VARIABLE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_DOF_TO_VARIABLE_MAP_TYPE) :: SOLVER_DOF_TO_VARIABLE_MAP !<The solver dof to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_SET_INDICES)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_SET_INDICES)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_DOF)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_DOF)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_COEFFICIENT)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_COEFFICIENT)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%ADDITIVE_CONSTANT)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%ADDITIVE_CONSTANT)
    
    CALL EXITS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver dof to variable mapping.
  SUBROUTINE SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE(SOLVER_DOF_TO_VARIABLE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_DOF_TO_VARIABLE_MAP_TYPE) :: SOLVER_DOF_TO_VARIABLE_MAP !<The solver dof to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_DOF_TO_VARIABLE_MAP%NUMBER_OF_EQUATIONS_SETS=0
    
    CALL EXITS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solver matrices for the solver mapping
  SUBROUTINE SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLVER_MAPPING,NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SOLVER_MATRICES !<The number of solver matrices for the solver.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,matrix_idx,MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_VARIABLE_TYPES(:,:,:)
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solver mappings has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
          MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES=1
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            EQUATIONS=>EQUATIONS_SET%EQUATIONS
            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
              IF(LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES) &
              & MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES=LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
            ENDIF
          ENDDO !equations_set_idx
          !Check number of matrices to set is valid
          IF(NUMBER_OF_SOLVER_MATRICES>0.AND.NUMBER_OF_SOLVER_MATRICES<=MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES) THEN
            !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
            IF(NUMBER_OF_SOLVER_MATRICES/=SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
              ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types",ERR,ERROR,*999)
              OLD_MATRIX_VARIABLE_TYPES=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
              DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
              ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types",ERR,ERROR,*999)
              IF(NUMBER_OF_SOLVER_MATRICES>SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
                SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,:,1:SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES)= &
                  & OLD_MATRIX_VARIABLE_TYPES(:,:,1:SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES)
                DO matrix_idx=SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES+1,NUMBER_OF_SOLVER_MATRICES
                  SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,:,matrix_idx)=1
                  SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1,:,matrix_idx)=FIELD_U_VARIABLE_TYPE
                  SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(2:FIELD_NUMBER_OF_VARIABLE_TYPES,:,matrix_idx)=0
                ENDDO !matrix_idx
              ELSE
                SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,:,1:NUMBER_OF_SOLVER_MATRICES)= &
                  & OLD_MATRIX_VARIABLE_TYPES(:,:,1:NUMBER_OF_SOLVER_MATRICES)
              ENDIF
              SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES=NUMBER_OF_SOLVER_MATRICES
              IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified number of solver matrices of "// &
              & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))// &
              & " is invalid. The number must be >= 1 and <= "// &
              & TRIM(NUMBER_TO_VSTRING(MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
    CALL ERRORS("SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises a variable to solver column map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE(SOLVER_ROW_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_ROW_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_ROW_TO_EQUATIONS_SET_MAP !<The solver row to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_SET)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_SET)
    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_ROW_NUMBER)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_ROW_NUMBER)
    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_SET_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_SET_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solver row to equations set map.
  SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE(SOLVER_ROW_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_ROW_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_ROW_TO_EQUATIONS_SET_MAP !<The solver row to equations set map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR,*999)

   SOLVER_ROW_TO_EQUATIONS_SET_MAP%NUMBER_OF_ROWS=0
        
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a variable to solver column map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE(VARIABLE_TO_SOLVER_COL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VARIABLE_TO_SOLVER_COL_MAP_TYPE) :: VARIABLE_TO_SOLVER_COL_MAP !<The variable to solver column map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(VARIABLE_TO_SOLVER_COL_MAP%COLUMN_NUMBERS)) DEALLOCATE(VARIABLE_TO_SOLVER_COL_MAP%COLUMN_NUMBERS)
    IF(ALLOCATED(VARIABLE_TO_SOLVER_COL_MAP%COUPLING_COEFFICIENTS)) DEALLOCATE(VARIABLE_TO_SOLVER_COL_MAP%COUPLING_COEFFICIENTS)
    IF(ALLOCATED(VARIABLE_TO_SOLVER_COL_MAP%ADDITIVE_CONSTANTS)) DEALLOCATE(VARIABLE_TO_SOLVER_COL_MAP%ADDITIVE_CONSTANTS)
        
    CALL EXITS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE

  !
  !================================================================================================================================
  !
  
END MODULE SOLVER_MAPPING_ROUTINES

