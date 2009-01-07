!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all solution mapping routines.
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

!> This module handles all solution mapping routines.
MODULE SOLUTION_MAPPING_ROUTINES

  USE BASE_ROUTINES
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

  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLUTION_MAPPING_CREATE_FINISH,SOLUTION_MAPPING_CREATE_START,SOLUTION_MAPPING_DESTROY, &
    & SOLUTION_MAPPING_EQUATIONS_SET_ADD,SOLUTION_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET, &
    & SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the solution mappings
  SUBROUTINE SOLUTION_MAPPING_CALCULATE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,dof_idx,equations_column,equations_matrix,equations_matrix_idx,equations_set_idx, &
      & global_dof,global_field_dof,global_row,jacobian_column,local_dof,LINEAR_EQUATIONS_MATRIX_OFFSET,local_row, &
      & matrix_number,myrank,myrank_local_dof,NUMBER_OF_COLUMNS,NUMBER_OF_LINEAR_EQUATIONS_MATRICES,NUMBER_OF_GLOBAL_SOLVER_COLS, &
      & LOCAL_SOLVER_DOF_OFFSET,NUMBER_OF_GHOST_SOLVER_DOFS,NUMBER_OF_GLOBAL_SOLVER_ROWS,NUMBER_OF_LOCAL_SOLVER_COLS, &
      & NUMBER_OF_LOCAL_SOLVER_DOFS,NUMBER_OF_LOCAL_SOLVER_ROWS,NUMBER_OF_VARIABLES,rank,rank_idx,row_idx,SOLVER_DOF, &
      & solver_matrix_idx,TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS,variable_idx,variable_type
    LOGICAL :: INCLUDE_ROW,MYRANK_DOF,RANK_DOF
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COL_DOMAIN_MAPPING,COL_DOFS_MAPPING,DEPENDENT_DOFS_MAPPING,ROW_DOMAIN_MAPPING, &
      & ROW_DOFS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION

    CALL ENTERS("SOLUTION_MAPPING_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
        SOLUTION=>SOLUTION_MAPPING%SOLUTION
        IF(ASSOCIATED(SOLUTION)) THEN          
          !
          !--- Row mappings ---
          !
          !Calculate the row mappings.
          !We do not have any couplings defined at the moment there is only a 1-1 mapping.
          myrank=COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER
          NUMBER_OF_GLOBAL_SOLVER_ROWS=0
          NUMBER_OF_LOCAL_SOLVER_ROWS=0                                 
          DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING                
                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                  LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                  NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                  RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                  SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
                  ROW_DOFS_MAPPING=>EQUATIONS_MAPPING%ROW_DOFS_MAPPING
                  IF(ASSOCIATED(ROW_DOFS_MAPPING)) THEN                    
                    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                      DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
                      IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
                        FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                        IF(ASSOCIATED(FIXED_CONDITIONS)) THEN                      
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
                                IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                                  DEPENDENT_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLE
                                  !This is wrong as we only have the mappings for the local rank not the global ranks.
                                  !For now assume 1-1 mapping between rows and dofs.
                                  global_dof=global_row                                  
                                  global_field_dof=DEPENDENT_VARIABLE%GLOBAL_DOF_OFFSET+global_dof
                                  INCLUDE_ROW= &
                                    & FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_field_dof)==EQUATIONS_SET_NOT_FIXED
                                ELSE IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                  !Loop over the variables in the equations set. Don't include the row in the solver matrices if
                                  !all the variable dofs associated with this equations row are fixed.
                                  DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                                    !This is wrong as we only have the mappings for the local rank not the global ranks.
                                    !For now assume 1-1 mapping between rows and dofs.
                                    !
                                    !local_dof=EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(local_row)% &
                                    !  & ROW_TO_DOFS_MAP(equations_matrix_idx)
                                    !variable_type=EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES(equations_matrix_idx)
                                    DEPENDENT_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)% &
                                      & VARIABLE
                                    !global_dof=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_dof)
                                    global_dof=global_row                                    
                                    global_field_dof=DEPENDENT_VARIABLE%GLOBAL_DOF_OFFSET+global_dof
                                    INCLUDE_ROW=INCLUDE_ROW.OR. &
                                      & FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_field_dof)==EQUATIONS_SET_NOT_FIXED
                                  ENDDO !matrix_idx
                                ENDIF
                                IF(INCLUDE_ROW) THEN
                                  NUMBER_OF_GLOBAL_SOLVER_ROWS=NUMBER_OF_GLOBAL_SOLVER_ROWS+1
                                  IF(rank==myrank) &
                                    & NUMBER_OF_LOCAL_SOLVER_ROWS=NUMBER_OF_LOCAL_SOLVER_ROWS+1 !1-1 mapping
                                ENDIF !include row
                              ENDIF !rank dof
                            ENDDO !global_row
                          ENDDO !rank
                        ELSE
                          CALL FLAG_ERROR("Equations set fixed conditions is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Dependent field domain mapping is not associated.",ERR,ERROR,*999)
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
          !Allocate memory for the rows mapping
          !Allocate equations set to solver map
          ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping equations set to solver map.",ERR,ERROR,*999)      
          !Allocate the solver rows to equations set maps
          ALLOCATE(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping solver row to equations set map.",ERR,ERROR,*999)
          !Set the number of rows
          SOLUTION_MAPPING%NUMBER_OF_ROWS=NUMBER_OF_LOCAL_SOLVER_ROWS
          SOLUTION_MAPPING%NUMBER_OF_GLOBAL_ROWS=NUMBER_OF_GLOBAL_SOLVER_ROWS
          !Allocate the solver rows domain mapping
          ALLOCATE(SOLUTION_MAPPING%ROW_DOFS_MAPPING,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping row dofs mapping.",ERR,ERROR,*999)
!!TODO: what is the real number of domains for a solution???
          CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(SOLUTION_MAPPING%ROW_DOFS_MAPPING,COMPUTATIONAL_ENVIRONMENT% &
            & NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)
          ROW_DOMAIN_MAPPING=>SOLUTION_MAPPING%ROW_DOFS_MAPPING
          ALLOCATE(ROW_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row dofs mapping global to local map.",ERR,ERROR,*999)
          ROW_DOMAIN_MAPPING%NUMBER_OF_GLOBAL=NUMBER_OF_GLOBAL_SOLVER_ROWS              
          !Calculate the row mappings
          NUMBER_OF_GLOBAL_SOLVER_ROWS=0
          !Loop over the ranks to  ensure that the lowest ranks have the lowest numbered solver variables
          DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            NUMBER_OF_LOCAL_SOLVER_ROWS=0
            DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
              !Note that pointers have been checked for association above
              EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
              ROW_DOFS_MAPPING=>EQUATIONS_MAPPING%ROW_DOFS_MAPPING
              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
              IF(rank==myrank) THEN
                !Initialise the equations set to solver map
                CALL SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                  & equations_set_idx),ERR,ERROR,*999)
                SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_SET_INDEX=equations_set_idx
                SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%SOLUTION_MAPPING=>SOLUTION_MAPPING
                SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS=>EQUATIONS
                !Allocate the equations set to solver maps for solver matrix (sm) indexing
                ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                  & SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set to solver map equations to solver matrix maps sm.", &
                  & ERR,ERROR,*999)
                IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                  !Allocate the equations set to solver maps for equations matrix (em) indexing
                  ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                    & LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                  IF(ERR/=0) &
                    & CALL FLAG_ERROR("Could not allocate equations set to solver map equations to solver matrix maps em.", &
                    & ERR,ERROR,*999)
                ENDIF
                IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                  !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
                  CALL SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE(SOLUTION_MAPPING% &
                    & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx),ERR,ERROR,*999)
                ENDIF
                !Allocate the equations row to solver rows maps
                ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                  & EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set to solver map equations row to solver rows maps.", &
                  & ERR,ERROR,*999)
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
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    DEPENDENT_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLE
                    !This is wrong as we only have the mappings for the local rank not the global ranks.
                    !For now assume 1-1 mapping between rows and dofs.
                    global_dof=global_row
                    global_field_dof=DEPENDENT_VARIABLE%GLOBAL_DOF_OFFSET+global_dof
                    INCLUDE_ROW=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_field_dof)==EQUATIONS_SET_NOT_FIXED
                  ELSE IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                    DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                      !This is wrong as we only have the mappings for the local rank not the global ranks.
                      !For now assume 1-1 mapping between rows and dofs.
                      !
                      !local_dof=EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(local_row)% &
                      !  & ROW_TO_DOFS_MAP(equations_matrix_idx)
                      !variable_type=EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES(equations_matrix_idx)
                      DEPENDENT_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)%VARIABLE
                      !global_dof=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_dof)
                      global_dof=global_row
                      global_field_dof=DEPENDENT_VARIABLE%GLOBAL_DOF_OFFSET+global_dof
                      INCLUDE_ROW=INCLUDE_ROW.OR. &
                        & FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_field_dof)==EQUATIONS_SET_NOT_FIXED
                    ENDDO !matrix_idx
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
                      CALL SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE(SOLUTION_MAPPING% &
                        & SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS),ERR,ERROR,*999)
                      !Allocate the solver row to equations row mapping arrays
                      ALLOCATE(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_SET(1), &
                        & STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row to equations set map equations set.",ERR,ERROR,*999)
                      ALLOCATE(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)% &
                        & EQUATIONS_ROW_NUMBER(1),STAT=ERR)
                      IF(ERR/=0) &
                        & CALL FLAG_ERROR("Could not allocate row to equations set map equations row number.",ERR,ERROR,*999)
                      ALLOCATE(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)% &
                        & COUPLING_COEFFICIENTS(1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row to equations set map coupling coefficients.", &
                        & ERR,ERROR,*999)
                      !Set the mappings
                      SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%NUMBER_OF_ROWS=1
                      SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_SET(1)= &
                        & equations_set_idx
                      SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_ROW_NUMBER(1)= &
                        & local_row
                      SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(NUMBER_OF_LOCAL_SOLVER_ROWS)%COUPLING_COEFFICIENTS(1)= &
                        & 1.0_DP
                      !Set up the equations row -> solver row mappings
                      !Initialise
                      CALL SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE(SOLUTION_MAPPING% &
                        & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(local_row), &
                        & ERR,ERROR,*999)
                      !Allocate the equations row to solver row mappings arrays
                      ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%SOLVER_ROWS(1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to solver rows maps solver rows.", &
                        & ERR,ERROR,*999)
                      ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to solver rows maps solver rows.", &
                        & ERR,ERROR,*999)
                      !Set the mappings
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%NUMBER_OF_SOLVER_ROWS=1
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%SOLVER_ROWS(1)=NUMBER_OF_LOCAL_SOLVER_ROWS
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                        & local_row)%COUPLING_COEFFICIENTS(1)=1.0_DP
                    ENDIF !rank==my rank
                  ELSE
                    IF(rank==myrank) THEN
                      !Set up the equations row -> solver row mappings
                      !Initialise
                      CALL SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE(SOLUTION_MAPPING% &
                        & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(local_row), &
                        & ERR,ERROR,*999)
                      !Set the mappings
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
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
          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping solver column to equations sets map.",ERR,ERROR,*999)
          !Calculate the column mappings for each solver matrix
          DO solver_matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
            !Calculate the number of solution variables
            NUMBER_OF_GLOBAL_SOLVER_COLS=0
            NUMBER_OF_LOCAL_SOLVER_COLS=0
            TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS=0
            !Initialise solver column to equations sets mapping array
            CALL SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE(SOLUTION_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx),ERR,ERROR,*999)
            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_MATRIX_NUMBER=solver_matrix_idx
            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLUTION_MAPPING=>SOLUTION_MAPPING
            !Allocate the solver col to equations set maps array
            ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver col to equations sets map solver col to equation set maps.", &
              & ERR,ERROR,*999)
            !Loop over the ranks
            DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
              !Loop over the equations sets
              DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                !The pointers below have been checked for association above.
                EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                  NUMBER_OF_VARIABLES=1
                ELSE
                  NUMBER_OF_VARIABLES=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,equations_set_idx, &
                    & solver_matrix_idx)
                ENDIF
                IF(rank==0) THEN
                  !Initialise equations set to solver map (sm)
                  CALL SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE(SOLUTION_MAPPING% &
                    & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx), &
                    & ERR,ERROR,*999)
                  SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                  !Allocate the equations set to solver map variables arrays
                  ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%VARIABLE_TYPES(NUMBER_OF_VARIABLES),STAT=ERR)
                  IF(ERR/=0)  &
                    & CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm variable types.",ERR,ERROR,*999)
                  ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%VARIABLES(NUMBER_OF_VARIABLES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm variables.",ERR,ERROR,*999)
                  ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(NUMBER_OF_VARIABLES),STAT=ERR)
                  IF(ERR/=0)  &
                    & CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm variables to solver col maps.", &
                    & ERR,ERROR,*999)                  
                  !Setup
                  SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%NUMBER_OF_VARIABLES=NUMBER_OF_VARIABLES
                  NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                ENDIF !rank==0
!!TODO: see how slow this is. At the moment we go through number of ranks*number of globals. We could presort the global nys into a list for each rank. This would take additional memory.
                DO variable_idx=1,NUMBER_OF_VARIABLES
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(equations_set_idx)
                  ELSE
                    variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,equations_set_idx, &
                      & solver_matrix_idx)
                  ENDIF
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(COL_DOFS_MAPPING)) THEN
                    IF(rank==0) THEN
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%VARIABLE_TYPES(variable_idx)=variable_type
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%VARIABLES(variable_idx)%PTR=>DEPENDENT_VARIABLE
                      !Allocate the variable to solver col maps arrays
                      ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)% &
                        & COLUMN_NUMBERS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps column numbers.", &
                        & ERR,ERROR,*999)
                      ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)% &
                        & COUPLING_COEFFICIENTS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps coupling coefficients.", &
                        & ERR,ERROR,*999)
                      ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)% &
                        & ADDITIVE_CONSTANTS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps additive constants.", &
                        & ERR,ERROR,*999)
                      !Setup
                      IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                        NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES+LINEAR_MAPPING% &
                          & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
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
                        global_field_dof=DEPENDENT_VARIABLE%GLOBAL_DOF_OFFSET+global_dof
                        IF(FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_field_dof)==EQUATIONS_SET_NOT_FIXED) THEN
                          NUMBER_OF_GLOBAL_SOLVER_COLS=NUMBER_OF_GLOBAL_SOLVER_COLS+1
                          IF(MYRANK_DOF) TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS=TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS+1
                          IF(rank==myrank) &
                            & NUMBER_OF_LOCAL_SOLVER_COLS=NUMBER_OF_LOCAL_SOLVER_COLS+1
                        ENDIF !field dof not fixed
                      ENDIF !rank dof
                    ENDDO !global_dof
                  ELSE
                    CALL FLAG_ERROR("Equations matrix columns degree of freedom mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !variable_idx
              ENDDO !equations_set_idx
            ENDDO !rank
            !Allocate memory for this solver matrix
            !Allocate solver columns to equations sets maps
            ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
              & TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps.",ERR,ERROR,*999)
            !Set the number of columns
            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS=NUMBER_OF_GLOBAL_SOLVER_COLS
            !Set the number of variables
            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_DOFS=NUMBER_OF_LOCAL_SOLVER_COLS
            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS= &
              & TOTAL_NUMBER_OF_LOCAL_SOLVER_COLS
            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_GLOBAL_DOFS=NUMBER_OF_GLOBAL_SOLVER_COLS
            !Allocate the columns domain mapping
            ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%COLUMN_DOFS_MAPPING,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver col to equations sets map column dofs mapping.",ERR,ERROR,*999)
!!TODO: what is the real number of domains for a solution???
            CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
              & COLUMN_DOFS_MAPPING,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)            
            COL_DOMAIN_MAPPING=>SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%COLUMN_DOFS_MAPPING
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
              DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                !The pointers below have been checked for association above.
                EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                  NUMBER_OF_VARIABLES=1
                ELSE
                  NUMBER_OF_VARIABLES=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,equations_set_idx, &
                    & solver_matrix_idx)
                ENDIF
                IF(rank==0) THEN
                  !Allocate memory
                  !Initialise solver columns to equations set map
                  CALL SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE(SOLUTION_MAPPING% &
                    & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                    & equations_set_idx),ERR,ERROR,*999)
                  !Allocate the solver columns to equations set map arrays
                  ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                    & equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS(NUMBER_OF_COLUMNS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver columns to equations map.",ERR,ERROR,*999)
                  !Set the solver column to equations set map
                  SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                    & equations_set_idx)%EQUATIONS=>EQUATIONS
                  !Allocate the equations to solver matrix maps sm equations to solver maps
                  ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                  IF(ERR/=0) &
                    & CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm equations to solver matrix maps.", &
                    & ERR,ERROR,*999)
                  !Set up linear arrays
                  SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    NULLIFY(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR)
                    ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR,STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not equations to solver matrix maps.",ERR,ERROR,*999)
                    CALL SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                      & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                      & equations_matrix_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !equations_matrix_idx
                  !Set up nonlinear arrays
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN                   
                    ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP,STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not Jacobian to solver matrix maps.",ERR,ERROR,*999)
                    CALL SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                      & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP, &
                      & ERR,ERROR,*999)
                  ENDIF
                ENDIF !rank==0
                !Loop over the variables
                LINEAR_EQUATIONS_MATRIX_OFFSET=0
                DO variable_idx=1,NUMBER_OF_VARIABLES
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(equations_set_idx)
                  ELSE
                    variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,equations_set_idx, &
                      & solver_matrix_idx)
                  ENDIF
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                    NUMBER_OF_LINEAR_EQUATIONS_MATRICES=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  ENDIF
                  IF(rank==0) THEN
                    SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%VARIABLES(variable_idx)%PTR=>DEPENDENT_VARIABLE
                    !Allocate linear equations to solver matrix maps equations column to solver columns maps
                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                      DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                        MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                          & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                        SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                          & equations_matrix_idx)%PTR%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                        SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                          & equations_matrix_idx)%PTR%EQUATIONS_MATRIX_NUMBER=MATRIX_NUMBER
                        SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                          & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                          & equations_matrix_idx)%PTR%EQUATIONS_MATRIX=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                          & MATRIX_NUMBER)%EQUATIONS_MATRIX
                        NUMBER_OF_COLUMNS=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                        ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                          & LINEAR_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                          & NUMBER_OF_COLUMNS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations column to solver columns map.",ERR,ERROR,*999)
                      ENDDO !equations_matrix_idx
                    ENDIF
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                      SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_MATRIX=>NONLINEAR_MAPPING% &
                        & JACOBIAN_TO_VAR_MAP%JACOBIAN
                      NUMBER_OF_COLUMNS=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
                      ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                        & JACOBIAN_COL_SOLVER_COLS_MAP(NUMBER_OF_COLUMNS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate jacobian column to solver columns map.",ERR,ERROR,*999)
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
                      global_field_dof=DEPENDENT_VARIABLE%GLOBAL_DOF_OFFSET+global_dof
                      IF(FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_field_dof)==EQUATIONS_SET_NOT_FIXED) THEN
                        !DOF is not fixed so map the variable/equation dof to a new solution dof
                        NUMBER_OF_GLOBAL_SOLVER_COLS=NUMBER_OF_GLOBAL_SOLVER_COLS+1
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
                          !Initialise_sm
                          CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                            & NUMBER_OF_GLOBAL_SOLVER_COLS),ERR,ERROR,*999)
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
                          !Initialise
                          CALL SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_INITIALISE(SOLUTION_MAPPING% &
                            & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                            & equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS(NUMBER_OF_GLOBAL_SOLVER_COLS),ERR,ERROR,*999)
                          !Allocate the solver column to equations column mapping arrays
                          !No coupling yet so the number of columns the solver column is mapped to is just the number of matrices
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                            & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_MATRIX_NUMBERS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrix numbers.",ERR,ERROR,*999)
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                            & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_COL_NUMBERS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations column numbers.",ERR,ERROR,*999)
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                            & NUMBER_OF_GLOBAL_SOLVER_COLS)%COUPLING_COEFFICIENTS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations coupling coefficients.",ERR,ERROR,*999)
                          !Set the mappings
                          SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                            & NUMBER_OF_GLOBAL_SOLVER_COLS)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                            !Loop over the linear equations matrices associated with the variable and set the column maps
                            DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                              !Set the column map
                              MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                              equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                              SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)=MATRIX_NUMBER
                              SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%EQUATIONS_COL_NUMBERS(equations_matrix_idx)=equations_column
                              SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                                & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS)%COUPLING_COEFFICIENTS(equations_matrix_idx)= &
                                & LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%MATRIX_COEFFICIENT
                            ENDDO !equations_matrix_idx
                          ENDIF
                          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                            jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(myrank_local_dof)
                            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%JACOBIAN_COL_NUMBER=jacobian_column
                            SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                              & NUMBER_OF_GLOBAL_SOLVER_COLS)%JACOBIAN_COUPLING_COEFFICIENT= &
                              & NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT
                          ENDIF
                          !Set up the solver dofs -> variable dofs map
                          !Initialise
                          CALL SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE(SOLUTION_MAPPING% &
                            & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF),ERR,ERROR,*999)
                          !Allocate the solver dofs to variable dofs arrays
                          !No coupling so there is only one equations set at the moment
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%EQUATIONS_SET_INDICES(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps equations set indices.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%VARIABLE(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable type.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%VARIABLE_DOF(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable dof.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%VARIABLE_COEFFICIENT(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable coefficient.", &
                            & ERR,ERROR,*999)
                          ALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                            & SOLVER_DOF_TO_VARIABLE_MAPS(SOLVER_DOF)%ADDITIVE_CONSTANT(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps additive constant.", &
                            & ERR,ERROR,*999)
                          !Setup
                          SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%NUMBER_OF_EQUATIONS_SETS=1
                          SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%EQUATIONS_SET_INDICES(1)=equations_set_idx
                          SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%VARIABLE(1)%PTR=>DEPENDENT_VARIABLE
                          SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%VARIABLE_DOF(1)=myrank_local_dof
                          SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%VARIABLE_COEFFICIENT(1)=1.0_DP
                          SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                            & SOLVER_DOF)%ADDITIVE_CONSTANT(1)=0.0_DP
                          !Set up the equations variables -> solver columns mapping
                          !No coupling yet so the mapping is 1-1
                          SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COLUMN_NUMBERS(myrank_local_dof)= &
                            & NUMBER_OF_GLOBAL_SOLVER_COLS
                          SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COUPLING_COEFFICIENTS( &
                            & myrank_local_dof)=1.0_DP
                          SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%ADDITIVE_CONSTANTS( &
                            & myrank_local_dof)=0.0_DP
                          !Set up the equations columns -> solver columns mapping
                          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                            DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                              MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                              equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                              !Allocate the equation to solver map column items.
                              !No coupling yet so the mapping is 1-1
                              ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                                & LINEAR_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                & equations_column)%SOLVER_COLS(1),STAT=ERR)
                              IF(ERR/=0) CALL  &
                                & FLAG_ERROR("Could not allocate equations column to solver columns map solver colums.", &
                                & ERR,ERROR,*999)
                              ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                                & LINEAR_EQUATIONS_MATRIX_OFFSET+equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                & equations_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                              IF(ERR/=0) CALL &
                                & FLAG_ERROR("Could not allocate equations column to solver columns map coupling coefficients.",&
                                & ERR,ERROR,*999)
                              SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=1
                              SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)%SOLVER_COLS(1)= &
                                & NUMBER_OF_GLOBAL_SOLVER_COLS
                              SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)% &
                                & COUPLING_COEFFICIENTS(1)=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)% &
                                & MATRIX_COEFFICIENT                              
                            ENDDO !equations_matrix_idx
                          ENDIF
                          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                            jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(myrank_local_dof)
                            !Allocate the Jacobian to solver map column items.
                            !No coupling yet so the mapping is 1-1
                            ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                              & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)%SOLVER_COLS(1),STAT=ERR)
                            IF(ERR/=0) CALL  &
                              & FLAG_ERROR("Could not allocate Jacobian column to solver columns map solver colums.", &
                              & ERR,ERROR,*999)
                            ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                              & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                            IF(ERR/=0) CALL &
                              & FLAG_ERROR("Could not allocate Jacobain column to solver columns map coupling coefficients.",&
                              & ERR,ERROR,*999)
                            SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                              & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)% &
                              & NUMBER_OF_SOLVER_COLS=1
                            SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                              & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)% &
                              & SOLVER_COLS(1)=NUMBER_OF_GLOBAL_SOLVER_COLS
                            SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                              & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column)% &
                              & COUPLING_COEFFICIENTS(1)=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT
                          ENDIF
                        ENDIF
                      ELSE
                        IF(MYRANK_DOF) THEN
                          !Set up the equations variables -> solver columns mapping
                          !No coupling yet so the mapping is 1-1
                          SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COLUMN_NUMBERS(myrank_local_dof)=0
                          SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%COUPLING_COEFFICIENTS( &
                            & myrank_local_dof)=0.0_DP
                          SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                            & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_type)%ADDITIVE_CONSTANTS( &
                            & myrank_local_dof)=0.0_DP
                          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                            !Set up the equations columns -> solver columns mapping
                            DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                              MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                              equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(myrank_local_dof)
                              !No coupling yet so the mapping is 1-1
                              SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(LINEAR_EQUATIONS_MATRIX_OFFSET+ &
                                & equations_matrix_idx)%PTR%EQUATIONS_COL_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=0
                            ENDDO !equations_matrix_idx
                          ENDIF
                          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                            jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(myrank_local_dof)
                            !No coupling yet so the mapping is 1-1
                            SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                              & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                              & jacobian_column)%NUMBER_OF_SOLVER_COLS=0
                          ENDIF
                        ENDIF
                      ENDIF !field dof is fixed
                    ENDIF !rank dof
                  ENDDO !global_dof
                  LINEAR_EQUATIONS_MATRIX_OFFSET=LINEAR_EQUATIONS_MATRIX_OFFSET+NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                ENDDO !variable_idx
              ENDDO !equations_set_idx
            ENDDO !rank
            CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(COL_DOMAIN_MAPPING,ERR,ERROR,*999)
          ENDDO !solver_matrix_idx
        ELSE
          CALL FLAG_ERROR("The solution mapping solution is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution mapping create values cache is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Solution mappings:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of equations sets = ",SOLUTION_MAPPING% &
        & NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
      DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equations_set_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Region user number = ",EQUATIONS_SET%REGION%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set user number = ",EQUATIONS_SET%USER_NUMBER, &
          & ERR,ERROR,*999)                
      ENDDO !equations_set_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Row mappings:",ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",SOLUTION_MAPPING%NUMBER_OF_ROWS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global rows = ",SOLUTION_MAPPING%NUMBER_OF_GLOBAL_ROWS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver rows to equations sets rows mappings:",ERR,ERROR,*999)      
      DO row_idx=1,SOLUTION_MAPPING%NUMBER_OF_ROWS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver row : ",row_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of rows mapped to = ",SOLUTION_MAPPING% &
          & SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%NUMBER_OF_ROWS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)% &
          & NUMBER_OF_ROWS,5,5,SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%EQUATIONS_SET, &
          & '("      Equations sets indices  :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)% &
          & NUMBER_OF_ROWS,5,5,SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%EQUATIONS_ROW_NUMBER, &
          & '("      Equations row numbers   :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)% &
          & NUMBER_OF_ROWS,5,5,SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS(row_idx)%COUPLING_COEFFICIENTS, &
          & '("      Coupling coefficients   :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999) 
      ENDDO !row_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equations sets rows to solver rows mappings:",ERR,ERROR,*999)
      DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equations_set_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations set rows = ",EQUATIONS_MAPPING% &
          & TOTAL_NUMBER_OF_ROWS,ERR,ERROR,*999)
        DO row_idx=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set row : ",row_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows mapped to = ",SOLUTION_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS, &
            & ERR,ERROR,*999)
          IF(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)% &
            & NUMBER_OF_SOLVER_ROWS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS,5,5,SOLUTION_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%SOLVER_ROWS, &
              & '("        Solver row numbers    :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS,5,5,SOLUTION_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%COUPLING_COEFFICIENTS, &
              & '("        Coupling coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !row_idx
      ENDDO !equations_set_idx            
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Column mappings:",ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of solver matrices = ",SOLUTION_MAPPING% &
        & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
      DO solver_matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver columns to equations sets columns mappings:",ERR,ERROR,*999)        
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",SOLUTION_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of DOFs = ",SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
          & solver_matrix_idx)%NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of DOFs = ",SOLUTION_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global DOFs = ",SOLUTION_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_GLOBAL_DOFS,ERR,ERROR,*999)
        DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Equations set index : ",equations_set_idx,ERR,ERROR,*999)          
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Column mappings:",ERR,ERROR,*999)        
          DO column_idx=1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS           
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver column : ",column_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of linear equations matrices mapped to = ", &
              & SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS(column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*999)
            IF(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS(column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                & column_idx)%EQUATIONS_MATRIX_NUMBERS,'("      Equation matrices numbers  :",5(X,I13))','(34X,5(X,I13))', &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                & column_idx)%EQUATIONS_COL_NUMBERS,'("      Equation column numbers    :",5(X,I13))','(34X,5(X,I13))', &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
                & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS( &
                & column_idx)%COUPLING_COEFFICIENTS,'("      Coupling coefficients      :",5(X,E13.6))','(34X,5(X,E13.6))', &
                & ERR,ERROR,*999)
            ENDIF
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian column number     = ", &
              & SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS(column_idx)%JACOBIAN_COL_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian coupling coeff    = ", &
              & SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & equations_set_idx)%SOLVER_COL_TO_EQUATIONS_MAPS(column_idx)%JACOBIAN_COUPLING_COEFFICIENT,ERR,ERROR,*999)
          ENDDO !column_idx
        ENDDO !equations_set_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Variable mappings:",ERR,ERROR,*999)        
        DO dof_idx=1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS     
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver dof : ",dof_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations sets mapped to = ",SOLUTION_MAPPING% &
            & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
            & NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
          IF(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
            & NUMBER_OF_EQUATIONS_SETS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLUTION_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%EQUATIONS_SET_INDICES, &
              & '("      Equations set indices  :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
!!TODO: write out the variables somehow.
            !CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
            !  & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(column_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLUTION_MAPPING% &
            !  & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%VARIABLE_TYPE, &
            !  & '("      Variable types         :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLUTION_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%VARIABLE_DOF, &
              & '("      Variable dofs          :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLUTION_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
              & VARIABLE_COEFFICIENT,'("      Variable coefficients  :",5(X,E13.6))','(28X,5(X,E13.6))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS_SETS,5,5,SOLUTION_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
              & ADDITIVE_CONSTANT,'("      Additive constants     :",5(X,E13.6))','(28X,5(X,E13.6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !dof_idx
      ENDDO !solver_matrix_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equations sets columns to solver columns mappings:",ERR,ERROR,*999)
      DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Equations set index : ",equations_set_idx,ERR,ERROR,*999)          
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix indexing:",ERR,ERROR,*999)
        DO solver_matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Linear equations matrix columns to solver matrix columns:", &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of liner equations matrices = ",SOLUTION_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*999)
            DO equations_matrix_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
              EQUATIONS_TO_SOLVER_MAP=>SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations matrix index : ",equations_matrix_idx,ERR,ERROR,*999)
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
            JACOBIAN_TO_SOLVER_MAP=>SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
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
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Variable dofs to solver matrix dofs:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",SOLUTION_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
            & NUMBER_OF_VARIABLES,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_VARIABLES,5,5,SOLUTION_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
            & VARIABLE_TYPES,'("      Variable types :",5(X,I13))','(21X,5(X,I13))',ERR,ERROR,*999)
          DO variable_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_VARIABLES
            DEPENDENT_VARIABLE=>SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLES(variable_idx)%PTR
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable index : ",variable_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of variable dofs = ",DEPENDENT_VARIABLE% &
              & NUMBER_OF_DOFS,ERR,ERROR,*999)
            DO local_dof=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Variable dof : ",local_dof,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Solver column number = ",SOLUTION_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COLUMN_NUMBERS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Coupling coefficient = ",SOLUTION_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COUPLING_COEFFICIENTS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Additive constant    = ",SOLUTION_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%ADDITIVE_CONSTANTS(local_dof),ERR,ERROR,*999)              
            ENDDO !local_dof
          ENDDO !variable_idx
        ENDDO !solver_matrix_idx
        IF(ASSOCIATED(LINEAR_MAPPING)) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Linear equations matrix indexing:",ERR,ERROR,*999)
          DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations matrix : ",equations_matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of solver matrices = ",SOLUTION_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)% &
              & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
            DO solver_matrix_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES
              EQUATIONS_TO_SOLVER_MAP=>SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS(solver_matrix_idx)%PTR
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
          JACOBIAN_TO_SOLVER_MAP=>SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM%JACOBIAN_TO_SOLVER_MATRIX_MAP
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",JACOBIAN_TO_SOLVER_MAP% &
            & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
        ENDIF
      ENDDO ! equations_set_idx       
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_CALCULATE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_CALCULATE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CALCULATE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solution mapping
  SUBROUTINE SOLUTION_MAPPING_CREATE_FINISH(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solution mapping has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
          CALL SOLUTION_MAPPING_CALCULATE(SOLUTION_MAPPING,ERR,ERROR,*999)
          CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
          SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED=.TRUE.            
        ELSE
          CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_FINISH")
    RETURN
999 CALL SOLUTION_MAPPING_FINALISE(SOLUTION_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solution mapping for a problem solution
  SUBROUTINE SOLUTION_MAPPING_CREATE_START(SOLUTION,SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to create the solution mapping from.
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(SOLUTION%SOLUTION_FINISHED) THEN
        IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
          CALL FLAG_ERROR("Solution mapping is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(SOLUTION_MAPPING)
          CALL SOLUTION_MAPPING_INITIALISE(SOLUTION,ERR,ERROR,*999)
          SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_CREATE_START",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_START")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises a solution mapping create values cache and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solution mapping create values cache 
  SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Solution mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping create values cache.",ERR,ERROR,*999)
        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
          & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping create values cache matrix variable types.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping create values cache residual variable type.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping create values cache RHS variable type.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping create values cache source variable type.", &
          & ERR,ERROR,*999)
        SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=0
        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE=0
        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=0
        SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Destroy a solution mapping.
  SUBROUTINE SOLUTION_MAPPING_DESTROY(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer the solution mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      CALL SOLUTION_MAPPING_FINALISE(SOLUTION_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_DESTROY")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_DESTROY",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a equations column to solver columns map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(EQUATIONS_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE) :: EQUATIONS_COL_SOLVER_COLS_MAP !<The equations col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_COL_SOLVER_COLS_MAP%SOLVER_COLS)) &
      & DEALLOCATE(EQUATIONS_COL_SOLVER_COLS_MAP%SOLVER_COLS)
    IF(ALLOCATED(EQUATIONS_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(EQUATIONS_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations column to solver columns map
  SUBROUTINE SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE(EQUATIONS_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE) :: EQUATIONS_COL_SOLVER_COLS_MAP !<The equations column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_COL_SOLVER_COLS_MAP%NUMBER_OF_SOLVER_COLS=0
    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the mapping of global variables to a solver matrix for the solution mapping
  SUBROUTINE SOLUTION_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET(SOLUTION_MAPPING,SOLVER_MATRIX,EQUATIONS_SET_INDEX, &
    & VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(IN) :: SOLVER_MATRIX !<The solver matrix number to set the equations variables for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_INDEX !<The equations set index in the solution mapping to specify the variable types for
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
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solution mappings has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(SOLVER_MATRIX>=1.AND.SOLVER_MATRIX<=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
            IF(EQUATIONS_SET_INDEX>=1.AND.EQUATIONS_SET_INDEX<=SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS) THEN
              EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(EQUATIONS_SET_INDEX)%PTR
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
                            LOCAL_ERROR="The variable type number of "// &
                              & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))// &
                              & " at position "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                              & " in the array is invalid. The number must be >=1 and <= "// &
                              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                          IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(VARIABLE_TYPES(variable_idx))% &
                            & NUMBER_OF_LINEAR_EQUATIONS_MATRICES==0) THEN
                            LOCAL_ERROR="The variable type number of "// &
                              & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))// &
                              & " at position "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                              & " in the array is invalid. That variable type is not mapped to any equations matrices"
                          ENDIF
                        ENDDO !variable_idx
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,EQUATIONS_SET_INDEX,SOLVER_MATRIX)= &
                          & SIZE(VARIABLE_TYPES,1)
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:SIZE(VARIABLE_TYPES,1),EQUATIONS_SET_INDEX, &
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
                & TRIM(NUMBER_TO_VSTRING(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The solver matrix number of "//TRIM(NUMBER_TO_VSTRING(SOLVER_MATRIX,"*",ERR,ERROR))// &
              & " is invalid. The number must be >= 1 and <= "// &
              & TRIM(NUMBER_TO_VSTRING(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises a equations row to solver rows map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE(EQUATIONS_ROW_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_ROW_TO_SOLVER_ROWS_MAP_TYPE) :: EQUATIONS_ROW_SOLVER_ROWS_MAP !<The equations row to solver rows map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_ROW_SOLVER_ROWS_MAP%SOLVER_ROWS)) &
      & DEALLOCATE(EQUATIONS_ROW_SOLVER_ROWS_MAP%SOLVER_ROWS)
    IF(ALLOCATED(EQUATIONS_ROW_SOLVER_ROWS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(EQUATIONS_ROW_SOLVER_ROWS_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations row to solver rows map
  SUBROUTINE SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE(EQUATIONS_ROW_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_ROW_TO_SOLVER_ROWS_MAP_TYPE) :: EQUATIONS_ROW_SOLVER_ROWS_MAP !<The equations row to solver rows map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_ROW_SOLVER_ROWS_MAP%NUMBER_OF_SOLVER_ROWS=0
    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds an equations set to a solution mapping
  SUBROUTINE SOLUTION_MAPPING_EQUATIONS_SET_ADD(SOLUTION_MAPPING,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer the solution mapping to add the equations set to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_INDEX !<On exit, the index of the equations set in the solution mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,matrix_idx,variable_idx,variable_type
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_VARIABLE_TYPES(:,:,:),OLD_RHS_VARIABLE_TYPE(:),OLD_RESIDUAL_VARIABLE_TYPE(:), &
      & OLD_SOURCE_VARIABLE_TYPE(:)
    LOGICAL :: MATRIX_DONE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_SET_PTR_TYPE), ALLOCATABLE :: OLD_EQUATIONS_SETS(:)
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLUTION_MAPPING_EQUATIONS_SET_ADD",ERR,ERROR,*999)

    EQUATIONS_SET_INDEX=0
    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solution mapping has been finished.",ERR,ERROR,*999)
      ELSE
        SOLUTION=>SOLUTION_MAPPING%SOLUTION
        IF(ASSOCIATED(SOLUTION)) THEN
          IF(ASSOCIATED(EQUATIONS_SET)) THEN
            IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN                
                  IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
                    !Check the equations set linearity matches the solution linearity
                    IF(EQUATIONS_SET%LINEARITY==SOLUTION%LINEARITY) THEN
                      IF(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS>0) THEN
                        ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES,SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types.",ERR,ERROR,*999)
                        ALLOCATE(OLD_RESIDUAL_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old residual variable type.",ERR,ERROR,*999)
                        ALLOCATE(OLD_RHS_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old RHS variable type.",ERR,ERROR,*999)
                        ALLOCATE(OLD_SOURCE_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old source variable type.",ERR,ERROR,*999)
                        OLD_MATRIX_VARIABLE_TYPES=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
                        OLD_RESIDUAL_VARIABLE_TYPE=SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE
                        OLD_RHS_VARIABLE_TYPE=SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE
                        OLD_SOURCE_VARIABLE_TYPE=SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE
                        ALLOCATE(OLD_EQUATIONS_SETS(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old equations sets.",ERR,ERROR,*999)
                        DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                          OLD_EQUATIONS_SETS(equations_set_idx)%PTR=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                        ENDDO !equations_set_idx
                        DEALLOCATE(SOLUTION_MAPPING%EQUATIONS_SETS)
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                          & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types.",ERR,ERROR,*999)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate residual variable type.",ERR,ERROR,*999)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate RHS variable type.",ERR,ERROR,*999)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate source variable type.",ERR,ERROR,*999)
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS,:)=OLD_MATRIX_VARIABLE_TYPES
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS)=OLD_RESIDUAL_VARIABLE_TYPE
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS)=OLD_RHS_VARIABLE_TYPE
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS)=OLD_SOURCE_VARIABLE_TYPE
                      ELSE IF(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS==0) THEN
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)
                        IF(ALLOCATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)) &
                          & DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                          & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types.",ERR,ERROR,*999)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate residual variable type.",ERR,ERROR,*999)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate RHS variable type.",ERR,ERROR,*999)
                        ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate source variable type.",ERR,ERROR,*999)
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS,:)=0
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS)=0
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS)=0
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(1:SOLUTION_MAPPING% &
                          & NUMBER_OF_EQUATIONS_SETS)=0
                      ELSE
                        CALL FLAG_ERROR("The number of equations sets is < 0.",ERR,ERROR,*999)
                      ENDIF
                      IF(ASSOCIATED(EQUATIONS_MAPPING%LINEAR_MAPPING)) THEN
                        !Linear matrices to map. 
                        !Map the first matrix variable found in the equations set to the first solver matrix, the second
                        !variable found to the second, etc.
                        variable_type=1
                        DO matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
                          MATRIX_DONE=.FALSE.
                          DO WHILE(variable_type<=FIELD_NUMBER_OF_VARIABLE_TYPES.AND..NOT.MATRIX_DONE)
                            IF(EQUATIONS_MAPPING%LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN                  
                              SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0, &
                                & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,matrix_idx)=1
                              SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1, &
                                & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1,matrix_idx)=variable_type
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
                            & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
                            LOCAL_ERROR="Variable type "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                              & " is mapped to a linear matrix but has not been mapped to any solver matrices."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                        ENDDO !variable_idx
                      ENDIF
                      IF(ASSOCIATED(EQUATIONS_MAPPING%NONLINEAR_MAPPING)) THEN
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                          & EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                      ENDIF
                      IF(ASSOCIATED(EQUATIONS_MAPPING%RHS_MAPPING)) THEN
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                          & EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE
                      ENDIF
                      IF(ASSOCIATED(EQUATIONS_MAPPING%SOURCE_MAPPING)) THEN
                        SOLUTION_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)= &
                          & EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_VARIABLE_TYPE
                      ENDIF
                      ALLOCATE(SOLUTION_MAPPING%EQUATIONS_SETS(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations sets.",ERR,ERROR,*999)
                      DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR=>OLD_EQUATIONS_SETS(equations_set_idx)%PTR
                      ENDDO !equations_set_idx
                      SOLUTION_MAPPING%EQUATIONS_SETS(SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)%PTR=>EQUATIONS_SET
                      SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS=SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS+1
                      EQUATIONS_SET_INDEX=SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
                      IF(ALLOCATED(OLD_RESIDUAL_VARIABLE_TYPE)) DEALLOCATE(OLD_RESIDUAL_VARIABLE_TYPE)
                      IF(ALLOCATED(OLD_RHS_VARIABLE_TYPE)) DEALLOCATE(OLD_RHS_VARIABLE_TYPE)
                      IF(ALLOCATED(OLD_SOURCE_VARIABLE_TYPE)) DEALLOCATE(OLD_SOURCE_VARIABLE_TYPE)
                      IF(ALLOCATED(OLD_EQUATIONS_SETS)) DEALLOCATE(OLD_EQUATIONS_SETS)
                    ELSE
                      LOCAL_ERROR="The specified equations set linearity type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%LINEARITY,"*",ERR,ERROR))// &
                        & " does not match the solution linearity type of "// &
                        & TRIM(NUMBER_TO_VSTRING(SOLUTION%LINEARITY,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Solutions mapping create values cache is not associated.",ERR,ERROR,*999)
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
          CALL FLAG_ERROR("Solution mapping solution is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_SET_ADD")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
    IF(ALLOCATED(OLD_RESIDUAL_VARIABLE_TYPE)) DEALLOCATE(OLD_RESIDUAL_VARIABLE_TYPE)
    IF(ALLOCATED(OLD_RHS_VARIABLE_TYPE)) DEALLOCATE(OLD_RHS_VARIABLE_TYPE)
    IF(ALLOCATED(OLD_SOURCE_VARIABLE_TYPE)) DEALLOCATE(OLD_SOURCE_VARIABLE_TYPE)
    IF(ALLOCATED(OLD_EQUATIONS_SETS)) DEALLOCATE(OLD_EQUATIONS_SETS)
    CALL ERRORS("SOLUTION_MAPPING_EQUATIONS_SET_ADD",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_SET_ADD")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATIONS_SET_ADD

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_matrix_idx,row_idx,solver_matrix_idx
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM)) THEN
      DO solver_matrix_idx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM,1)
        CALL SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx),ERR,ERROR,*999)
      ENDDO !solver_matrix_idx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM)
    ENDIF
    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM)) THEN
      DO equations_matrix_idx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM,1)
        CALL SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx),ERR,ERROR,*999)
      ENDDO !equations_matrix_idx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM)
    ENDIF
    CALL SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,ERR,ERROR,*999)
    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS)) THEN
      DO row_idx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS,1)
        CALL SOLUTION_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx),ERR,ERROR,*999)
      ENDDO !row_idx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS)
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver map.
  SUBROUTINE SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_SET_INDEX=0
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%SOLUTION_MAPPING)
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS)
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)
        
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE(EQUATIONS_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
      IF(ALLOCATED(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP)) THEN
        DO column_idx=1,SIZE(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP,1)
          CALL SOLUTION_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(EQUATIONS_TO_SOLVER_MAP% &
            & EQUATIONS_COL_SOLVER_COLS_MAP(column_idx),ERR,ERROR,*999)
        ENDDO !column_idx
        DEALLOCATE(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP)
      ENDIF
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver maps
  SUBROUTINE SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE(EQUATIONS_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP !<The equations to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
      EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX_NUMBER=0
      EQUATIONS_TO_SOLVER_MAP%SOLVER_MATRIX_NUMBER=0
      NULLIFY(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX)
      NULLIFY(EQUATIONS_TO_SOLVER_MAP%SOLVER_MATRIX)
    ELSE
      CALL FLAG_ERROR("Equations to solver map is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map em and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM !<The equations set to solver matrix maps em to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_TO_SOLVER_MATRIX_MAPS)) THEN
      DO matrix_idx=1,SIZE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_TO_SOLVER_MATRIX_MAPS,1)
        CALL SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS(matrix_idx)%PTR,ERR,ERROR,*999)        
      ENDDO !variable_idx
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_TO_SOLVER_MATRIX_MAPS)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps em.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM !<The equations to solver matrix maps em to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%EQUATIONS_MATRIX_NUMBER=0
    EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM%NUMBER_OF_SOLVER_MATRICES=0
        
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map jm and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM !<The equations set to solver matrix maps jm to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)) THEN
      CALL SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM%JACOBIAN_TO_SOLVER_MATRIX_MAP, &
        & ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps jm.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to initialise the euqations_to_solver_matrix_map_jm for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM)) THEN
      CALL FLAG_ERROR("Equations to solver matrix maps jm is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps jm.",ERR,ERROR,*999)
      NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM%JACOBIAN_TO_SOLVER_MATRIX_MAP)
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map sm and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM !<The equations set to solver matrix maps sm to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,variable_idx
    
    CALL ENTERS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TYPES)) DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TYPES)
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLES)) DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLES)
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TO_SOLVER_COL_MAPS)) THEN
      DO variable_idx=1,SIZE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TO_SOLVER_COL_MAPS,1)
        CALL SOLUTION_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM% &
          & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx),ERR,ERROR,*999)        
      ENDDO !variable_idx
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%VARIABLE_TO_SOLVER_COL_MAPS)
    ENDIF
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%EQUATIONS_TO_SOLVER_MATRIX_MAPS)) THEN
      DO matrix_idx=1,SIZE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%EQUATIONS_TO_SOLVER_MATRIX_MAPS,1)
        CALL SOLUTION_MAPPING_EQUATIONS_TO_SOLVER_MAPS_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS(matrix_idx)%PTR,ERR,ERROR,*999)        
      ENDDO !variable_idx
      DEALLOCATE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%EQUATIONS_TO_SOLVER_MATRIX_MAPS)
    ENDIF
    CALL SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%JACOBIAN_TO_SOLVER_MATRIX_MAP, &
      & ERR,ERROR,*999)
    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps sm.
  SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM_TYPE) :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM !<The equations to solver matrix maps sm to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%SOLVER_MATRIX_NUMBER=0
    EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_VARIABLES=0
    EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
    NULLIFY(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM%JACOBIAN_TO_SOLVER_MATRIX_MAP)
        
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solution mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_FINALISE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,row_idx,solver_matrix_idx

    CALL ENTERS("SOLUTION_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ALLOCATED(SOLUTION_MAPPING%EQUATIONS_SETS)) DEALLOCATE(SOLUTION_MAPPING%EQUATIONS_SETS)        
      IF(ALLOCATED(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)) THEN
        DO equations_set_idx=1,SIZE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP,1)
          CALL SOLUTION_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
            & equations_set_idx),ERR,ERROR,*999)
        ENDDO !equations_set_idx
        DEALLOCATE(SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)
      ENDIF
      IF(ALLOCATED(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP)) THEN
        DO solver_matrix_idx=1,SIZE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP,1)
          CALL SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP( &
            & solver_matrix_idx),ERR,ERROR,*999)
        ENDDO !solver_matrix_idx
        DEALLOCATE(SOLUTION_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP)
      ENDIF
      IF(ALLOCATED(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS)) THEN
        DO row_idx=1,SIZE(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS,1)
          CALL SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS( &
            & row_idx),ERR,ERROR,*999)
        ENDDO !row_idx
        DEALLOCATE(SOLUTION_MAPPING%SOLVER_ROW_TO_EQUATIONS_SET_MAPS)
      ENDIF
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(SOLUTION_MAPPING%ROW_DOFS_MAPPING,ERR,ERROR,*999)
      CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(SOLUTION_MAPPING)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solution mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_INITIALISE(SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(ASSOCIATED(SOLUTION%SOLUTION_MAPPING)) THEN
        CALL FLAG_ERROR("Solution solution mapping is already associated",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLUTION%SOLUTION_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution solution mapping",ERR,ERROR,*999)
        SOLUTION%SOLUTION_MAPPING%SOLUTION=>SOLUTION
        SOLUTION%SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED=.FALSE.
        SOLUTION%SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES=1
        SOLUTION%SOLUTION_MAPPING%NUMBER_OF_ROWS=0
        SOLUTION%SOLUTION_MAPPING%NUMBER_OF_GLOBAL_ROWS=0
        SOLUTION%SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS=0
        NULLIFY(SOLUTION%SOLUTION_MAPPING%ROW_DOFS_MAPPING)
        NULLIFY(SOLUTION%SOLUTION_MAPPING%CREATE_VALUES_CACHE)
        CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE(SOLUTION%SOLUTION_MAPPING,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_FINALISE(SOLUTION%SOLUTION_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a Jacobian column to solver columns map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE(JACOBIAN_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE) :: JACOBIAN_COL_SOLVER_COLS_MAP !<The Jacobian col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(JACOBIAN_COL_SOLVER_COLS_MAP%SOLVER_COLS)) &
      & DEALLOCATE(JACOBIAN_COL_SOLVER_COLS_MAP%SOLVER_COLS)
    IF(ALLOCATED(JACOBIAN_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(JACOBIAN_COL_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an Jacobian column to solver columns map
  SUBROUTINE SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE(JACOBIAN_COL_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE) :: JACOBIAN_COL_SOLVER_COLS_MAP !<The Jacobian column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR,*999)

    JACOBIAN_COL_SOLVER_COLS_MAP%NUMBER_OF_SOLVER_COLS=0
    
    CALL EXITS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(JACOBIAN_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP !<The jacobian to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx
    
    CALL ENTERS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
      IF(ALLOCATED(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP)) THEN
        DO column_idx=1,SIZE(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP,1)
          CALL SOLUTION_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE(JACOBIAN_TO_SOLVER_MAP% &
            & JACOBIAN_COL_SOLVER_COLS_MAP(column_idx),ERR,ERROR,*999)
        ENDDO !column_idx
        DEALLOCATE(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP)
      ENDIF
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a Jacobian to solver maps
  SUBROUTINE SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE(JACOBIAN_TO_SOLVER_MATRIX_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MATRIX_MAP !<The Jacobian to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MATRIX_MAP)) THEN
      JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX_NUMBER=0
      NULLIFY(JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_MATRIX)
      NULLIFY(JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX)
    ELSE
      CALL FLAG_ERROR("Jacobian to solver matrix map is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_FINALISE(JACOBIAN_TO_SOLVER_MATRIX_MAP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_JACOBIAN_TO_SOLVER_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_MAP !<The solver column to equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_MAP%EQUATIONS_MATRIX_NUMBERS)) &
      & DEALLOCATE(SOLVER_COL_TO_EQUATIONS_MAP%EQUATIONS_MATRIX_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_MAP%EQUATIONS_COL_NUMBERS)) &
      & DEALLOCATE(SOLVER_COL_TO_EQUATIONS_MAP%EQUATIONS_COL_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(SOLVER_COL_TO_EQUATIONS_MAP%COUPLING_COEFFICIENTS)
       
    CALL EXITS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_INITIALISE(SOLVER_COL_TO_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_MAP !<The solver column to equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_EQUATIONS_MAP%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
    SOLVER_COL_TO_EQUATIONS_MAP%JACOBIAN_COL_NUMBER=0    
    SOLVER_COL_TO_EQUATIONS_MAP%JACOBIAN_COUPLING_COEFFICIENT=0.0_DP
    
    CALL EXITS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_COL_TO_EQUATIONS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations set map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SET_MAP !<The solver column to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_EQUATIONS_MAPS)) &
      & DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_EQUATIONS_MAPS)
       
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations set mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE(SOLVER_COL_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SET_MAP !<The solver column to equations set map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR,*999)

    NULLIFY(SOLVER_COL_TO_EQUATIONS_SET_MAP%EQUATIONS)
    
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations sets map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SETS_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SETS_MAP !<The solver column to equations sets map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: col_idx,equations_set_idx
    
    CALL ENTERS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_COL_TO_EQUATIONS_SET_MAPS)) THEN
      DO equations_set_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_COL_TO_EQUATIONS_SET_MAPS,1)
        CALL SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP% &
          SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx),ERR,ERROR,*999)
      ENDDO !equations_set_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_COL_TO_EQUATIONS_SET_MAPS)
    ENDIF
    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS)) THEN
      DO col_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS,1)
        CALL SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS( &
          & col_idx),ERR,ERROR,*999)
      ENDDO !col_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_DOF_TO_VARIABLE_MAPS)
    ENDIF
    CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP%COLUMN_DOFS_MAPPING,ERR,ERROR,*999)
    
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations sets mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE(SOLVER_COL_TO_EQUATIONS_SETS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_SETS_MAP_TYPE) :: SOLVER_COL_TO_EQUATIONS_SETS_MAP !<The solver column to equations sets map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_MATRIX_NUMBER=0
    NULLIFY(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLVER_MATRIX)
    NULLIFY(SOLVER_COL_TO_EQUATIONS_SETS_MAP%SOLUTION_MAPPING)
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%NUMBER_OF_COLUMNS=0
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%NUMBER_OF_DOFS=0
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%TOTAL_NUMBER_OF_DOFS=0
    SOLVER_COL_TO_EQUATIONS_SETS_MAP%NUMBER_OF_GLOBAL_DOFS=0
    NULLIFY(SOLVER_COL_TO_EQUATIONS_SETS_MAP%COLUMN_DOFS_MAPPING)
    
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOL_COL_TO_EQUATS_SETS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver dof to variable mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE(SOLVER_DOF_TO_VARIABLE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_DOF_TO_VARIABLE_MAP_TYPE) :: SOLVER_DOF_TO_VARIABLE_MAP !<The solver dof to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_SET_INDICES)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_SET_INDICES)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_DOF)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_DOF)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_COEFFICIENT)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%VARIABLE_COEFFICIENT)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%ADDITIVE_CONSTANT)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%ADDITIVE_CONSTANT)
    
    CALL EXITS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver dof to variable mapping.
  SUBROUTINE SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE(SOLVER_DOF_TO_VARIABLE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_DOF_TO_VARIABLE_MAP_TYPE) :: SOLVER_DOF_TO_VARIABLE_MAP !<The solver dof to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_DOF_TO_VARIABLE_MAP%NUMBER_OF_EQUATIONS_SETS=0
    
    CALL EXITS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solver matrices for the solution mapping
  SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLUTION_MAPPING,NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SOLVER_MATRICES !<The number of solver matrices for the solution.
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

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solution mappings has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
          MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES=1
          DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
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
            IF(NUMBER_OF_SOLVER_MATRICES/=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
              ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types",ERR,ERROR,*999)
              OLD_MATRIX_VARIABLE_TYPES=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
              DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
              ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                & SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS,NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types",ERR,ERROR,*999)
              IF(NUMBER_OF_SOLVER_MATRICES>SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
                SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,:,1:SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES)= &
                  & OLD_MATRIX_VARIABLE_TYPES(:,:,1:SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES)
                DO matrix_idx=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES+1,NUMBER_OF_SOLVER_MATRICES
                  SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,:,matrix_idx)=1
                  SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1,:,matrix_idx)=FIELD_STANDARD_VARIABLE_TYPE
                  SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(2:FIELD_NUMBER_OF_VARIABLE_TYPES,:,matrix_idx)=0
                ENDDO !matrix_idx
              ELSE
                SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,:,1:NUMBER_OF_SOLVER_MATRICES)= &
                  & OLD_MATRIX_VARIABLE_TYPES(:,:,1:NUMBER_OF_SOLVER_MATRICES)
              ENDIF
              SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES=NUMBER_OF_SOLVER_MATRICES
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
          CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
    CALL ERRORS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises a variable to solver column map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE(SOLVER_ROW_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_ROW_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_ROW_TO_EQUATIONS_SET_MAP !<The solver row to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_SET)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_SET)
    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_ROW_NUMBER)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_SET_MAP%EQUATIONS_ROW_NUMBER)
    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_SET_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_SET_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE")
    RETURN 1
    
  END SUBROUTINE SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solver row to equations set map.
  SUBROUTINE SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE(SOLVER_ROW_TO_EQUATIONS_SET_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_ROW_TO_EQUATIONS_SET_MAP_TYPE) :: SOLVER_ROW_TO_EQUATIONS_SET_MAP !<The solver row to equations set map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR,*999)

   SOLVER_ROW_TO_EQUATIONS_SET_MAP%NUMBER_OF_ROWS=0
        
    CALL EXITS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN 1
    
  END SUBROUTINE SOLUTION_MAPPING_SOL_ROW_TO_EQUATS_SET_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a variable to solver column map and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE(VARIABLE_TO_SOLVER_COL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VARIABLE_TO_SOLVER_COL_MAP_TYPE) :: VARIABLE_TO_SOLVER_COL_MAP !<The variable to solver column map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLUTION_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(VARIABLE_TO_SOLVER_COL_MAP%COLUMN_NUMBERS)) DEALLOCATE(VARIABLE_TO_SOLVER_COL_MAP%COLUMN_NUMBERS)
    IF(ALLOCATED(VARIABLE_TO_SOLVER_COL_MAP%COUPLING_COEFFICIENTS)) DEALLOCATE(VARIABLE_TO_SOLVER_COL_MAP%COUPLING_COEFFICIENTS)
    IF(ALLOCATED(VARIABLE_TO_SOLVER_COL_MAP%ADDITIVE_CONSTANTS)) DEALLOCATE(VARIABLE_TO_SOLVER_COL_MAP%ADDITIVE_CONSTANTS)
        
    CALL EXITS("SOLUTION_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE

  !
  !================================================================================================================================
  !
  
END MODULE SOLUTION_MAPPING_ROUTINES

