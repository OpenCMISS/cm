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
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup SOLVER_MAPPING_EquationsMatrixTypes SOLVER_MAPPING::EquationsMatrixTypes
  !> \brief Equations matrix types
  !> \see SOLVER_MAPPING
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX=1 !<The equations matrix in the solver mapping is a dynamic equations matrix \see SOLVER_MAPPING_EquationsMatrixTypes,SOLVER_MAPPING
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX=2 !<The equations matrix in the solver mapping is a linear equations matrix \see SOLVER_MAPPING_EquationsMatrixTypes,SOLVER_MAPPING
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX=3 !<The equations matrix in the solver mapping is a nonlinear equations (Jacobian) matrix \see SOLVER_MAPPING_EquationsMatrixTypes,SOLVER_MAPPING
  !>@}
 
  !> \addtogroup SOLVER_MAPPING_EquationsTypes SOLVER_MAPPING::EquationsTypes
  !> \brief Equations Matrix types
  !> \see SOLVER_MAPPING
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET=1 !<The equations in the solver mapping is from an equations set \see SOLVER_MAPPING_EquationsTypes,SOLVER_MAPPING
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION=2 !<The equations in the solver mapping is from an interface condition \see SOLVER_MAPPING_EquationsTypes,SOLVER_MAPPING
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE=3 !<The equations in the solver mapping is from a transposed interface condition \see SOLVER_MAPPING_EquationsTypes,SOLVER_MAPPING
  !>@}
 
  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX,SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX

  PUBLIC SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION

  PUBLIC SOLVER_MAPPING_CREATE_FINISH,SOLVER_MAPPING_CREATE_START

  PUBLIC SOLVER_MAPPING_DESTROY
  
  PUBLIC SOLVER_MAPPING_EQUATIONS_SET_ADD

  PUBLIC SOLVER_MAPPING_INTERFACE_CONDITION_ADD

  PUBLIC SOLVER_MAPPING_EQUATS_VARS_TO_SOLVER_MATRIX_SET

  PUBLIC SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
CONTAINS

  !
  !=================================================================================================================================
  !
  
  !>Calculates the solver mappings
  SUBROUTINE SOLVER_MAPPING_CALCULATE(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,COLUMN_LIST_ITEM(4),COLUMN_RANK,DEPENDENT_VARIABLE_TYPE,dof_idx,dof_type,equation_type, &
      & equations_column,equations_idx,equations_idx2,equations_matrix,equations_matrix_idx,equations_row_number, &
      & equations_set_idx,EQUATIONS_VARIABLE_LIST_ITEM(3),global_column,global_dof,global_dof_idx,GLOBAL_DOFS_OFFSET, &
      & global_row,global_row_idx,interface_column,interface_col_number,interface_condition_idx,interface_condition_idx2, &
      & INTERFACE_EQUATIONS_LIST_ITEM(2),interface_idx,interface_matrix_idx,interface_row,interface_row_number,jacobian_column, &
      & local_column,local_dof,LOCAL_DOFS_OFFSET,local_row,matrices_type,matrix_number,matrix_type,matrix_type_idx, &
      & matrix_variable_idx,myrank,NUMBER_OF_COLUMNS,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,NUMBER_OF_EQUATIONS_COLUMNS, &
      & NUMBER_OF_EQUATIONS_SETS,NUMBER_OF_EQUATIONS_VARIABLES,NUMBER_OF_INTERFACES,NUMBER_OF_INTERFACE_COLUMNS, &
      & NUMBER_OF_INTERFACE_ROWS,NUMBER_OF_INTERFACE_VARIABLES,NUMBER_OF_GLOBAL_SOLVER_DOFS,NUMBER_OF_GLOBAL_SOLVER_ROWS, &
      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,NUMBER_OF_LOCAL_SOLVER_DOFS,NUMBER_OF_LOCAL_SOLVER_ROWS,NUMBER_OF_RANK_COLS, &
      & NUMBER_OF_RANK_ROWS,NUMBER_OF_VARIABLES,rank,rank_idx,row_idx,ROW_LIST_ITEM(3),ROW_RANK,solver_global_dof, &
      & solver_local_dof,solver_matrix_idx,solver_variable_idx,TOTAL_NUMBER_OF_LOCAL_SOLVER_DOFS,variable_idx, &
      & VARIABLE_LIST_ITEM(3),variable_position_idx,variable_type
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_SET_VARIABLES(:,:),EQUATIONS_VARIABLES(:,:),INTERFACE_EQUATIONS_LIST(:,:), &
      & INTERFACE_VARIABLES(:,:),RANK_GLOBAL_ROWS_LIST(:,:),RANK_GLOBAL_COLS_LIST(:,:)
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(:),NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(:), &
      & TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(:),SUB_MATRIX_INFORMATION(:,:,:),SUB_MATRIX_LIST(:,:,:),VARIABLE_TYPES(:)
    LOGICAL :: FOUND,INCLUDE_COLUMN,INCLUDE_ROW
    LOGICAL, ALLOCATABLE :: VARIABLE_PROCESSED(:),VARIABLE_RANK_PROCESSED(:,:)
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
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,LAGRANGE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,LAGRANGE_VARIABLE,VARIABLE
    TYPE(INTEGER_INTG_PTR_TYPE), POINTER :: DOF_MAP(:)
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_TO_SOLVER_MAPS_TYPE), POINTER :: INTERFACE_TO_SOLVER_MAP
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(LIST_TYPE), POINTER :: EQUATIONS_SET_VARIABLE_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: INTERFACE_EQUATIONS_LISTS(:),RANK_GLOBAL_ROWS_LISTS(:,:), &
      & RANK_GLOBAL_COLS_LISTS(:,:,:,:),VARIABLES_LIST(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MAPPING_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
        SOLVER_EQUATIONS=>SOLVER_MAPPING%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          !
          !--- Equations set <-> interface conditions  ---
          !
          !Allocate equations set to solver map
          ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping equations set to solver map.",ERR,ERROR,*999)      
          !Allocate interface condition to solver map
          ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping interface condition to solver map.",ERR,ERROR,*999)

          ALLOCATE(INTERFACE_EQUATIONS_LISTS(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set list.",ERR,ERROR,*999)
          DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
            INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                CALL SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_INITIALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                  & interface_condition_idx),ERR,ERROR,*999)
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_CONDITION_INDEX= &
                  & interface_condition_idx
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%SOLVER_MAPPING=>SOLVER_MAPPING
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
                NULLIFY(INTERFACE_EQUATIONS_LISTS(interface_condition_idx)%PTR)
                CALL LIST_CREATE_START(INTERFACE_EQUATIONS_LISTS(interface_condition_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(INTERFACE_EQUATIONS_LISTS(interface_condition_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_DATA_DIMENSION_SET(INTERFACE_EQUATIONS_LISTS(interface_condition_idx)%PTR,2,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(INTERFACE_EQUATIONS_LISTS(interface_condition_idx)%PTR,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !interface_condition_idx
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                CALL SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                  & equations_set_idx),ERR,ERROR,*999)
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_SET_INDEX=equations_set_idx
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%SOLVER_MAPPING=>SOLVER_MAPPING
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS=>EQUATIONS
                !Set up list of interface conditions affecting this equations set
                CALL LIST_DETACH_AND_DESTROY(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR, &
                  & NUMBER_OF_INTERFACES,INTERFACE_EQUATIONS_LIST,ERR,ERROR,*999)
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(NUMBER_OF_INTERFACES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver maps interface.",ERR,ERROR,*999)
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%NUMBER_OF_INTERFACE_CONDITIONS=NUMBER_OF_INTERFACES
                DO interface_idx=1,NUMBER_OF_INTERFACES
                  CALL SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                    & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(interface_idx),ERR,ERROR,*999)
                  interface_condition_idx=INTERFACE_EQUATIONS_LIST(1,interface_idx)
                  interface_matrix_idx=INTERFACE_EQUATIONS_LIST(2,interface_idx)
                  INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
                  IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE( &
                      & interface_idx)%INTERFACE_CONDITION_INDEX=interface_condition_idx
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE( &
                      & interface_idx)%INTERFACE_CONDITION=>INTERFACE_CONDITION
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE( &
                      & interface_idx)%INTERFACE_MATRIX_NUMBER=interface_matrix_idx
                    INTERFACE_EQUATIONS_LIST_ITEM(1)=equations_set_idx
                    INTERFACE_EQUATIONS_LIST_ITEM(2)=interface_matrix_idx
                    CALL LIST_ITEM_ADD(INTERFACE_EQUATIONS_LISTS(interface_condition_idx)%PTR,INTERFACE_EQUATIONS_LIST_ITEM, &
                      & ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !interface_condition_idx
                IF(ALLOCATED(INTERFACE_EQUATIONS_LIST)) DEALLOCATE(INTERFACE_EQUATIONS_LIST)
              ELSE
                CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !equations_set_idx
          DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
            CALL LIST_DETACH_AND_DESTROY(INTERFACE_EQUATIONS_LISTS(interface_condition_idx)%PTR,NUMBER_OF_EQUATIONS_SETS, &
              INTERFACE_EQUATIONS_LIST,ERR,ERROR,*999)
            ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface to solver maps equations.",ERR,ERROR,*999)
            SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%NUMBER_OF_EQUATIONS_SETS= &
              & NUMBER_OF_EQUATIONS_SETS
            DO equations_idx=1,NUMBER_OF_EQUATIONS_SETS
              CALL SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_INITIALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_idx),ERR,ERROR,*999)
              equations_set_idx=INTERFACE_EQUATIONS_LIST(1,equations_idx)
              interface_matrix_idx=INTERFACE_EQUATIONS_LIST(2,equations_idx)
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_idx)%EQUATIONS_SET_INDEX=equations_set_idx
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_idx)%EQUATIONS_SET=>EQUATIONS_SET
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_idx)%INTERFACE_MATRIX_INDEX=interface_matrix_idx
              ELSE
                CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !equations_idx
            IF(ALLOCATED(INTERFACE_EQUATIONS_LIST)) DEALLOCATE(INTERFACE_EQUATIONS_LIST)
          ENDDO !interface_condition_idx
          !
          !--- Row mappings ---
          !
          !Calculate the row mappings.
          !We do not have any couplings defined at the moment there is only a 1-1 mapping.
          myrank=COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER
          NUMBER_OF_GLOBAL_SOLVER_ROWS=0
          NUMBER_OF_LOCAL_SOLVER_ROWS=0
          !Add in the rows from any equations sets that have been added to the solver equations
          !Presort the row numbers by rank.
          ALLOCATE(RANK_GLOBAL_ROWS_LISTS(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+SOLVER_MAPPING% &
            & NUMBER_OF_INTERFACE_CONDITIONS,0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocat rank global rows lists.",ERR,ERROR,*999)
          DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            equations_idx=0
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              equations_idx=equations_idx+1
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                EQUATIONS=>EQUATIONS_SET%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING                
                  IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                    NULLIFY(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR)
                    CALL LIST_CREATE_START(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,INT(EQUATIONS_MAPPING% &
                      & NUMBER_OF_GLOBAL_ROWS/COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,INTG), &
                      & ERR,ERROR,*999)
                    CALL LIST_DATA_DIMENSION_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,3,ERR,ERROR,*999)
                    CALL LIST_KEY_DIMENSION_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,1,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,ERR,ERROR,*999)                    
                  ELSE
                    CALL FLAG_ERROR("Equations equations mapping is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !equations_set_idx
            DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
              equations_idx=equations_idx+1
              INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
              IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                SELECT CASE(INTERFACE_CONDITION%METHOD)
                CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                 
                  INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
                  IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                    INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                      NULLIFY(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR)
                      CALL LIST_CREATE_START(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR, &
                        & INT(INTERFACE_MAPPING%NUMBER_OF_GLOBAL_COLUMNS/COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES, &
                        & INTG),ERR,ERROR,*999)
                      CALL LIST_DATA_DIMENSION_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,3,ERR,ERROR,*999)
                      CALL LIST_KEY_DIMENSION_SET(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,1,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,ERR,ERROR,*999)
                      
                     ELSE
                      CALL FLAG_ERROR("Interface equations interface mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
                  ENDIF
                CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The interface condition method of "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT                
              ELSE
                CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !interface_condition_idx
          ENDDO !rank
          !Calculate the number of local and global rows.
          equations_idx=0
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            equations_idx=equations_idx+1
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
                        !Loop over the global rows for this equations set
                        DO global_row=1,EQUATIONS_MAPPING%NUMBER_OF_GLOBAL_ROWS
                          !Find the rank that owns this global row
                          ROW_RANK=-1
                          DO rank_idx=1,ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%NUMBER_OF_DOMAINS
                            IF(ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%LOCAL_TYPE(rank_idx)/=DOMAIN_LOCAL_GHOST) THEN
                              ROW_RANK=ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%DOMAIN_NUMBER(rank_idx)
                              local_row=ROW_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_row)%LOCAL_NUMBER(rank_idx)
                              EXIT
                            ENDIF
                          ENDDO !rank_idx
                          IF(ROW_RANK>=0) THEN
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
                                  & BOUNDARY_CONDITION_FREE.OR.BOUNDARY_CONDITIONS_VARIABLE% & 
                                  & GLOBAL_BOUNDARY_CONDITIONS(global_dof)==BOUNDARY_CONDITION_FREE_WALL
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
                                    & BOUNDARY_CONDITION_FREE.OR.BOUNDARY_CONDITIONS_VARIABLE% & 
                                    & GLOBAL_BOUNDARY_CONDITIONS(global_dof)==BOUNDARY_CONDITION_FREE_WALL 
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
                                      & global_dof)==BOUNDARY_CONDITION_FREE.OR.BOUNDARY_CONDITIONS_VARIABLE% & 
                                      & GLOBAL_BOUNDARY_CONDITIONS(global_dof)==BOUNDARY_CONDITION_FREE_WALL 
                                  ELSE
                                    CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !matrix_idx
                              ENDIF
                            ENDIF
                            ROW_LIST_ITEM(1)=global_row
                            ROW_LIST_ITEM(2)=local_row
                            IF(INCLUDE_ROW) THEN
                              ROW_LIST_ITEM(3)=1
                              NUMBER_OF_GLOBAL_SOLVER_ROWS=NUMBER_OF_GLOBAL_SOLVER_ROWS+1
                              !Don't need to worry about ghosted rows.
                              IF(ROW_RANK==myrank) NUMBER_OF_LOCAL_SOLVER_ROWS=NUMBER_OF_LOCAL_SOLVER_ROWS+1 !1-1 mapping
                            ELSE
                              ROW_LIST_ITEM(3)=0
                            ENDIF !include row
                            CALL LIST_ITEM_ADD(RANK_GLOBAL_ROWS_LISTS(equations_idx,ROW_RANK)%PTR,ROW_LIST_ITEM,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Global row is not owned by a domain.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !global_row
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
          !Now add in rows from any interface matrices
          DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
            equations_idx=equations_idx+1
            INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
                IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                  INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                  IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                    COL_DOFS_MAPPING=>INTERFACE_MAPPING%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(COL_DOFS_MAPPING)) THEN                    
                      DO global_column=1,INTERFACE_MAPPING%NUMBER_OF_GLOBAL_COLUMNS
                        !Find the rank that owns this global column
                        COLUMN_RANK=-1
                        DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_column)%NUMBER_OF_DOMAINS
                          IF(COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_column)%LOCAL_TYPE(rank_idx)/=DOMAIN_LOCAL_GHOST) THEN
                            COLUMN_RANK=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_column)%DOMAIN_NUMBER(rank_idx)
                            local_column=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_column)%LOCAL_NUMBER(rank_idx)
                            EXIT
                          ENDIF
                        ENDDO !rank_idx
                        IF(COLUMN_RANK>=0) THEN
                          INCLUDE_COLUMN=.TRUE.
                          ROW_LIST_ITEM(1)=global_column
                          ROW_LIST_ITEM(2)=local_column
                          IF(INCLUDE_COLUMN) THEN
                            ROW_LIST_ITEM(3)=1
                            NUMBER_OF_GLOBAL_SOLVER_ROWS=NUMBER_OF_GLOBAL_SOLVER_ROWS+1
                            !Don't need to worry about ghosted rows.
                            IF(COLUMN_RANK==myrank) NUMBER_OF_LOCAL_SOLVER_ROWS=NUMBER_OF_LOCAL_SOLVER_ROWS+1 !1-1 mapping
                          ELSE
                            ROW_LIST_ITEM(3)=0                          
                          ENDIF !include column
                          CALL LIST_ITEM_ADD(RANK_GLOBAL_ROWS_LISTS(equations_idx,COLUMN_RANK)%PTR,ROW_LIST_ITEM,ERR,ERROR,*999)
                        ELSE
                          CALL FLAG_ERROR("Global row is not owned by a domain.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO !global_column
                    ELSE
                      CALL FLAG_ERROR("Interface condition column degree of freedom mappings is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface equations interface mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
                ENDIF
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !interface_condition_idx

          IF(NUMBER_OF_LOCAL_SOLVER_ROWS==0) &
            & CALL FLAG_ERROR("Invalid problem setup. The number of local solver rows is zero.",ERR,ERROR,*999)
          IF(NUMBER_OF_GLOBAL_SOLVER_ROWS==0) &
            & CALL FLAG_ERROR("Invalid problem setup. The number of global solver rows is zero.",ERR,ERROR,*999)

          !Allocate memory for the rows mapping
          !Allocate the solver rows to equations set maps
          ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping solver row to equation rows map.",ERR,ERROR,*999)
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

          !Initialise the equations sets to solver maps
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS

            !Note that pointers have been checked for association above
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            EQUATIONS=>EQUATIONS_SET%EQUATIONS
            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
            
            !Allocate the equations set to solver maps for solver matrix (sm) indexing
            ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
              & SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set to solver map equations to solver matrix maps sm.", &
              & ERR,ERROR,*999)
            DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
              CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_SM_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx),ERR,ERROR,*999)
            ENDDO !solver_matrix_idx
            
            !Allocate the equations row to solver rows maps
            ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
              & EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set to solver map equations row to solver rows maps.", &
              & ERR,ERROR,*999)
            DO equations_row_number=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
              !Initialise
              CALL SOLVER_MAPPING_EQUATS_ROW_TO_SOL_ROWS_MAP_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number),ERR,ERROR,*999)
            ENDDO
            
          ENDDO !equations_set_idx
          
          !Initialise the interface condition to solver maps
          DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
            
            !Note that pointers have been checked for association above
            INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
            INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
            
            !Allocate the interface to solver maps for solver matrix (sm) indexing
            ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
              & SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface to solver map interface to solver matrix maps sm.", &
              & ERR,ERROR,*999)
            DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
              CALL SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx),ERR,ERROR,*999)
            ENDDO !solver_matrix_idx
            
            !Allocate the interface to solver maps for interface matrix (im) indexing
            ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
            IF(ERR/=0) &
              & CALL FLAG_ERROR("Could not allocate interface to solver map equations to solver matrix maps im.", &
              & ERR,ERROR,*999)
            DO interface_matrix_idx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
              CALL SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_INITIALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx),ERR,ERROR,*999)
                      
              !Allocate the interfafce row to solver row maps
              ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP( &
                & INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%TOTAL_NUMBER_OF_ROWS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface condition to solver map interface row to solver row map.", &
                & ERR,ERROR,*999)
              DO interface_row_number=1,INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx) &
                & %TOTAL_NUMBER_OF_ROWS
                CALL SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_INITIALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                  & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                  & INTERFACE_ROW_TO_SOLVER_ROWS_MAP(interface_row_number),ERR,ERROR,*999)                
              ENDDO !interface_row_number
              
            ENDDO !interface_matrix_idx

            !Allocate the interface column to solver row maps
            ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS),STAT=ERR)
            IF(ERR/=0)  &
              & CALL FLAG_ERROR("Could not allocate interface condition to solver map interface column to solver row map.", &
              & ERR,ERROR,*999)
            DO interface_col_number=1,INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS
              !Initialise
              CALL SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_INITIALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                & interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(interface_col_number),ERR,ERROR,*999)
            ENDDO !interface_col_number
            
          ENDDO !interface condition_idx
          
          !Calculate the row mappings
          NUMBER_OF_GLOBAL_SOLVER_ROWS=0
          !Loop over the ranks to  ensure that the lowest ranks have the lowest numbered solver variables
          DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            NUMBER_OF_LOCAL_SOLVER_ROWS=0
            equations_idx=0
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              equations_idx=equations_idx+1

              !Get rows list
              CALL LIST_SORT(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,NUMBER_OF_RANK_ROWS, &
                & RANK_GLOBAL_ROWS_LIST,ERR,ERROR,*999)

              !Note that pointers have been checked for association above
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
              DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING

              !Loop over the global rows for this rank.
              DO global_row_idx=1,NUMBER_OF_RANK_ROWS
                global_row=RANK_GLOBAL_ROWS_LIST(1,global_row_idx)
                local_row=RANK_GLOBAL_ROWS_LIST(2,global_row_idx)
                INCLUDE_ROW=RANK_GLOBAL_ROWS_LIST(3,global_row_idx)==1
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
                  IF(rank==myrank) THEN
                    !If this is my rank then set up the solver->equations and equations->solver row mappings
                    
                    !Set up the solver row -> equations row mappings. 1-1 mapping as no coupling at the moment. Will need to look
                    !At the interface conditions for this equations set later.
                    !Initialise
                    CALL SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP( &
                      & NUMBER_OF_LOCAL_SOLVER_ROWS),ERR,ERROR,*999)

                    ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_INDEX(1), &
                      & STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver row to equations rows equations index.",ERR,ERROR,*999)
                    ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%ROWCOL_NUMBER(1),STAT=ERR)
                    IF(ERR/=0) &
                      & CALL FLAG_ERROR("Could not allocate solver row to equations rows row/col number.",ERR,ERROR,*999)
                    ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)% &
                      & COUPLING_COEFFICIENTS(1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver row to equations rows coupling coefficients.", &
                      & ERR,ERROR,*999)
                    !Set the mappings
                    SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%NUMBER_OF_EQUATIONS_SETS=1
                    SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%EQUATIONS_INDEX(1)= &
                      & equations_set_idx
                    SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%ROWCOL_NUMBER(1)= &
                      & local_row
                    SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%COUPLING_COEFFICIENTS(1)= &
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
                    !Now set up any interface condition rows to solver rows that affect this equations set.
                    DO interface_condition_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                      & NUMBER_OF_INTERFACE_CONDITIONS
                      interface_condition_idx2=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(interface_condition_idx)%INTERFACE_CONDITION_INDEX
                      interface_matrix_idx=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(interface_condition_idx)%INTERFACE_MATRIX_NUMBER
                      !Set the mappings
                      SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx2)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP( &
                        & local_row)%NUMBER_OF_SOLVER_ROWS=1
                      SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx2)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP( &
                        & local_row)%SOLVER_ROW=NUMBER_OF_LOCAL_SOLVER_ROWS
                      SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx2)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP( &
                        & local_row)%COUPLING_COEFFICIENT=1.0_DP                      
                    ENDDO !interface_condition_idx
                    
                  ENDIF !rank==my rank
                ELSE
                  IF(rank==myrank) THEN
                    !Set up the equations row -> solver row mappings
                    !Set the mappings
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                      & local_row)%NUMBER_OF_SOLVER_ROWS=0
                    !Now set up any interface condition rows to solver rows that affect this equations set.
                    DO interface_condition_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                      & NUMBER_OF_INTERFACE_CONDITIONS
                      interface_condition_idx2=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(interface_condition_idx)%INTERFACE_CONDITION_INDEX
                      interface_matrix_idx=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(interface_condition_idx)%INTERFACE_MATRIX_NUMBER
                      !Set the mappings
                      SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx2)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP( &
                        & local_row)%NUMBER_OF_SOLVER_ROWS=0
                    ENDDO !interface_condition_idx
                  ENDIF !rank==my rank
                ENDIF !include row
              ENDDO !global_row_idx
              IF(ALLOCATED(RANK_GLOBAL_ROWS_LIST)) DEALLOCATE(RANK_GLOBAL_ROWS_LIST)
            ENDDO !equations_set_idx
            DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
              equations_idx=equations_idx+1

              !Get rows list
              CALL LIST_SORT(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(RANK_GLOBAL_ROWS_LISTS(equations_idx,rank)%PTR,NUMBER_OF_RANK_ROWS, &
                & RANK_GLOBAL_ROWS_LIST,ERR,ERROR,*999)

              !Note that pointers have been checked for association above
              INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING

              !Loop over the global rows for this rank.
              DO global_row_idx=1,NUMBER_OF_RANK_ROWS
                global_column=RANK_GLOBAL_ROWS_LIST(1,global_row_idx)
                local_column=RANK_GLOBAL_ROWS_LIST(2,global_row_idx)
                INCLUDE_COLUMN=RANK_GLOBAL_ROWS_LIST(3,global_row_idx)==1
                IF(INCLUDE_COLUMN) THEN
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
                    !Set up the solver row -> interface column mappings.
                    !Initialise
                    CALL SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP( &
                      & NUMBER_OF_LOCAL_SOLVER_ROWS),ERR,ERROR,*999)

                    ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%ROWCOL_NUMBER(1),STAT=ERR)
                    IF(ERR/=0) &
                      & CALL FLAG_ERROR("Could not allocate solver row to equations rows row/col number.",ERR,ERROR,*999)
                    ALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)% &
                      & COUPLING_COEFFICIENTS(1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver row to equations rows coupling coefficients.", &
                      & ERR,ERROR,*999)
                    !Set the mappings
                    SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%INTERFACE_CONDITION_INDEX= &
                      & interface_condition_idx
                    SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%ROWCOL_NUMBER(1)= &
                      & local_column
                    SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(NUMBER_OF_LOCAL_SOLVER_ROWS)%COUPLING_COEFFICIENTS(1)= &
                      & 1.0_DP
                    !Set up the interface col -> solver row mappings
                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                      & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(local_column)%NUMBER_OF_SOLVER_ROWS=1
                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                      & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(local_column)%SOLVER_ROW=NUMBER_OF_LOCAL_SOLVER_ROWS
                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                      & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(local_column)%COUPLING_COEFFICIENT=1.0_DP
                  ENDIF !rank==my rank
                ELSE
                  IF(rank==myrank) THEN
                    !Set the interface column -> solver row mappings
                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                      & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(local_column)%NUMBER_OF_SOLVER_ROWS=1
                  ENDIF
                ENDIF
              ENDDO !global_row_idx
              IF(ALLOCATED(RANK_GLOBAL_ROWS_LIST)) DEALLOCATE(RANK_GLOBAL_ROWS_LIST)              
            ENDDO !interface_condition_idx            
          ENDDO !rank
          CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(ROW_DOMAIN_MAPPING,ERR,ERROR,*999)
          !
          !--- Column mappings ---
          !
          !Allocate solver column to equations sets mapping array
          ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping solver column to equations column maps.",ERR,ERROR,*999)
          ALLOCATE(SOLVER_MAPPING%VARIABLES_LIST(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping variables list.",ERR,ERROR,*999)
          
          !Calculate the column mappings for each solver matrix
          DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES

            !Initialise the variables list
            CALL SOLVER_MAPPING_VARIABLES_INITIALISE(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx),ERR,ERROR,*999)
            !Compute the order of variables for the solver matrices          
            CALL LIST_DETACH_AND_DESTROY(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR, &
              & NUMBER_OF_EQUATIONS_VARIABLES,EQUATIONS_VARIABLES,ERR,ERROR,*999)
            CALL LIST_DETACH_AND_DESTROY(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR, &
              & NUMBER_OF_INTERFACE_VARIABLES,INTERFACE_VARIABLES,ERR,ERROR,*999)
            ALLOCATE(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(NUMBER_OF_EQUATIONS_VARIABLES+ &
              & NUMBER_OF_INTERFACE_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables list variables.",ERR,ERROR,*999)
            SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES=NUMBER_OF_EQUATIONS_VARIABLES+ &
              & NUMBER_OF_INTERFACE_VARIABLES
            ALLOCATE(VARIABLES_LIST(NUMBER_OF_EQUATIONS_VARIABLES+NUMBER_OF_INTERFACE_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables list.",ERR,ERROR,*999)
            ALLOCATE(VARIABLE_PROCESSED(NUMBER_OF_EQUATIONS_VARIABLES+NUMBER_OF_INTERFACE_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable processed.",ERR,ERROR,*999)
            VARIABLE_PROCESSED=.FALSE.
            solver_variable_idx=0
            DO variable_idx=1,NUMBER_OF_EQUATIONS_VARIABLES
              solver_variable_idx=solver_variable_idx+1
              CALL SOLVER_MAPPING_VARIABLE_INITIALISE(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)% &
                & VARIABLES(solver_variable_idx),ERR,ERROR,*999)
              equations_set_idx=EQUATIONS_VARIABLES(1,variable_idx)
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  variable_type=EQUATIONS_VARIABLES(2,variable_idx)
                  VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  IF(ASSOCIATED(VARIABLE)) THEN
                    SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)%VARIABLE=> &
                      & VARIABLE
                    SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)%VARIABLE_TYPE= &
                      & VARIABLE_TYPE
                    NULLIFY(VARIABLES_LIST(solver_variable_idx)%PTR)
                    CALL LIST_CREATE_START(VARIABLES_LIST(solver_variable_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(VARIABLES_LIST(solver_variable_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_DATA_DIMENSION_SET(VARIABLES_LIST(solver_variable_idx)%PTR,3,ERR,ERROR,*999)
                    CALL LIST_KEY_DIMENSION_SET(VARIABLES_LIST(solver_variable_idx)%PTR,1,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(VARIABLES_LIST(solver_variable_idx)%PTR,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !variable_idx
            IF(ALLOCATED(EQUATIONS_VARIABLES)) DEALLOCATE(EQUATIONS_VARIABLES)
            DO variable_idx=1,NUMBER_OF_INTERFACE_VARIABLES
              solver_variable_idx=solver_variable_idx+1
              CALL SOLVER_MAPPING_VARIABLE_INITIALISE(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)% &
                & VARIABLES(solver_variable_idx),ERR,ERROR,*999)
              interface_condition_idx=INTERFACE_VARIABLES(1,variable_idx)
              INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
              IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
                  LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
                  IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                    variable_type=INTERFACE_VARIABLES(2,variable_idx)
                    VARIABLE=>LAGRANGE_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                    IF(ASSOCIATED(VARIABLE)) THEN
                      SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)%VARIABLE=> &
                        & VARIABLE
                      SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)%VARIABLE_TYPE= &
                        & VARIABLE_TYPE
                      NULLIFY(VARIABLES_LIST(solver_variable_idx)%PTR)
                      CALL LIST_CREATE_START(VARIABLES_LIST(solver_variable_idx)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(VARIABLES_LIST(solver_variable_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_DATA_DIMENSION_SET(VARIABLES_LIST(solver_variable_idx)%PTR,3,ERR,ERROR,*999)
                      CALL LIST_KEY_DIMENSION_SET(VARIABLES_LIST(solver_variable_idx)%PTR,1,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(VARIABLES_LIST(solver_variable_idx)%PTR,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Lagrange field variable is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface condition Lagrange field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !variable_idx
            IF(ALLOCATED(INTERFACE_VARIABLES)) DEALLOCATE(INTERFACE_VARIABLES)
            
            !Initialise solver column to equations sets mapping array
            CALL SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_INITIALISE(SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx),ERR,ERROR,*999)
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_MATRIX_NUMBER=solver_matrix_idx
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_MAPPING=>SOLVER_MAPPING
            !Allocate the solver col to equations set maps array
            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
              & SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver col to equations map solver col to equation set maps.", &
              & ERR,ERROR,*999)
            !Allocate the solver col to interface maps array
            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_INTERFACE_MAPS( &
              & SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver col to equations map solver col to interface maps.", &
              & ERR,ERROR,*999)
            !Presort the column numbers by rank.
            ALLOCATE(RANK_GLOBAL_COLS_LISTS(2,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+SOLVER_MAPPING% &
              & NUMBER_OF_INTERFACE_CONDITIONS,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES, &
              & 0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate rank global columns lists.",ERR,ERROR,*999)
            DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
              DO solver_variable_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
                DO equations_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
                  DO dof_type=1,2
                    NULLIFY(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR)
                    CALL LIST_CREATE_START(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR, &
                      & ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR, &
                      & LIST_INTG_TYPE,ERR,ERROR,*999)
                    !Set approximate size for the number of columns per variable.
                    CALL LIST_INITIAL_SIZE_SET(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR, &
                      & SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)%VARIABLE% &
                      & TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
                    CALL LIST_DATA_DIMENSION_SET(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR,4, &
                      & ERR,ERROR,*999)
                    CALL LIST_KEY_DIMENSION_SET(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR,1, &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR, &
                      & ERR,ERROR,*999)
                  ENDDO !dof_type
                ENDDO !equations_idx
              ENDDO !solver_variable_set_idx
            ENDDO !rank

            !Allocate sub-matrix information
            ALLOCATE(SUB_MATRIX_INFORMATION(3,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+ &
              & SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)% &
              & NUMBER_OF_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub matrix information.",ERR,ERROR,*999)
            SUB_MATRIX_INFORMATION=0
            !Allocate sub-matrix list information
            ALLOCATE(SUB_MATRIX_LIST(0:3,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+ &
              & SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)% &
              & NUMBER_OF_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub matrix list.",ERR,ERROR,*999)
            SUB_MATRIX_LIST=0
            
            !Calculate the number of solver dofs
            ALLOCATE(NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES), &
              & STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of global solver dofs.",ERR,ERROR,*999)
            ALLOCATE(NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES), &
              & STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of local solver dofs.",ERR,ERROR,*999)
            ALLOCATE(TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)% &
              & NUMBER_OF_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate total number of local solver dofs.",ERR,ERROR,*999)

            NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS=0
            NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS=0
            TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS=0
            
            equations_idx=0
            !Loop over the equations sets
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              equations_idx=equations_idx+1
              !The pointers below have been checked for association above.
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
              DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
               LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
              NULLIFY(EQUATIONS_SET_VARIABLE_LIST)
              CALL LIST_CREATE_START(EQUATIONS_SET_VARIABLE_LIST,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(EQUATIONS_SET_VARIABLE_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_DATA_DIMENSION_SET(EQUATIONS_SET_VARIABLE_LIST,3,ERR,ERROR,*999)
              CALL LIST_KEY_DIMENSION_SET(EQUATIONS_SET_VARIABLE_LIST,1,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(EQUATIONS_SET_VARIABLE_LIST,ERR,ERROR,*999)
              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                EQUATIONS_VARIABLE_LIST_ITEM(1)=SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(equations_set_idx)
                EQUATIONS_VARIABLE_LIST_ITEM(2)=SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
                EQUATIONS_VARIABLE_LIST_ITEM(3)=0
                CALL LIST_ITEM_ADD(EQUATIONS_SET_VARIABLE_LIST,EQUATIONS_VARIABLE_LIST_ITEM,ERR,ERROR,*999)
              ENDIF
              IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                DO variable_idx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,equations_set_idx,solver_matrix_idx)
                  EQUATIONS_VARIABLE_LIST_ITEM(1)=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES( &
                    & variable_idx,equations_set_idx,solver_matrix_idx)
                  EQUATIONS_VARIABLE_LIST_ITEM(2)=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
                  EQUATIONS_VARIABLE_LIST_ITEM(3)=variable_idx                
                  CALL LIST_ITEM_ADD(EQUATIONS_SET_VARIABLE_LIST,EQUATIONS_VARIABLE_LIST_ITEM,ERR,ERROR,*999)
                ENDDO
              ENDIF
              IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                EQUATIONS_VARIABLE_LIST_ITEM(1)=SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(equations_set_idx)
                EQUATIONS_VARIABLE_LIST_ITEM(2)=SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
                EQUATIONS_VARIABLE_LIST_ITEM(3)=0                
                CALL LIST_ITEM_ADD(EQUATIONS_SET_VARIABLE_LIST,EQUATIONS_VARIABLE_LIST_ITEM,ERR,ERROR,*999)
              ENDIF
              CALL LIST_REMOVE_DUPLICATES(EQUATIONS_SET_VARIABLE_LIST,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(EQUATIONS_SET_VARIABLE_LIST,NUMBER_OF_VARIABLES,EQUATIONS_SET_VARIABLES,ERR,ERROR,*999)
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
              !Loop over the variables in this equations set.
              DO variable_idx=1,NUMBER_OF_VARIABLES
                variable_type=EQUATIONS_SET_VARIABLES(1,variable_idx)
                matrices_type=EQUATIONS_SET_VARIABLES(2,variable_idx)
                matrix_variable_idx=EQUATIONS_SET_VARIABLES(3,variable_idx)
                DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                  !Find the variable in the list of solver variables
                  FOUND=.FALSE.
                  DO variable_position_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
                    IF(ASSOCIATED(DEPENDENT_VARIABLE,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES( &
                      & variable_position_idx)%VARIABLE)) THEN
                      FOUND=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO !variable_position_idx
                  IF(FOUND) THEN
                    !Add the equations set variable to the list of equations involving the solver variable
                    VARIABLE_LIST_ITEM(1)=equations_set_idx
                    VARIABLE_LIST_ITEM(2)=variable_type
                    VARIABLE_LIST_ITEM(3)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                    CALL LIST_ITEM_ADD(VARIABLES_LIST(variable_position_idx)%PTR,VARIABLE_LIST_ITEM,ERR,ERROR,*999)
                    COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                    IF(ASSOCIATED(COL_DOFS_MAPPING)) THEN
                      BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
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
                        !Set the sub-matrix information
                        SUB_MATRIX_INFORMATION(1,equations_idx,variable_position_idx)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                        SUB_MATRIX_INFORMATION(2,equations_idx,variable_position_idx)=equations_set_idx
                        !Set the sub-matrix lists
                        IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                          NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+ &
                            & DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES
                          IF(DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                            SUB_MATRIX_LIST(0,equations_idx,variable_position_idx)= &
                              SUB_MATRIX_LIST(0,equations_idx,variable_position_idx)+1
                            SUB_MATRIX_LIST(SUB_MATRIX_LIST(0,equations_idx,variable_position_idx), &
                              & equations_idx,variable_position_idx)=SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
                          ENDIF
                        ENDIF
                        IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                          NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES+ &
                            & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES
                          IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                            SUB_MATRIX_LIST(0,equations_idx,variable_position_idx)= &
                              SUB_MATRIX_LIST(0,equations_idx,variable_position_idx)+1
                            SUB_MATRIX_LIST(SUB_MATRIX_LIST(0,equations_idx,variable_position_idx), &
                              & equations_idx,variable_position_idx)=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
                          ENDIF
                        ENDIF
                        IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                          IF(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%VARIABLE_TYPE==variable_type) THEN
                            SUB_MATRIX_LIST(0,equations_idx,variable_position_idx)= &
                              SUB_MATRIX_LIST(0,equations_idx,variable_position_idx)+1
                            SUB_MATRIX_LIST(SUB_MATRIX_LIST(0,equations_idx,variable_position_idx), &
                              & equations_idx,variable_position_idx)=SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
                          ENDIF
                        ENDIF
                        !Loop over the global dofs for this variable.
                        DO global_dof=1,DEPENDENT_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%NUMBER_OF_DOMAINS
                            local_dof=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_NUMBER(rank_idx)
                            dof_type=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_TYPE(rank_idx)
                            COLUMN_RANK=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%DOMAIN_NUMBER(rank_idx)
                            INCLUDE_COLUMN=BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)== &
                              & BOUNDARY_CONDITION_FREE.OR.BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS( &
                              & global_dof)==BOUNDARY_CONDITION_FREE_WALL                       
                            COLUMN_LIST_ITEM(1)=global_dof
                            COLUMN_LIST_ITEM(2)=local_dof
                            IF(dof_type/=DOMAIN_LOCAL_GHOST) THEN
                              !DOF is not a ghost dof
                              IF(INCLUDE_COLUMN) THEN
                                COLUMN_LIST_ITEM(3)=1
                                IF(.NOT.VARIABLE_PROCESSED(variable_position_idx)) THEN
                                  NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(variable_position_idx)= &
                                    & NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(variable_position_idx)+1
                                  IF(COLUMN_RANK==myrank) THEN
                                    NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                      & NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                    TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                      & TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                  ENDIF
                                ENDIF
                              ELSE
                                COLUMN_LIST_ITEM(3)=0
                              ENDIF
                              COLUMN_LIST_ITEM(4)=variable_idx
                              CALL LIST_ITEM_ADD(RANK_GLOBAL_COLS_LISTS(1,equations_idx,variable_position_idx,COLUMN_RANK)%PTR, &
                                & COLUMN_LIST_ITEM,ERR,ERROR,*999)
                            ELSE
                              !DOF is a ghost dof
                              IF(INCLUDE_COLUMN) THEN
                                COLUMN_LIST_ITEM(3)=1
                                IF(.NOT.VARIABLE_PROCESSED(variable_position_idx)) THEN
                                  IF(COLUMN_RANK==myrank) TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                    & TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                ENDIF
                              ELSE
                                COLUMN_LIST_ITEM(3)=0
                              ENDIF
                              COLUMN_LIST_ITEM(4)=variable_idx
                              CALL LIST_ITEM_ADD(RANK_GLOBAL_COLS_LISTS(2,equations_idx,variable_position_idx,COLUMN_RANK)%PTR, &
                                & COLUMN_LIST_ITEM,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !rank_idx
                        ENDDO !global_dof
                      ELSE
                        CALL FLAG_ERROR("Boundary condition variable not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations matrix columns degree of freedom mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    VARIABLE_PROCESSED(variable_position_idx)=.TRUE.
                  ELSE
                    CALL FLAG_ERROR("Dependent variable does not exist in the list of solver variables.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Dependent variable is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !variable_idx
              IF(ALLOCATED(EQUATIONS_SET_VARIABLES)) DEALLOCATE(EQUATIONS_SET_VARIABLES)
            ENDDO !equations_set_idx
            DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
              equations_idx=equations_idx+1
              !The pointers below have been checked for association above.
              INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
                INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
                !Initialise interface condition to solver map (sm)
                CALL SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                  & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx), &
                  & ERR,ERROR,*999)
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                !Allocate the interface to solver map variables arrays
                ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TYPES( &
                  & INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                IF(ERR/=0)  &
                  & CALL FLAG_ERROR("Could not allocate interface to solver matrix maps sm dependent variable types.", &
                  & ERR,ERROR,*999)
                ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLES(INTERFACE_MAPPING% &
                  & NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to solver matrix maps sm variables.",ERR,ERROR,*999)
                ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS( &
                  & INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface to solver matrix maps sm dependent variables "// &
                  & "to solver col maps.",ERR,ERROR,*999)                  
                !First add in the Lagrange to solver variables
                LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
!!TODO: sort out Lagrange variable type
                variable_type=1
                LAGRANGE_VARIABLE=>LAGRANGE_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                IF(ASSOCIATED(LAGRANGE_VARIABLE)) THEN
                  !Setup
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%LAGRANGE_VARIABLE_TYPE=variable_type
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%LAGRANGE_VARIABLE=>LAGRANGE_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  !Find the variable in the list of solver variables
                  FOUND=.FALSE.
                  DO variable_position_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
                    IF(ASSOCIATED(LAGRANGE_VARIABLE,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES( &
                      & variable_position_idx)%VARIABLE)) THEN
                      FOUND=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO !variable_position_idx
                  IF(FOUND) THEN
                    !Add the interface condition variable to the list of equations involving the solver variable
                    VARIABLE_LIST_ITEM(1)=interface_condition_idx
                    VARIABLE_LIST_ITEM(2)=variable_type
                    VARIABLE_LIST_ITEM(3)=SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                    CALL LIST_ITEM_ADD(VARIABLES_LIST(variable_position_idx)%PTR,VARIABLE_LIST_ITEM,ERR,ERROR,*999)
                    COL_DOFS_MAPPING=>LAGRANGE_VARIABLE%DOMAIN_MAPPING
                    IF(ASSOCIATED(COL_DOFS_MAPPING)) THEN
                      SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE_TYPE=variable_type
                      SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE=>LAGRANGE_VARIABLE
                      !Allocate the variable to solver col maps arrays
                      ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP% &
                        & COLUMN_NUMBERS(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps column numbers.", &
                        & ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP% &
                        & COUPLING_COEFFICIENTS(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps coupling coefficients.", &
                        & ERR,ERROR,*999)
                      ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP% &
                        & ADDITIVE_CONSTANTS(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps additive constants.", &
                        & ERR,ERROR,*999)
                      DO equations_idx2=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                        & NUMBER_OF_EQUATIONS_SETS
                        equations_set_idx=SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                          & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_idx2)%EQUATIONS_SET_INDEX
                        interface_matrix_idx=SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                          & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_idx2)%INTERFACE_MATRIX_INDEX
                        !Set the sub-matrix information
                        SUB_MATRIX_INFORMATION(1,equations_set_idx,variable_position_idx)= &
                          & SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                        SUB_MATRIX_INFORMATION(2,equations_set_idx,variable_position_idx)=interface_condition_idx
                        SUB_MATRIX_INFORMATION(3,equations_set_idx,variable_position_idx)=interface_matrix_idx
                        !Loop over the global dofs for this variable.
                        DO global_dof=1,LAGRANGE_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%NUMBER_OF_DOMAINS
                            local_dof=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_NUMBER(rank_idx)
                            dof_type=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_TYPE(rank_idx)
                            COLUMN_RANK=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%DOMAIN_NUMBER(rank_idx)
                            !For now include Lagrange column
                            INCLUDE_COLUMN=.TRUE.
                            COLUMN_LIST_ITEM(1)=global_dof
                            COLUMN_LIST_ITEM(2)=local_dof
                            IF(dof_type/=DOMAIN_LOCAL_GHOST) THEN
                              !DOF is not a ghost dof
                              IF(INCLUDE_COLUMN) THEN
                                COLUMN_LIST_ITEM(3)=1
                                IF(.NOT.VARIABLE_PROCESSED(variable_position_idx)) THEN
                                  NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(variable_position_idx)= &
                                    & NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(variable_position_idx)+1
                                  IF(COLUMN_RANK==myrank) THEN
                                    NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                      & NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                    TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                      & TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                  ENDIF
                                ENDIF
                              ELSE
                                COLUMN_LIST_ITEM(3)=0
                              ENDIF
                              COLUMN_LIST_ITEM(4)=variable_idx
                              CALL LIST_ITEM_ADD(RANK_GLOBAL_COLS_LISTS(1,equations_set_idx,variable_position_idx, &
                                & COLUMN_RANK)%PTR,COLUMN_LIST_ITEM,ERR,ERROR,*999)
                            ELSE
                              !DOF is a ghost dof
                              IF(INCLUDE_COLUMN) THEN
                                COLUMN_LIST_ITEM(3)=1
                                IF(.NOT.VARIABLE_PROCESSED(variable_position_idx)) THEN
                                  IF(COLUMN_RANK==myrank) TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                    & TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                ENDIF
                              ELSE
                                COLUMN_LIST_ITEM(3)=0
                              ENDIF
                              COLUMN_LIST_ITEM(4)=variable_idx
                              CALL LIST_ITEM_ADD(RANK_GLOBAL_COLS_LISTS(2,equations_set_idx,variable_position_idx, &
                                & COLUMN_RANK)%PTR,COLUMN_LIST_ITEM,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !rank_idx
                        ENDDO !global_dof
                        VARIABLE_PROCESSED(variable_position_idx)=.TRUE.
                      ENDDO !equations_idx2
                    ELSE
                      CALL FLAG_ERROR("Columns degree of freedom mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Lagrange variable does not exist in the list of solver variables.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Lagrange variable is not associated.",ERR,ERROR,*999)
                ENDIF
                !Now add in the Dependent variables
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%NUMBER_OF_DEPENDENT_VARIABLES=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                DO interface_matrix_idx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                  DEPENDENT_VARIABLE=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE
                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                    variable_type=DEPENDENT_VARIABLE%VARIABLE_TYPE
                    !Find the variable in the list of solver variables
                    FOUND=.FALSE.
                    DO variable_position_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
                      IF(ASSOCIATED(DEPENDENT_VARIABLE,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES( &
                        & variable_position_idx)%VARIABLE)) THEN
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !variable_position_idx
                    IF(FOUND) THEN
                      EQUATIONS_SET=>INTERFACE_DEPENDENT%EQUATIONS_SETS(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS( &
                        & interface_matrix_idx)%MESH_INDEX)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                          COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                          IF(ASSOCIATED(COL_DOFS_MAPPING)) THEN
                            BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                              & variable_type)%PTR
                            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                              !Setup
                              SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TYPES( &
                                & interface_matrix_idx)=variable_type
                              SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLES(interface_matrix_idx)% &
                                & PTR=>DEPENDENT_VARIABLE
                              ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS( &
                                & interface_matrix_idx)%COLUMN_NUMBERS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps column numbers.", &
                                & ERR,ERROR,*999)
                              ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS( &
                                & interface_matrix_idx)%COUPLING_COEFFICIENTS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                              IF(ERR/=0)  &
                                & CALL FLAG_ERROR("Could not allocate variables to solver column maps coupling coefficients.", &
                                & ERR,ERROR,*999)
                              ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS( &
                                & interface_matrix_idx)%ADDITIVE_CONSTANTS(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variables to solver column maps additive constants.", &
                                & ERR,ERROR,*999)
                              !Set the sub-matrix information
                              SUB_MATRIX_INFORMATION(1,equations_idx,variable_position_idx)= &
                                & SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE
                              SUB_MATRIX_INFORMATION(2,equations_idx,variable_position_idx)=interface_condition_idx
                              SUB_MATRIX_INFORMATION(3,equations_idx,variable_position_idx)=interface_matrix_idx
                              !Loop over the global dofs for this variable.
                              DO global_dof=1,DEPENDENT_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                                DO rank_idx=1,COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%NUMBER_OF_DOMAINS
                                  local_dof=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_NUMBER(rank_idx)
                                  dof_type=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%LOCAL_TYPE(rank_idx)
                                  COLUMN_RANK=COL_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_dof)%DOMAIN_NUMBER(rank_idx)
                                  INCLUDE_COLUMN=BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_dof)== &
                                    & BOUNDARY_CONDITION_FREE.OR.BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS( &
                                    & global_dof)==BOUNDARY_CONDITION_FREE_WALL                    
                                  COLUMN_LIST_ITEM(1)=global_dof
                                  COLUMN_LIST_ITEM(2)=local_dof
                                  IF(dof_type/=DOMAIN_LOCAL_GHOST) THEN
                                    !DOF is not a ghost dof
                                    IF(INCLUDE_COLUMN) THEN
                                      COLUMN_LIST_ITEM(3)=1
                                      IF(.NOT.VARIABLE_PROCESSED(variable_position_idx)) THEN
                                        NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(variable_position_idx)= &
                                          & NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(variable_position_idx)+1
                                        IF(COLUMN_RANK==myrank) THEN
                                          NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                            & NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                          TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                            & TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                        ENDIF
                                      ENDIF
                                    ELSE
                                      COLUMN_LIST_ITEM(3)=0
                                    ENDIF
                                    COLUMN_LIST_ITEM(4)=variable_idx
                                    CALL LIST_ITEM_ADD(RANK_GLOBAL_COLS_LISTS(1,equations_idx,variable_position_idx, &
                                      & COLUMN_RANK)%PTR,COLUMN_LIST_ITEM,ERR,ERROR,*999)
                                  ELSE
                                    !DOF is a ghost dof
                                    IF(INCLUDE_COLUMN) THEN
                                      COLUMN_LIST_ITEM(3)=1
                                      IF(.NOT.VARIABLE_PROCESSED(variable_position_idx)) THEN
                                        IF(COLUMN_RANK==myrank) TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)= &
                                          & TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(variable_position_idx)+1
                                      ENDIF
                                    ELSE
                                      COLUMN_LIST_ITEM(3)=0
                                    ENDIF
                                    COLUMN_LIST_ITEM(4)=variable_idx
                                    CALL LIST_ITEM_ADD(RANK_GLOBAL_COLS_LISTS(2,equations_idx,variable_position_idx, &
                                      & COLUMN_RANK)%PTR,COLUMN_LIST_ITEM,ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !rank_idx
                              ENDDO !global_dof
                            ELSE
                              CALL FLAG_ERROR("Boundary condition variable not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Interface matrix columns degree of freedom mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations set boundary condititions is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Interface depdendent equations set is not associated.",ERR,ERROR,*999)
                      ENDIF
                      VARIABLE_PROCESSED(variable_position_idx)=.TRUE.
                    ELSE
                      CALL FLAG_ERROR("Dependent variable does not exist in the list of solver variables.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Dependent variable is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDDO !interface_idx

            IF(ALLOCATED(VARIABLE_PROCESSED)) DEALLOCATE(VARIABLE_PROCESSED)

            NUMBER_OF_LOCAL_SOLVER_DOFS=0
            TOTAL_NUMBER_OF_LOCAL_SOLVER_DOFS=0
            NUMBER_OF_GLOBAL_SOLVER_DOFS=0
            DO solver_variable_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
              NUMBER_OF_LOCAL_SOLVER_DOFS=NUMBER_OF_LOCAL_SOLVER_DOFS+NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS(solver_variable_idx)
              TOTAL_NUMBER_OF_LOCAL_SOLVER_DOFS=TOTAL_NUMBER_OF_LOCAL_SOLVER_DOFS+TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS( &
                & solver_variable_idx)
              NUMBER_OF_GLOBAL_SOLVER_DOFS=NUMBER_OF_GLOBAL_SOLVER_DOFS+NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS(solver_variable_idx)
            ENDDO !solver_variable_idx

            IF(NUMBER_OF_LOCAL_SOLVER_DOFS==0) THEN
              LOCAL_ERROR="Invalid problem setup. The number of local solver DOFs for solver matrix "// &
                & TRIM(NUMBER_TO_VSTRING(solver_matrix_idx,"*",ERR,ERROR))//" is zero."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
            IF(NUMBER_OF_GLOBAL_SOLVER_DOFS==0) THEN
              LOCAL_ERROR="Invalid problem setup. The number of global solver DOFs for solver matrix "// &
                & TRIM(NUMBER_TO_VSTRING(solver_matrix_idx,"*",ERR,ERROR))//" is zero."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF

            !Allocate memory for this solver matrix
            !Allocate solver columns to equations sets maps
            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
              & TOTAL_NUMBER_OF_LOCAL_SOLVER_DOFS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps.",ERR,ERROR,*999)
            !Set the number of columns
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS=NUMBER_OF_GLOBAL_SOLVER_DOFS
            !Set the number of variables
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%NUMBER_OF_DOFS=NUMBER_OF_LOCAL_SOLVER_DOFS
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS= &
              & TOTAL_NUMBER_OF_LOCAL_SOLVER_DOFS
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%NUMBER_OF_GLOBAL_DOFS=NUMBER_OF_GLOBAL_SOLVER_DOFS
            !Allocate the columns domain mapping
            ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%COLUMN_DOFS_MAPPING,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver col to equations sets map column dofs mapping.",ERR,ERROR,*999)
!!TODO: what is the real number of domains for a solver???
            CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
              & COLUMN_DOFS_MAPPING,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)            
            COL_DOMAIN_MAPPING=>SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%COLUMN_DOFS_MAPPING
            ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column dofs mapping global to local.",ERR,ERROR,*999)
            COL_DOMAIN_MAPPING%NUMBER_OF_GLOBAL=NUMBER_OF_GLOBAL_SOLVER_DOFS
            ALLOCATE(VARIABLE_RANK_PROCESSED(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES, &
              & 0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable rank processed.",ERR,ERROR,*999)
            VARIABLE_RANK_PROCESSED=.FALSE.
            !Calculate the column mappings
            NUMBER_OF_COLUMNS=NUMBER_OF_GLOBAL_SOLVER_DOFS

            !Initialise            
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS

              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
              DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
              
              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                !Allocate the equations set to solver maps for equations matrix (em) indexing
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM( &
                  & DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) &
                  & CALL FLAG_ERROR("Could not allocate equations set to solver map equations to solver matrix maps em.", &
                  & ERR,ERROR,*999)
                DO equations_matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_EM_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                    & equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx),ERR,ERROR,*999)
                ENDDO !equations_matrix_idx
                IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                  !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
                  CALL SOLVER_MAPPING_EQUATS_TO_SOL_MAT_MAPS_JM_INITIALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                    & equations_set_idx),ERR,ERROR,*999)
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
              
              !Initialise solver columns to equations set map
              CALL SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE(SOLVER_MAPPING% &
                & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                & equations_set_idx),ERR,ERROR,*999)
              
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
              
              !Allocate the solver columns to equations set map arrays
              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%HAVE_DYNAMIC=.TRUE.
                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(NUMBER_OF_COLUMNS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver columns to dynamic equations map.",ERR,ERROR,*999)
              ELSE
                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%HAVE_STATIC=.TRUE.
                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(NUMBER_OF_COLUMNS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver columns to static equations map.",ERR,ERROR,*999)
              ENDIF
              !Set the solver column to equations set map
              SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
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
              DO variable_idx=1,NUMBER_OF_VARIABLES
                IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                  variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(equations_set_idx)
                  NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                    & NUMBER_OF_EQUATIONS_MATRICES
                ELSE
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(equations_set_idx)
                  ELSE
                    variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,equations_set_idx, &
                      & solver_matrix_idx)
                    NUMBER_OF_LINEAR_EQUATIONS_MATRICES=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                      & NUMBER_OF_EQUATIONS_MATRICES
                  ENDIF
                ENDIF

                DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%VARIABLES(variable_idx)%PTR=>DEPENDENT_VARIABLE
                IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                  !Allocate dynamic equations to solver matrix maps equations column to solver columns maps
                  DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                    MATRIX_NUMBER=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                      & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                      & SOLVER_MATRIX_NUMBER=solver_matrix_idx
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                      & EQUATIONS_MATRIX_NUMBER=MATRIX_NUMBER
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                      & EQUATIONS_MATRIX=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%EQUATIONS_MATRIX
                    NUMBER_OF_EQUATIONS_COLUMNS=DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                      & equations_matrix_idx)%PTR%EQUATIONS_COL_TO_SOLVER_COLS_MAP(NUMBER_OF_EQUATIONS_COLUMNS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic equations column to solver columns map.", &
                      & ERR,ERROR,*999)
                  ENDDO !equations_matrix_idx
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX_NUMBER=solver_matrix_idx
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%JACOBIAN_MATRIX=>NONLINEAR_MAPPING% &
                      & JACOBIAN_TO_VAR_MAP%JACOBIAN
                    NUMBER_OF_EQUATIONS_COLUMNS=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                      & JACOBIAN_COL_TO_SOLVER_COLS_MAP(NUMBER_OF_EQUATIONS_COLUMNS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Jacobian column to solver columns map.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  !Allocate linear equations to solver matrix maps equations column to solver columns maps
                  IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                    DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                      MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                        & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                        & SOLVER_MATRIX_NUMBER=solver_matrix_idx
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                        & EQUATIONS_MATRIX_NUMBER=MATRIX_NUMBER
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                        & solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                        & EQUATIONS_MATRIX=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%EQUATIONS_MATRIX
                      NUMBER_OF_EQUATIONS_COLUMNS=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                      ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                        & equations_matrix_idx)%PTR%EQUATIONS_COL_TO_SOLVER_COLS_MAP(NUMBER_OF_EQUATIONS_COLUMNS),STAT=ERR)
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
                    NUMBER_OF_EQUATIONS_COLUMNS=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                      & JACOBIAN_COL_TO_SOLVER_COLS_MAP(NUMBER_OF_EQUATIONS_COLUMNS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Jacobian column to solver columns map.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !variable_idx
            ENDDO !equations_set_idx
            
            DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
              
              INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING

              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                                
                !Initialise solver columns to interface condition map
                CALL SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_INITIALISE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                  & solver_matrix_idx)%SOLVER_COL_TO_INTERFACE_MAPS(interface_condition_idx),ERR,ERROR,*999)
                
                !Allocate the solver columns to equations set map arrays
                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_INTERFACE_MAPS( &
                  & interface_condition_idx)%SOLVER_COL_TO_INTERFACE_EQUATIONS_MAPS(NUMBER_OF_COLUMNS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver columns to interface equations map.",ERR,ERROR,*999)
                                
                !Allocate the interface to solver matrix maps sm interface to solver maps
                ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  &INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface to solver matrix maps sm interface equations "// &
                  & "to solver matrix maps.",ERR,ERROR,*999)
                
                !Set up interface arrays
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%NUMBER_OF_INTERFACE_MATRICES=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES                
                DO interface_matrix_idx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                  NULLIFY(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                    & interface_matrix_idx)%PTR)
                  ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                    & interface_matrix_idx)%PTR,STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface to solver matrix maps.",ERR,ERROR,*999)
                  CALL SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_INITIALISE(SOLVER_MAPPING% &
                    & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR,ERR,ERROR,*999)
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                    & interface_matrix_idx)%INTERFACE_MATRIX_NUMBER=interface_matrix_idx
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                    & interface_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=1
                  
                  DEPENDENT_VARIABLE=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%DEPENDENT_VARIABLES(interface_matrix_idx)%PTR=>DEPENDENT_VARIABLE
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                    & SOLVER_MATRIX_NUMBER=solver_matrix_idx
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                    & INTERFACE_MATRIX_NUMBER=interface_matrix_idx
                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                    & INTERFACE_MATRIX=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%INTERFACE_MATRIX
                  NUMBER_OF_INTERFACE_ROWS=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%NUMBER_OF_ROWS
                  ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                    & interface_matrix_idx)%PTR%INTERFACE_ROW_TO_SOLVER_COLS_MAP(NUMBER_OF_INTERFACE_ROWS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface column to solver columns map.",ERR,ERROR,*999)
                ENDDO !interface_matrix_idx
                NUMBER_OF_INTERFACE_COLUMNS=INTERFACE_MAPPING%NUMBER_OF_COLUMNS
                ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_COL_TO_SOLVER_COLS_MAP( &
                  & NUMBER_OF_INTERFACE_COLUMNS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface column to solver columns map.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDDO !interface_condition_idx
            
            !Loop over the ranks to ensure that the lowest ranks have the lowest numbered solver variables

            !Allocate dof map to record column reordering
            ALLOCATE(DOF_MAP(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof map.",ERR,ERROR,*999)
            DO solver_variable_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
              ALLOCATE(DOF_MAP(solver_variable_idx)%PTR(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)% &
                & VARIABLES(solver_variable_idx)%VARIABLE%NUMBER_OF_GLOBAL_DOFS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof map global dof map.",ERR,ERROR,*999)
              DOF_MAP(solver_variable_idx)%PTR=0
            ENDDO !solver_variable_idx
            
            NUMBER_OF_GLOBAL_SOLVER_DOFS=0
            solver_global_dof=0
            solver_local_dof=0
            DO dof_type=1,2
              DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
                
                DO solver_variable_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
                  
                  GLOBAL_DOFS_OFFSET=solver_global_dof
                  LOCAL_DOFS_OFFSET=solver_local_dof
                  
                  variable_type=SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)%VARIABLE_TYPE
                  
                  DO equations_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
                    
                    solver_global_dof=GLOBAL_DOFS_OFFSET
                    solver_local_dof=LOCAL_DOFS_OFFSET
                    
                    !Get columns list
                    CALL LIST_SORT(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(RANK_GLOBAL_COLS_LISTS(dof_type,equations_idx,solver_variable_idx,rank)%PTR, &
                      & NUMBER_OF_RANK_COLS,RANK_GLOBAL_COLS_LIST,ERR,ERROR,*999)
                    
                    IF(NUMBER_OF_RANK_COLS>0) THEN

                      equation_type=SUB_MATRIX_INFORMATION(1,equations_idx,solver_variable_idx)
                      SELECT CASE(equation_type)
                      CASE(SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET)
                        
                        equations_set_idx=SUB_MATRIX_INFORMATION(2,equations_idx,solver_variable_idx)
                       
                        !The pointers below have been checked for association above.
                        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                        EQUATIONS=>EQUATIONS_SET%EQUATIONS
                        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                        DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                        LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS

                        NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
                        NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                        IF(ASSOCIATED(DYNAMIC_MAPPING)) NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=DYNAMIC_MAPPING% &
                          & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES
                        IF(ASSOCIATED(LINEAR_MAPPING)) NUMBER_OF_LINEAR_EQUATIONS_MATRICES=LINEAR_MAPPING% &
                          & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES
                        
                        !Loop over the variables
                        
                        DEPENDENT_VARIABLE=>SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)% &
                          & VARIABLE
                        COL_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                        BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR
                        
                        DO global_dof_idx=1,NUMBER_OF_RANK_COLS
                          global_dof=RANK_GLOBAL_COLS_LIST(1,global_dof_idx)
                          local_dof=RANK_GLOBAL_COLS_LIST(2,global_dof_idx)
                          !dof_type=RANK_GLOBAL_COLS_LIST(3,global_dof_idx)
                          INCLUDE_COLUMN=RANK_GLOBAL_COLS_LIST(3,global_dof_idx)==1                      
                          variable_idx=RANK_GLOBAL_COLS_LIST(4,global_dof_idx)
                          
                          IF(INCLUDE_COLUMN) THEN
                            !DOF is not fixed so map the variable/equation dof to a new solver dof
                            
                            IF(dof_type==2) THEN
                              solver_global_dof=DOF_MAP(solver_variable_idx)%PTR(global_dof)
                            ELSE
                              solver_global_dof=solver_global_dof+1
                              DOF_MAP(solver_variable_idx)%PTR(global_dof)=solver_global_dof
                            ENDIF
                            
                            IF(rank==myrank) THEN
                              
                              solver_local_dof=solver_local_dof+1
                             
                              IF(.NOT.VARIABLE_RANK_PROCESSED(solver_variable_idx,rank)) THEN
                                
                                !Set up the column domain mappings.
                                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                                  & solver_global_dof),ERR,ERROR,*999)
                                !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                                !local map.
                                !Allocate the global to local map arrays
                                ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_NUMBER(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column domain global to local map local number.", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%DOMAIN_NUMBER(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column domain global to local map domain number.", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_TYPE(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column domain global to local map domain number.", &
                                  & ERR,ERROR,*999)
                                !Set up the global to local mappings
                                COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%NUMBER_OF_DOMAINS=1
                                COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_NUMBER(1)= &
                                  & solver_local_dof
                                COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%DOMAIN_NUMBER(1)=rank
                                COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                                
                                !Set up the solver column -> equations column mappings. 1-1 as no coupling yet
!!TODO
                                !Set up the solver dofs -> variable dofs map
                                !Initialise
                                CALL SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE(SOLVER_MAPPING% &
                                  & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof),ERR,ERROR,*999)
                                !Allocate the solver dofs to variable dofs arrays
!!TODO: allow for multiple equations sets in the column
                                !No coupling so there is only one equations set at the moment
                                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                  & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%EQUATIONS_TYPES(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps equations types.", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                  & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%EQUATIONS_INDICES(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps equations indices.", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                  & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%VARIABLE(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable type.", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                  & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%VARIABLE_DOF(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable dof.", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                  & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%VARIABLE_COEFFICIENT(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable coefficient.", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                  & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%ADDITIVE_CONSTANT(1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps additive constant.", &
                                  & ERR,ERROR,*999)
                                !Setup
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof)%NUMBER_OF_EQUATIONS=1
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof)%EQUATIONS_TYPES(1)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof)%EQUATIONS_INDICES(1)=equations_set_idx
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof)%VARIABLE(1)%PTR=>DEPENDENT_VARIABLE
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof)%VARIABLE_DOF(1)=local_dof
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof)%VARIABLE_COEFFICIENT(1)=1.0_DP
                                SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                  & solver_local_dof)%ADDITIVE_CONSTANT(1)=0.0_DP
                              ENDIF
                              !Set up the equations variables -> solver columns mapping
                              !No coupling yet so the mapping is 1-1
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COLUMN_NUMBERS(local_dof)= &
                                & solver_global_dof
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COUPLING_COEFFICIENTS( &
                                & local_dof)=1.0_DP
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%ADDITIVE_CONSTANTS( &
                                & local_dof)=0.0_DP
                              !Set up the equations columns -> solver columns mapping
                              DO matrix_type_idx=1,SUB_MATRIX_LIST(0,equations_idx,solver_variable_idx)
                                matrix_type=SUB_MATRIX_LIST(matrix_type_idx,equations_idx,solver_variable_idx)
                                SELECT CASE(matrix_type)
                                CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                                  !Dynamic matrix
                                  DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                                    MATRIX_NUMBER=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                    equations_column=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(local_dof)
                                    !Allocate the equation to solver map column items.
                                    !No coupling yet so the mapping is 1-1
                                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%SOLVER_COLS(1), &
                                      & STAT=ERR)
                                    IF(ERR/=0) &
                                      & CALL  FLAG_ERROR("Could not allocate dynamic equations column to solver columns map "// &
                                      & "solver colums.",ERR,ERROR,*999)
                                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                                    IF(ERR/=0) &
                                      & CALL FLAG_ERROR("Could not allocate dynamic equations column to solver columns map "// &
                                      & "coupling coefficients.",ERR,ERROR,*999)
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=1
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%SOLVER_COLS(1)=solver_global_dof
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%COUPLING_COEFFICIENTS(1)= &
                                      & DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%MATRIX_COEFFICIENT
                                  ENDDO !equations_matrix_idx
                                CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                                  DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                    MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                    equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(local_dof)
                                    !Allocate the equation to solver map column items.
                                    !No coupling yet so the mapping is 1-1
                                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%SOLVER_COLS(1), &
                                      & STAT=ERR)
                                    IF(ERR/=0) &
                                      & CALL FLAG_ERROR("Could not allocate linear equations column to solver columns map "// &
                                      & "solver colums.",ERR,ERROR,*999)
                                    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                                    IF(ERR/=0) &
                                      & CALL FLAG_ERROR("Could not allocate linear equations column to solver columns map "// &
                                      & "coupling coefficients.",ERR,ERROR,*999)
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=1
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%SOLVER_COLS(1)=solver_global_dof
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%COUPLING_COEFFICIENTS(1)= &
                                      & LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%MATRIX_COEFFICIENT
                                  ENDDO !equations_matrix_idx
                                CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                                  jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(local_dof)
                                  !Allocate the Jacobian to solver map column items.
                                  !No coupling yet so the mapping is 1-1
                                  ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                    & JACOBIAN_COL_TO_SOLVER_COLS_MAP(jacobian_column)%SOLVER_COLS(1),STAT=ERR)
                                  IF(ERR/=0) &
                                    & CALL FLAG_ERROR("Could not allocate Jacobian column to solver columns map solver colums.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                    & JACOBIAN_COL_TO_SOLVER_COLS_MAP(jacobian_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                                  IF(ERR/=0) CALL &
                                    & FLAG_ERROR("Could not allocate Jacobain column to solver columns map coupling coefficients.",&
                                    & ERR,ERROR,*999)
                                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                    & JACOBIAN_COL_TO_SOLVER_COLS_MAP(jacobian_column)%NUMBER_OF_SOLVER_COLS=1
                                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                    & JACOBIAN_COL_TO_SOLVER_COLS_MAP(jacobian_column)%SOLVER_COLS(1)=solver_global_dof
                                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                    & JACOBIAN_COL_TO_SOLVER_COLS_MAP(jacobian_column)%COUPLING_COEFFICIENTS(1)= &
                                    & NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT
                                CASE DEFAULT
                                  LOCAL_ERROR="The equations matrix type of "// &
                                    & TRIM(NUMBER_TO_VSTRING(matrix_type,"*",ERR,ERROR))// &
                                    & " is invalid."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                END SELECT
                              ENDDO !matrix_type_idx
                            ENDIF
                          ELSE
                            IF(rank==myrank) THEN
                              !Set up the equations variables -> solver columns mapping
                              !No coupling yet so the mapping is 1-1
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COLUMN_NUMBERS(local_dof)=0
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COUPLING_COEFFICIENTS( &
                                & local_dof)=0.0_DP
                              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                                & solver_matrix_idx)%VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%ADDITIVE_CONSTANTS( &
                                & local_dof)=0.0_DP
                              DO matrix_type_idx=1,SUB_MATRIX_LIST(0,equations_idx,solver_variable_idx)
                                matrix_type=SUB_MATRIX_LIST(matrix_type_idx,equations_idx,solver_variable_idx)
                                SELECT CASE(matrix_type)
                                CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                                  !Set up the equations columns -> solver columns mapping
                                  DO equations_matrix_idx=1,NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                                    MATRIX_NUMBER=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                    equations_column=DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(local_dof)
                                    !No coupling yet so the mapping is 1-1
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=0
                                  ENDDO !equations_matrix_idx
                                CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                                  !Set up the equations columns -> solver columns mapping
                                  DO equations_matrix_idx=1,NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                    MATRIX_NUMBER=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                    equations_column=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                      & DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(local_dof)
                                    !No coupling yet so the mapping is 1-1
                                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                      & LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR% &
                                      & EQUATIONS_COL_TO_SOLVER_COLS_MAP(equations_column)%NUMBER_OF_SOLVER_COLS=0
                                  ENDDO !equations_matrix_idx
                                CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                                  jacobian_column=NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(local_dof)
                                  !No coupling yet so the mapping is 1-1
                                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP% &
                                    & JACOBIAN_COL_TO_SOLVER_COLS_MAP(jacobian_column)%NUMBER_OF_SOLVER_COLS=0
                                CASE DEFAULT
                                  LOCAL_ERROR="The equations matrix type of "// &
                                    & TRIM(NUMBER_TO_VSTRING(matrix_type,"*",ERR,ERROR))//" is invalid."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                END SELECT
                              ENDDO !matrix_type_idx
                            ENDIF !rank==myrank
                          ENDIF !include_column
                        ENDDO !global_dof

                      CASE(SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION,SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE)
                        
                        !Now handle the interface condition rows and columns
                        
                        interface_condition_idx=SUB_MATRIX_INFORMATION(2,equations_idx,solver_variable_idx)
                        interface_matrix_idx=SUB_MATRIX_INFORMATION(3,equations_idx,solver_variable_idx)
                        
                        !The pointers below have been checked for association above.
                        INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
                        
                        SELECT CASE(INTERFACE_CONDITION%METHOD)
                        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
                          INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING

                          !Loop over the variables
                          
                          LAGRANGE_VARIABLE=>SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES(solver_variable_idx)% &
                            & VARIABLE
                         
                          DO global_dof_idx=1,NUMBER_OF_RANK_COLS
                            global_dof=RANK_GLOBAL_COLS_LIST(1,global_dof_idx)
                            local_dof=RANK_GLOBAL_COLS_LIST(2,global_dof_idx)
                            !dof_type=RANK_GLOBAL_COLS_LIST(3,global_dof_idx)
                            INCLUDE_COLUMN=RANK_GLOBAL_COLS_LIST(3,global_dof_idx)==1
                             
                            IF(INCLUDE_COLUMN) THEN
                              !DOF is not fixed so map the variable/equation dof to a new solver dof
                              
                              IF(dof_type==2) THEN
                                solver_global_dof=DOF_MAP(solver_variable_idx)%PTR(global_dof)
                              ELSE
                                solver_global_dof=solver_global_dof+1
                                DOF_MAP(solver_variable_idx)%PTR(global_dof)=solver_global_dof
                              ENDIF
                              
                              IF(rank==myrank) THEN
                                
                                solver_local_dof=solver_local_dof+1
                                 
                                IF(.NOT.VARIABLE_RANK_PROCESSED(solver_variable_idx,rank)) THEN
                                  
                                  !Set up the column domain mappings.
                                  CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                                    & solver_global_dof),ERR,ERROR,*999)
                                  !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                                  !local map.
                                  !Allocate the global to local map arrays
                                  ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_NUMBER(1),STAT=ERR)
                                  IF(ERR/=0) &
                                    & CALL FLAG_ERROR("Could not allocate column domain global to local map local number.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%DOMAIN_NUMBER(1),STAT=ERR)
                                  IF(ERR/=0) &
                                    & CALL FLAG_ERROR("Could not allocate column domain global to local map domain number.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_TYPE(1),STAT=ERR)
                                  IF(ERR/=0) &
                                    & CALL FLAG_ERROR("Could not allocate column domain global to local map domain number.", &
                                    & ERR,ERROR,*999)
                                  !Set up the global to local mappings
                                  COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%NUMBER_OF_DOMAINS=1
                                  COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_NUMBER(1)= &
                                    & solver_local_dof
                                  COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%DOMAIN_NUMBER(1)=rank
                                  COL_DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(solver_global_dof)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                                  
                                  !Set up the solver column -> equations column mappings. 1-1 as no coupling yet
!!TODO
                                  !Set up the solver dofs -> variable dofs map
                                  !Initialise
                                  CALL SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_INITIALISE(SOLVER_MAPPING% &
                                    & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof),ERR,ERROR,*999)
                                  !Allocate the solver dofs to variable dofs arrays
!!TODO: allow for multiple equations sets in the column
                                  !No coupling so there is only one equations set at the moment
                                  ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                    & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%EQUATIONS_TYPES(1),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps equations types.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                    & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%EQUATIONS_INDICES(1),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps equations indices.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                    & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%VARIABLE(1),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable type.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                    & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%VARIABLE_DOF(1),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable dof.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                    & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%VARIABLE_COEFFICIENT(1),STAT=ERR)
                                  IF(ERR/=0) &
                                    & CALL FLAG_ERROR("Could not allocate solver dof to variable maps variable coefficient.", &
                                    & ERR,ERROR,*999)
                                  ALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
                                    & SOLVER_DOF_TO_VARIABLE_MAPS(solver_local_dof)%ADDITIVE_CONSTANT(1),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dof to variable maps additive constant.", &
                                    & ERR,ERROR,*999)
                                  !Setup
                                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof)%NUMBER_OF_EQUATIONS=1
                                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof)%EQUATIONS_TYPES(1)=SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof)%EQUATIONS_INDICES(1)=interface_condition_idx
                                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof)%VARIABLE(1)%PTR=>LAGRANGE_VARIABLE
                                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof)%VARIABLE_DOF(1)=local_dof
                                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof)%VARIABLE_COEFFICIENT(1)=1.0_DP
                                  SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
                                    & solver_local_dof)%ADDITIVE_CONSTANT(1)=0.0_DP
                                ENDIF
                                IF(equation_type==SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION) THEN
                                  IF(.NOT.VARIABLE_RANK_PROCESSED(solver_variable_idx,rank)) THEN
                                    !Set up the equations variables -> solver columns mapping
                                    !No coupling yet so the mapping is 1-1
                                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP% &
                                      & COLUMN_NUMBERS(local_dof)=solver_global_dof
                                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP% &
                                      & COUPLING_COEFFICIENTS(local_dof)=1.0_DP
                                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP% &
                                      & ADDITIVE_CONSTANTS(local_dof)=0.0_DP
                                    !Set up the equations columns -> solver columns mapping
                                    interface_column=INTERFACE_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(local_dof)
                                    !Allocate the equation to solver map column items.
                                    !No coupling yet so the mapping is 1-1
                                    ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_COL_TO_SOLVER_COLS_MAP( &
                                      & interface_column)%SOLVER_COLS(1),STAT=ERR)
                                    IF(ERR/=0) CALL  FLAG_ERROR("Could not allocate interface column to solver columns map "// &
                                      & "solver colums.",ERR,ERROR,*999)
                                    ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_COL_TO_SOLVER_COLS_MAP( &
                                      & interface_column)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface column to solver columns map "// &
                                      & "coupling coefficients.",ERR,ERROR,*999)
                                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_COL_TO_SOLVER_COLS_MAP( &
                                      & interface_column)%NUMBER_OF_SOLVER_COLS=1
                                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_COL_TO_SOLVER_COLS_MAP( &
                                      & interface_column)%SOLVER_COLS(1)=solver_global_dof
                                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_COL_TO_SOLVER_COLS_MAP( &
                                      & interface_column)%COUPLING_COEFFICIENTS(1)=1.0_DP
                                  ENDIF
                                ELSE
                                  !Set up the equations variables -> solver columns mapping
                                  !No coupling yet so the mapping is 1-1
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS( &
                                    & interface_matrix_idx)%COLUMN_NUMBERS(local_dof)=solver_global_dof
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS( &
                                    & interface_matrix_idx)%COUPLING_COEFFICIENTS(local_dof)=1.0_DP
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS( &
                                    & interface_matrix_idx)%ADDITIVE_CONSTANTS(local_dof)=0.0_DP
                                  !Set up the equations columns -> solver columns mapping
                                  interface_row=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)% &
                                    & VARIABLE_DOF_TO_ROW_MAP(local_dof)
                                  !Allocate the equation to solver map column items.
                                  !No coupling yet so the mapping is 1-1
                                  ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                                    & INTERFACE_ROW_TO_SOLVER_COLS_MAP(interface_row)%SOLVER_COLS(1),STAT=ERR)
                                  IF(ERR/=0) &
                                    & CALL FLAG_ERROR("Could not allocate interface equations row to solver columns map "// &
                                    & "solver colums.",ERR,ERROR,*999)
                                  ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                                    & INTERFACE_ROW_TO_SOLVER_COLS_MAP(interface_row)%COUPLING_COEFFICIENTS(1),STAT=ERR)
                                  IF(ERR/=0) &
                                    & CALL FLAG_ERROR("Could not allocate interface equations row to solver columns map "// &
                                    & "coupling coefficients.",ERR,ERROR,*999)
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                                    & INTERFACE_ROW_TO_SOLVER_COLS_MAP(interface_row)%NUMBER_OF_SOLVER_COLS=1
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                                    & INTERFACE_ROW_TO_SOLVER_COLS_MAP(interface_row)%SOLVER_COLS(1)=solver_global_dof
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                                    & INTERFACE_ROW_TO_SOLVER_COLS_MAP(interface_row)%COUPLING_COEFFICIENTS(1)= &
                                    & INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%MATRIX_COEFFICIENT
                                ENDIF
                              ENDIF
                            ELSE
                              IF(rank==myrank) THEN
                                IF(equation_type==SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION) THEN
                                  !Set up the equations variables -> solver columns mapping
                                  !No coupling yet so the mapping is 1-1
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP%COLUMN_NUMBERS(local_dof)=0
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP%COUPLING_COEFFICIENTS(local_dof)=0.0_DP
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP%ADDITIVE_CONSTANTS(local_dof)=0.0_DP
                                  interface_column=INTERFACE_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(local_dof)
                                  !No coupling yet so the mapping is 1-1
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_COL_TO_SOLVER_COLS_MAP( &
                                    & interface_column)%NUMBER_OF_SOLVER_COLS=0
                                ELSE
                                  !Set up the equations variables -> solver columns mapping
                                  !No coupling yet so the mapping is 1-1
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(interface_matrix_idx)%COLUMN_NUMBERS(local_dof)=0
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(interface_matrix_idx)% &
                                    & COUPLING_COEFFICIENTS(local_dof)=0.0_DP
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(interface_matrix_idx)% &
                                    & ADDITIVE_CONSTANTS(local_dof)=0.0_DP
                                  !Set up the equations columns -> solver columns mapping
                                  interface_row=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)% &
                                    & VARIABLE_DOF_TO_ROW_MAP(local_dof)
                                  !No coupling yet so the mapping is 1-1
                                  SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                    & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                                    & INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                                    & INTERFACE_ROW_TO_SOLVER_COLS_MAP(interface_row)%NUMBER_OF_SOLVER_COLS=0
                                ENDIF
                              ENDIF !rank==myrank
                            ENDIF !include_column
                          ENDDO !global_dof                  
                        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The interface condition method of "// &
                            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                            & " is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE DEFAULT
                        LOCAL_ERROR="The equation type of "//TRIM(NUMBER_TO_VSTRING(equation_type,"*",ERR,ERROR))// &
                          & " is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                      VARIABLE_RANK_PROCESSED(solver_variable_idx,rank)=.TRUE.
                    ENDIF !Number of rank columns > 0
                    IF(ALLOCATED(RANK_GLOBAL_COLS_LIST)) DEALLOCATE(RANK_GLOBAL_COLS_LIST)              
                  ENDDO !equation_idx
                  
                ENDDO !solver_variable_idx
              ENDDO !rank
            ENDDO !dof_type
            
            !Deallocate dof map
            DO solver_variable_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
              DEALLOCATE(DOF_MAP(solver_variable_idx)%PTR)
            ENDDO !solver_variable_idx
            DEALLOCATE(DOF_MAP)
            
            CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(COL_DOMAIN_MAPPING,ERR,ERROR,*999)
            
            IF(ALLOCATED(SUB_MATRIX_INFORMATION)) DEALLOCATE(SUB_MATRIX_INFORMATION)
            IF(ALLOCATED(SUB_MATRIX_LIST)) DEALLOCATE(SUB_MATRIX_LIST)
            IF(ALLOCATED(VARIABLE_RANK_PROCESSED)) DEALLOCATE(VARIABLE_RANK_PROCESSED)
            IF(ALLOCATED(NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS)) DEALLOCATE(NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS)
            IF(ALLOCATED(NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)) DEALLOCATE(NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)
            IF(ALLOCATED(TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)) DEALLOCATE(TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)
            
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
          DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
            INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            SELECT CASE(INTERFACE_CONDITION%METHOD)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
              DO interface_matrix_idx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS(SOLVER_MAPPING% &
                  & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                  & interface_matrix_idx)%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface to solver matrix maps.",ERR,ERROR,*999)
                SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                  & interface_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=0
              ENDDO !interface_matrix_idx
              DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                DO interface_matrix_idx=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                  & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_INTERFACE_MATRICES
                  IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                    & solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR% &
                    & SOLVER_MATRIX_NUMBER==solver_matrix_idx) THEN
                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                      & interface_matrix_idx)%NUMBER_OF_SOLVER_MATRICES=SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                      & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                      & NUMBER_OF_SOLVER_MATRICES+1
                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                      & interface_matrix_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                      & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%&
                      & NUMBER_OF_SOLVER_MATRICES)%PTR=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                      & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                      & interface_matrix_idx)%PTR
                  ENDIF
                ENDDO !interface_matrix_idx
              ENDDO !solver_matrix_idx             
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_PENALTY_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface condition method of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
                        
          ENDDO !interface_condition_idx
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
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of solver matrices = ",SOLVER_MAPPING% &
        & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)               
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of equations sets = ",SOLVER_MAPPING% &
        & NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equations_set_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Region user number        = ",EQUATIONS_SET%REGION%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set user number = ",EQUATIONS_SET%USER_NUMBER, &
          & ERR,ERROR,*999)                
      ENDDO !equations_set_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of interface conditions = ",SOLVER_MAPPING% &
        & NUMBER_OF_INTERFACE_CONDITIONS,ERR,ERROR,*999)
      DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
        INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition index : ",interface_condition_idx, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Parent region user number       = ",INTERFACE_CONDITION%INTERFACE% &
          & PARENT_REGION%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition user number = ",INTERFACE_CONDITION%USER_NUMBER, &
          & ERR,ERROR,*999)                
      ENDDO !equations_set_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equations variables list:",ERR,ERROR,*999)
      DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)        
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",SOLVER_MAPPING%VARIABLES_LIST( &
          & solver_matrix_idx)%NUMBER_OF_VARIABLES,ERR,ERROR,*999)        
        DO variable_idx=1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%NUMBER_OF_VARIABLES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Variable : ",variable_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",SOLVER_MAPPING%VARIABLES_LIST( &
            & solver_matrix_idx)%VARIABLES(variable_idx)%VARIABLE_TYPE,ERR,ERROR,*999)        
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations = ",SOLVER_MAPPING%VARIABLES_LIST( &
            & solver_matrix_idx)%VARIABLES(variable_idx)%NUMBER_OF_EQUATIONS,ERR,ERROR,*999)        
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES( &
            & variable_idx)%NUMBER_OF_EQUATIONS,5,5,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES( &
            & variable_idx)%EQUATION_TYPES,'("        Equation types   :",5(X,I13))','(26X,5(X,I13))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES( &
            & variable_idx)%NUMBER_OF_EQUATIONS,5,5,SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx)%VARIABLES( &
            & variable_idx)%EQUATION_TYPES,'("        Equation indices :",5(X,I13))','(26X,5(X,I13))',ERR,ERROR,*999)
        ENDDO !variable_idx
      ENDDO !solver_matrix_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver row to equations rows mappings:",ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of rows = ",SOLVER_MAPPING%NUMBER_OF_ROWS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global rows = ",SOLVER_MAPPING%NUMBER_OF_GLOBAL_ROWS, &
        & ERR,ERROR,*999)
      IF(DIAGNOSTICS2) THEN
        DO row_idx=1,SOLVER_MAPPING%NUMBER_OF_ROWS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver row : ",row_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations sets mapped to = ",SOLVER_MAPPING% &
            & SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition index          = ",SOLVER_MAPPING% &
            & SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%INTERFACE_CONDITION_INDEX,ERR,ERROR,*999)
          IF(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%INTERFACE_CONDITION_INDEX==0) THEN
            !Row is an equations set row
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)% &
              & NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%EQUATIONS_INDEX, &
              & '("      Equations indices      :",5(X,I13))','(30X,5(X,I13))',ERR,ERROR,*999) 
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)% &
              & NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%ROWCOL_NUMBER, &
              & '("      Equations row numbers  :",5(X,I13))','(30X,5(X,I13))',ERR,ERROR,*999) 
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)% &
              & NUMBER_OF_EQUATIONS_SETS,5,5,SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%COUPLING_COEFFICIENTS, &
              & '("      Coupling coefficients  :",5(X,E13.6))','(30X,5(X,E13.6))',ERR,ERROR,*999)
          ELSE
            !Row is an interface condition row
!!TODO: format better
            CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Interface col numbers : ",SOLVER_MAPPING% &
              & SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%ROWCOL_NUMBER(1),"(I13)",ERR,ERROR,*999) 
            CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficients : ",SOLVER_MAPPING% &
              & SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(row_idx)%COUPLING_COEFFICIENTS(1),"(E13.6)",ERR,ERROR,*999)
          ENDIF
        ENDDO !row_idx
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver column to equations column mappings:",ERR,ERROR,*999)      
      DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",SOLVER_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
        IF(DIAGNOSTICS2) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Solver column to equations set columns mappings:",ERR,ERROR,*999)
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set index : ",equations_set_idx,ERR,ERROR,*999)
            DO column_idx=1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%NUMBER_OF_COLUMNS           
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "          Solver column : ",column_idx,ERR,ERROR,*999)
              IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                & equations_set_idx)%HAVE_DYNAMIC) THEN
               CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES, &
                  & ERR,ERROR,*999)
                IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES>0) THEN
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                    & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                    & column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                    & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                    & column_idx)%EQUATIONS_MATRIX_NUMBERS,'("            Equation matrices numbers :",5(X,I13))', &
                    & '(39X,5(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                    & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                    & column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                    & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                    & column_idx)%EQUATIONS_COL_NUMBERS,'("            Equation column numbers   :",5(X,I13))', &
                    & '(39X,5(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                    & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                    & column_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                    & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS( &
                    & column_idx)%COUPLING_COEFFICIENTS,'("            Coupling coefficients     :",5(X,E13.6))', &
                    & '(39X,5(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ELSE
                 CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & 0_INTG,ERR,ERROR,*999)
              ENDIF
!!TODO what about dynamic nonlinear mappings???
              IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                & equations_set_idx)%HAVE_STATIC) THEN
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "            Number of linear equations matrices mapped to  = ", &
                  & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES, &
                  & ERR,ERROR,*999)
                IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
                  !CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                  !  & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                  !  & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                  !  & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                  !  & column_idx)%EQUATIONS_MATRIX_NUMBERS,'("            Equation matrices numbers :",5(X,I13))', &
                  !  & '(36X,5(X,I13))',ERR,ERROR,*999)
                  !CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                  !  & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                  !  & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                  !  & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                  !  & column_idx)%EQUATIONS_COL_NUMBERS,'("            Equation column numbers   :",5(X,I13))', &
                  !  & '(36X,5(X,I13))',ERR,ERROR,*999)
                  !CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                  !  & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                  !  & column_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,5,5,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
                  !  & solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS( &
                  !  & column_idx)%COUPLING_COEFFICIENTS,'("            Coupling coefficients     :",5(X,E13.6))', &
                  !  & '(36X,5(X,E13.6))',ERR,ERROR,*999)
                ENDIF
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian column number     = ", &
                  & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%JACOBIAN_COL_NUMBER,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian coupling coeff    = ", &
                  & SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_COL_TO_EQUATIONS_SET_MAPS( &
                  & equations_set_idx)%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(column_idx)%JACOBIAN_COUPLING_COEFFICIENT,ERR,ERROR,*999)
              ELSE
                 CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of static equations matrices mapped to  = ", &
                  & 0_INTG,ERR,ERROR,*999)
              ENDIF
            ENDDO !column_idx
          ENDDO !equations_set_idx
        ENDIF
      ENDDO !solver_matrix_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver DOF to field DOFs mappings:",ERR,ERROR,*999)
      DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of DOFs = ",SOLVER_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",SOLVER_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of global DOFs = ",SOLVER_MAPPING% &
          & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%NUMBER_OF_GLOBAL_DOFS,ERR,ERROR,*999)
        ALLOCATE(VARIABLE_TYPES(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable types.",ERR,ERROR,*999)
        DO dof_idx=1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%TOTAL_NUMBER_OF_DOFS     
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Solver local DOF : ",dof_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations mapped to     = ",SOLVER_MAPPING% &
            & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
            & NUMBER_OF_EQUATIONS,ERR,ERROR,*999)
          IF(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
            & NUMBER_OF_EQUATIONS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%EQUATIONS_INDICES, &
              & '("        Equations types       :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%EQUATIONS_INDICES, &
              & '("        Equations indices     :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
            DO variable_idx=1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS( &
              & dof_idx)%NUMBER_OF_EQUATIONS
              VARIABLE_TYPES(variable_idx)=SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)% &
              & SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%VARIABLE(variable_idx)%PTR%VARIABLE_TYPE
            ENDDO !variable_idx
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS,5,5,VARIABLE_TYPES, &
              & '("        Variable types        :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%VARIABLE_DOF, &
              & '("        Variable DOFs         :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
              & VARIABLE_COEFFICIENT, & 
              & '("        Variable coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
              & solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)%NUMBER_OF_EQUATIONS,5,5,SOLVER_MAPPING% &
              & SOLVER_COL_TO_EQUATIONS_COLS_MAP(solver_matrix_idx)%SOLVER_DOF_TO_VARIABLE_MAPS(dof_idx)% &
              & ADDITIVE_CONSTANT, &
              & '("        Additive constants    :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !dof_idx
        IF(ALLOCATED(VARIABLE_TYPES)) DEALLOCATE(VARIABLE_TYPES)
      ENDDO !solver_matrix_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets to solver mappings:",ERR,ERROR,*999)
      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
        LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
        RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
        SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equations_set_idx,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Equations sets rows to solver rows mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "        Number of equations set rows = ",EQUATIONS_MAPPING% &
         & TOTAL_NUMBER_OF_ROWS,ERR,ERROR,*999)
        DO row_idx=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set row : ",row_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to   = ",SOLVER_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS, &
            & ERR,ERROR,*999)
          IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)% &
            & NUMBER_OF_SOLVER_ROWS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS,5,5,SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%SOLVER_ROWS, &
              & '("          Solver row numbers    :",5(X,I13))','(33X,5(X,I13))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%NUMBER_OF_SOLVER_ROWS,5,5,SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(row_idx)%COUPLING_COEFFICIENTS, &
              & '("          Coupling coefficients :",5(X,E13.6))','(33X,5(X,E13.6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !row_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix indexing:",ERR,ERROR,*999)
        DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"          Interface conditions affecting:", &
              & ERR,ERROR,*999)
           CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of interface conditions = ",SOLVER_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%NUMBER_OF_INTERFACE_CONDITIONS,ERR,ERROR,*999)
          DO interface_condition_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%NUMBER_OF_INTERFACE_CONDITIONS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Interface condition : ",interface_condition_idx, &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Interface condition index = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE( &
              & interface_condition_idx)%INTERFACE_CONDITION_INDEX,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Interface matrix number   = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE( &
              & interface_condition_idx)%INTERFACE_MATRIX_NUMBER,ERR,ERROR,*999)
          ENDDO !interface_condition_idx                                 
          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
           CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"          Dynamic equations matrix columns to solver matrix columns:", &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
              & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,ERR,ERROR,*999)
            DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                & equations_matrix_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix index : ",equations_matrix_idx, &
                & ERR,ERROR,*999)
              equations_matrix=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX_NUMBER
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Equations matrix number = ",equations_matrix, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Solver matrix number    = ",EQUATIONS_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
              DO column_idx=1,DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix)%NUMBER_OF_COLUMNS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Equations matrix column : ",column_idx, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                  Number of solver columns mapped to = ", &
                  & EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS,ERR,ERROR,*999)
                IF(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS>0) THEN
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP(column_idx)% &
                    & SOLVER_COLS,'("                  Solver columns         :",5(X,I13))','(42X,5(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP(column_idx)% &
                    & COUPLING_COEFFICIENTS,'("                  Coupling coefficients  :",5(X,E13.6))','(42X,5(X,E13.6))', &
                    & ERR,ERROR,*999)
                ENDIF
              ENDDO !column_idx
            ENDDO !equations_matrix_idx
          ELSE
            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"          Linear equations matrix columns to solver matrix columns:", &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of liner equations matrices = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*999)
              DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                  & equations_matrix_idx)%PTR
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix index : ",equations_matrix_idx, &
                  & ERR,ERROR,*999)
                equations_matrix=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX_NUMBER
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix number = ",equations_matrix, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Solver matrix number    = ", &
                  & EQUATIONS_TO_SOLVER_MAP%SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
                DO column_idx=1,LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix)%NUMBER_OF_COLUMNS
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix column : ",column_idx, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Number of solver columns mapped to = ", &
                    & EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS,ERR,ERROR,*999)
                  IF(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS>0) THEN
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP( &
                      & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP( &
                      & column_idx)%SOLVER_COLS,'("                Solver columns         :",5(X,I13))','(40X,5(X,I13))', &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP( &
                      & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP( &
                      & column_idx)%COUPLING_COEFFICIENTS, &
                      & '("                Coupling coefficients  :",5(X,E13.6))','(40X,5(X,E13.6))',ERR,ERROR,*999)
                  ENDIF
                ENDDO !column_idx
              ENDDO !equations_matrix_idx
            ENDIF
            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
               CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"          Jacobian equations matrix columns to solver matrix columns:", &
                & ERR,ERROR,*999)
              JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",JACOBIAN_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
              DO column_idx=1,NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix column : ",column_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                  & JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS,ERR,ERROR,*999)
                IF(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP(column_idx)%NUMBER_OF_SOLVER_COLS>0) THEN
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP(column_idx)% &
                    & SOLVER_COLS,'("              Solver columns         :",5(X,I13))','(38X,5(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP( &
                    & column_idx)%NUMBER_OF_SOLVER_COLS,5,5,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP(column_idx)% &
                    & COUPLING_COEFFICIENTS,'("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))', &
                    & ERR,ERROR,*999)
                ENDIF
              ENDDO !column_idx
            ENDIF
          ENDIF
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"          Variable DOFs to solver matrix DOFs:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of variables = ",SOLVER_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
            & NUMBER_OF_VARIABLES,ERR,ERROR,*999) 
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_VARIABLES,5,5,SOLVER_MAPPING% &
            & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
            & VARIABLE_TYPES,'("            Variable types :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
          DO variable_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_VARIABLES
            DEPENDENT_VARIABLE=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%VARIABLES(variable_idx)%PTR
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Variable index : ",variable_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Number of variable DOFs = ",DEPENDENT_VARIABLE% &
              & NUMBER_OF_DOFS,ERR,ERROR,*999)
            DO local_dof=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Variable DOF : ",local_dof,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Solver column number = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COLUMN_NUMBERS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COUPLING_COEFFICIENTS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
                & VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%ADDITIVE_CONSTANTS(local_dof),ERR,ERROR,*999)              
            ENDDO !local_dof
          ENDDO !variable_idx
        ENDDO !solver_matrix_idx
        IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Dynamic equations matrix indexing:",ERR,ERROR,*999)
          DO equations_matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equations_matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)% &
              & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
            DO solver_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES
              EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                & solver_matrix_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solver_matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ",EQUATIONS_TO_SOLVER_MAP% &
                & EQUATIONS_MATRIX_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",EQUATIONS_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)            
            ENDDO !solver_matrix_idx
          ENDDO !equations_matrix_idx
        ELSE
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Linear equations matrix indexing:",ERR,ERROR,*999)
            DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equations_matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",SOLVER_MAPPING% &
                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)% &
                & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
              DO solver_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%NUMBER_OF_SOLVER_MATRICES
                EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                  & solver_matrix_idx)%PTR
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solver_matrix_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ",EQUATIONS_TO_SOLVER_MAP% &
                  & EQUATIONS_MATRIX_NUMBER,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",EQUATIONS_TO_SOLVER_MAP% &
                  & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)            
              ENDDO !solver_matrix_idx
            ENDDO !equations_matrix_idx
          ENDIF
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian matrix indexing:",ERR,ERROR,*999)
            JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM%JACOBIAN_TO_SOLVER_MATRIX_MAP
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",JACOBIAN_TO_SOLVER_MAP% &
              & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDDO !equations_set_idx
      IF(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS>0) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions to solver mappings:",ERR,ERROR,*999)
        DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
          INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Interface to equations sets mapping:",ERR,ERROR,*999)
          DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations sets = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
            DO equations_set_idx=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & NUMBER_OF_EQUATIONS_SETS
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set : ",equations_set_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Equations set index     = ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS( &
                & equations_set_idx)%EQUATIONS_SET_INDEX,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS( &
                & equations_set_idx)%INTERFACE_MATRIX_INDEX,ERR,ERROR,*999)
            ENDDO !equations_set_idx
          ENDDO !solver_matrix_idx
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition rows to solver rows mappings:",ERR,ERROR,*999)
          DO interface_matrix_idx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix idx : ",interface_matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of interface rows = ",INTERFACE_MAPPING% &
              & INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%TOTAL_NUMBER_OF_ROWS,ERR,ERROR,*999)
            DO row_idx=1,INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%TOTAL_NUMBER_OF_ROWS
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Interface row : ",row_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to = ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                & interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP(row_idx)%NUMBER_OF_SOLVER_ROWS,ERR,ERROR,*999)
              IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                & interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP(row_idx)%NUMBER_OF_SOLVER_ROWS>0) THEN
                CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Solver row numbers    : ",SOLVER_MAPPING% &
                  & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                  & interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP(row_idx)%SOLVER_ROW,"(I13)",ERR,ERROR,*999)
                CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Coupling coefficients : ",SOLVER_MAPPING% &
                  & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
                  & interface_matrix_idx)%INTERFACE_ROW_TO_SOLVER_ROWS_MAP(row_idx)%COUPLING_COEFFICIENT,"(E13.6)",ERR,ERROR,*999)
              ENDIF
            ENDDO !row_idx
          ENDDO !interface_matrix_idx
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition column to solver rows mappings:", &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of interface condition columns = ",INTERFACE_MAPPING% &
            & TOTAL_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
          DO column_idx=1,INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition column : ",column_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows mapped to = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(column_idx)% &
              & NUMBER_OF_SOLVER_ROWS,ERR,ERROR,*999)
            IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(column_idx)%NUMBER_OF_SOLVER_ROWS>0) THEN
              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver row number    : ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(column_idx)% &
                & SOLVER_ROW, "(I13)",ERR,ERROR,*999)
              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficient : ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(column_idx)% &
                & COUPLING_COEFFICIENT, "(E13.6)",ERR,ERROR,*999)
            ENDIF
          ENDDO !column_idx
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix indexing:",ERR,ERROR,*999)
          DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)        
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"        Interface equations matrix rows to solver matrix columns:", &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of interface matrices = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
              & NUMBER_OF_INTERFACE_MATRICES,ERR,ERROR,*999)
            DO interface_matrix_idx=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_INTERFACE_MATRICES
              INTERFACE_TO_SOLVER_MAP=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                & interface_matrix_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix index : ",interface_matrix_idx, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix number = ",INTERFACE_TO_SOLVER_MAP% &
                & INTERFACE_MATRIX_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",INTERFACE_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)
              DO row_idx=1,INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%NUMBER_OF_ROWS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix row : ",row_idx, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                  & INTERFACE_TO_SOLVER_MAP%INTERFACE_ROW_TO_SOLVER_COLS_MAP(row_idx)%NUMBER_OF_SOLVER_COLS,ERR,ERROR,*999)
                IF(INTERFACE_TO_SOLVER_MAP%INTERFACE_ROW_TO_SOLVER_COLS_MAP(row_idx)%NUMBER_OF_SOLVER_COLS>0) THEN
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,INTERFACE_TO_SOLVER_MAP%INTERFACE_ROW_TO_SOLVER_COLS_MAP( &
                    & row_idx)%NUMBER_OF_SOLVER_COLS,5,5,INTERFACE_TO_SOLVER_MAP%INTERFACE_ROW_TO_SOLVER_COLS_MAP(row_idx)% &
                    & SOLVER_COLS,'("              Solver columns         :",5(X,I13))','(38X,5(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,INTERFACE_TO_SOLVER_MAP%INTERFACE_ROW_TO_SOLVER_COLS_MAP( &
                    & row_idx)%NUMBER_OF_SOLVER_COLS,5,5,INTERFACE_TO_SOLVER_MAP%INTERFACE_ROW_TO_SOLVER_COLS_MAP(row_idx)% &
                    & COUPLING_COEFFICIENTS,'("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))', &
                    & ERR,ERROR,*999)
                ENDIF
              ENDDO !row_idx
            ENDDO !interface_matrix_idx
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"        Variable dofs to solver matrix dofs:",ERR,ERROR,*999)
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"          Lagrange variables:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Lagrange variable type = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
              & LAGRANGE_VARIABLE_TYPE,ERR,ERROR,*999)
            LAGRANGE_VARIABLE=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LAGRANGE_VARIABLE
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of Lagrange variable dofs = ",LAGRANGE_VARIABLE% &
              & NUMBER_OF_DOFS,ERR,ERROR,*999)
            DO local_dof=1,LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Variable dof : ",local_dof,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Solver column number = ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                & solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP%COLUMN_NUMBERS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Coupling coefficient = ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                & solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP%COUPLING_COEFFICIENTS(local_dof),ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Additive constant    = ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                & solver_matrix_idx)%LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP%ADDITIVE_CONSTANTS(local_dof),ERR,ERROR,*999)              
            ENDDO !local_dof
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"          Dependent variables:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dependent variables = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)% &
              & NUMBER_OF_DEPENDENT_VARIABLES,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
              & interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DEPENDENT_VARIABLES, &
              & 5,5,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
              & solver_matrix_idx)%DEPENDENT_VARIABLE_TYPES,'("            Dependent variable types :",5(X,I13))', &
              & '(38X,5(X,I13))',ERR,ERROR,*999) 
            DO variable_idx=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DEPENDENT_VARIABLES
              DEPENDENT_VARIABLE=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DEPENDENT_VARIABLES(variable_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Dependent variable index : ",variable_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Number of dependent variable dofs = ", &
                & DEPENDENT_VARIABLE%NUMBER_OF_DOFS,ERR,ERROR,*999)
              DO local_dof=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"              Variable dof : ",local_dof,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Solver column number = ",SOLVER_MAPPING% &
                  & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COLUMN_NUMBERS(local_dof), &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ",SOLVER_MAPPING% &
                  & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%COUPLING_COEFFICIENTS(local_dof), &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ",SOLVER_MAPPING% &
                  & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM( &
                  & solver_matrix_idx)%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(variable_idx)%ADDITIVE_CONSTANTS(local_dof), &
                  & ERR,ERROR,*999)              
              ENDDO !local_dof
            ENDDO !variable_idx
          ENDDO !solver_matrix_idx        
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Interface equations matrix indexing:",ERR,ERROR,*999)
          DO interface_matrix_idx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix : ",interface_matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver matrices = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM( &
              & interface_matrix_idx)%NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
            DO solver_matrix_idx=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%NUMBER_OF_SOLVER_MATRICES
              INTERFACE_TO_SOLVER_MAP=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)%INTERFACE_TO_SOLVER_MATRIX_MAPS( &
                & solver_matrix_idx)%PTR
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix index : ",solver_matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ",INTERFACE_TO_SOLVER_MAP% &
                & INTERFACE_MATRIX_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix number    = ",INTERFACE_TO_SOLVER_MAP% &
                & SOLVER_MATRIX_NUMBER,ERR,ERROR,*999)            
            ENDDO !solver_matrix_idx
          ENDDO !equations_matrix_idx
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Interface column to solver rows mapping:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",INTERFACE_MAPPING%NUMBER_OF_COLUMNS, &
            & ERR,ERROR,*999)            
          DO column_idx=1,INTERFACE_MAPPING%NUMBER_OF_COLUMNS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Column : ",column_idx, ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver rows = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(column_idx)% &
              & NUMBER_OF_SOLVER_ROWS,ERR,ERROR,*999)
            IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS( &
              & column_idx)%NUMBER_OF_SOLVER_ROWS>0) THEN
              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver row             : ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS( &
                & column_idx)%SOLVER_ROW,"(I13)",ERR,ERROR,*999)
              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficients  : ",SOLVER_MAPPING% &
                & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS( &
                & column_idx)%COUPLING_COEFFICIENT,"(E13.6)", &
                & ERR,ERROR,*999)
            ENDIF
          ENDDO !column_idx
        ENDDO !interface_condition_idx
      ENDIF
    ENDIF

    CALL EXITS("SOLVER_MAPPING_CALCULATE")
    RETURN
999 IF(ALLOCATED(SUB_MATRIX_INFORMATION)) DEALLOCATE(SUB_MATRIX_INFORMATION)
    IF(ALLOCATED(SUB_MATRIX_LIST)) DEALLOCATE(SUB_MATRIX_LIST)
    IF(ALLOCATED(VARIABLE_RANK_PROCESSED)) DEALLOCATE(VARIABLE_RANK_PROCESSED)
    IF(ALLOCATED(NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS)) DEALLOCATE(NUMBER_OF_VARIABLE_GLOBAL_SOLVER_DOFS)
    IF(ALLOCATED(NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)) DEALLOCATE(NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)
    IF(ALLOCATED(TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)) DEALLOCATE(TOTAL_NUMBER_OF_VARIABLE_LOCAL_SOLVER_DOFS)    
    CALL ERRORS("SOLVER_MAPPING_CALCULATE",ERR,ERROR)
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
    INTEGER(INTG) :: equations_set_idx,solver_matrix_idx

    CALL ENTERS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ASSOCIATED(CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST)) THEN
        DO solver_matrix_idx=1,SIZE(CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST,1)
          IF(ASSOCIATED(CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR)) & 
            & CALL LIST_DESTROY(CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !solver_matrix_idx
        DEALLOCATE(CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST)
      ENDIF
      IF(ALLOCATED(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE)
      IF(ASSOCIATED(CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST)) THEN
        DO solver_matrix_idx=1,SIZE(CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST,1)
          IF(ASSOCIATED(CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR)) & 
            & CALL LIST_DESTROY(CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !solver_matrix_idx
        DEALLOCATE(CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST)
      ENDIF
      IF(ASSOCIATED(CREATE_VALUES_CACHE%INTERFACE_INDICES)) THEN
        DO equations_set_idx=1,SIZE(CREATE_VALUES_CACHE%INTERFACE_INDICES,1)
          IF(ASSOCIATED(CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR)) &
            & CALL LIST_DESTROY(CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR, &
            & ERR,ERROR,*999)
        ENDDO !equaitons_set_idx
        DEALLOCATE(CREATE_VALUES_CACHE%INTERFACE_INDICES)
      ENDIF
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
    INTEGER(INTG) :: DUMMY_ERR,equations_set_idx,solver_matrix_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Solver mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache.",ERR,ERROR,*999)
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache equations variable list.", &
          & ERR,ERROR,*999)
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
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache interface variable list.", &
          & ERR,ERROR,*999)
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver mapping create values cache interface condition indices.", &
          & ERR,ERROR,*999)
        DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
          NULLIFY(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR)
          CALL LIST_CREATE_START(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR,ERR,ERROR,*999)
          CALL LIST_DATA_TYPE_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR, &
            & LIST_INTG_TYPE,ERR,ERROR,*999)
          CALL LIST_DATA_DIMENSION_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR, &
            & 2,ERR,ERROR,*999)
          CALL LIST_CREATE_FINISH(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR,ERR,ERROR,*999)
          NULLIFY(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR)
          CALL LIST_CREATE_START(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR,ERR,ERROR,*999)
          CALL LIST_DATA_TYPE_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR, &
            & LIST_INTG_TYPE,ERR,ERROR,*999)
          CALL LIST_DATA_DIMENSION_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR, &
            & 2,ERR,ERROR,*999)
          CALL LIST_CREATE_FINISH(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !solver_idx
        DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
          NULLIFY(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR)
          CALL LIST_CREATE_START(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR,ERR,ERROR,*999)
          CALL LIST_DATA_TYPE_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR,LIST_INTG_TYPE, &
            & ERR,ERROR,*999)
          CALL LIST_DATA_DIMENSION_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR,2, &
            & ERR,ERROR,*999)
          CALL LIST_KEY_DIMENSION_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR,1, &
            & ERR,ERROR,*999)
          CALL LIST_CREATE_FINISH(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)%PTR, &
            & ERR,ERROR,*999)            
        ENDDO !equations_set_idx
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

  !>Adds a variable type from an equations set dependent field to the list of variables for a particular solver matrix of a solver mapping.
  SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD(SOLVER_MAPPING,solver_matrix_idx,equations_set_idx, &
    & variable_type,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to add to the var list for.
    INTEGER(INTG), INTENT(IN) :: solver_matrix_idx !<The solver matrix index of the variable list
    INTEGER(INTG), INTENT(IN) :: equations_set_idx !<The equations set index of the variable to add
    INTEGER(INTG), INTENT(IN) :: variable_type !<The variable type to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx2,NUMBER_OF_VARIABLES,variable_idx,VARIABLE_ITEM(2)
    LOGICAL :: VARIABLE_FOUND
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET,VAR_EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,VAR_DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(equations_set_idx>0.AND.equations_set_idx<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS) THEN
        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
            IF(variable_type/=0) THEN
              VARIABLE_FOUND=.FALSE.
              CALL LIST_NUMBER_OF_ITEMS_GET(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR, &
                & NUMBER_OF_VARIABLES,ERR,ERROR,*999)
              DO variable_idx=1,NUMBER_OF_VARIABLES
                CALL LIST_ITEM_GET(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR, &
                  & variable_idx,VARIABLE_ITEM,ERR,ERROR,*999)
                equations_set_idx2=VARIABLE_ITEM(1)
                VAR_EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx2)%PTR
                IF(ASSOCIATED(VAR_EQUATIONS_SET)) THEN
                  VAR_DEPENDENT_FIELD=>VAR_EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(VAR_DEPENDENT_FIELD)) THEN
                    IF(ASSOCIATED(DEPENDENT_FIELD,VAR_DEPENDENT_FIELD)) THEN
                      IF(variable_type==VARIABLE_ITEM(2)) VARIABLE_FOUND=.TRUE.                               
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Variable dependent field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Variable equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !variable_idx
              IF(.NOT.VARIABLE_FOUND) THEN
                VARIABLE_ITEM(1)=equations_set_idx
                VARIABLE_ITEM(2)=variable_type
                CALL LIST_ITEM_ADD(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_VARIABLE_LIST(solver_matrix_idx)%PTR, &
                  & VARIABLE_ITEM,ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified equations set index of "//TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
          & " is invalid. The index must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,"*",ERR,ERROR))//"."        
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD

  !
  !================================================================================================================================
  !

  !>Adds a variable type from an interface condition Lagrange field to the list of variables for a particular solver matrix of a solver mapping.
  SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_INTERF_VAR_LIST_ADD(SOLVER_MAPPING,solver_matrix_idx,interface_condition_idx, &
    & variable_type,ERR,ERROR,*)
    
    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to add to the var list for.
    INTEGER(INTG), INTENT(IN) :: solver_matrix_idx !<The solver matrix index of the variable list
    INTEGER(INTG), INTENT(IN) :: interface_condition_idx !<The interface condition index of the variable to add
    INTEGER(INTG), INTENT(IN) :: variable_type !<The variable type to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_condition_idx2,NUMBER_OF_VARIABLES,variable_idx,VARIABLE_ITEM(2)
    LOGICAL :: VARIABLE_FOUND
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION,VAR_INTERFACE_CONDITION
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD,VAR_LAGRANGE_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INTERF_VAR_LIST_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(interface_condition_idx>0.AND.interface_condition_idx<=SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS) THEN
        INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
            IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
              LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
              IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                IF(variable_type/=0) THEN
                  VARIABLE_FOUND=.FALSE.
                  CALL LIST_NUMBER_OF_ITEMS_GET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR, &
                    & NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                  DO variable_idx=1,NUMBER_OF_VARIABLES
                    CALL LIST_ITEM_GET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR, &
                      & variable_idx,VARIABLE_ITEM,ERR,ERROR,*999)
                    interface_condition_idx2=VARIABLE_ITEM(1)
                    VAR_INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx2)%PTR
                    IF(ASSOCIATED(VAR_INTERFACE_CONDITION)) THEN
                      SELECT CASE(VAR_INTERFACE_CONDITION%METHOD)
                      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                        IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
                          VAR_LAGRANGE_FIELD=>VAR_INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
                          IF(ASSOCIATED(VAR_LAGRANGE_FIELD)) THEN
                            IF(ASSOCIATED(LAGRANGE_FIELD,VAR_LAGRANGE_FIELD)) THEN
                              IF(variable_type==VARIABLE_ITEM(2)) VARIABLE_FOUND=.TRUE.                               
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Variable Lagrange field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Variable interface Lagrange is not associated.",ERR,ERROR,*999)
                        ENDIF
                      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The interface condition method of "// &
                          & TRIM(NUMBER_TO_VSTRING(VAR_INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                          & " is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Variable equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !variable_idx
                  IF(.NOT.VARIABLE_FOUND) THEN
                    VARIABLE_ITEM(1)=interface_condition_idx
                    VARIABLE_ITEM(2)=variable_type
                    CALL LIST_ITEM_ADD(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_VARIABLE_LIST(solver_matrix_idx)%PTR, &
                      & VARIABLE_ITEM,ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                CALL FLAG_ERROR("Lagrange field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface condition Lagrange is not asssociated.",ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified interface condition index of "// &
          & TRIM(NUMBER_TO_VSTRING(interface_condition_idx,"*",ERR,ERROR))// &
          & " is invalid. The index must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS,"*",ERR,ERROR))//"."        
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INTERF_VAR_LIST_ADD")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INTERF_VAR_LIST_ADD",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_INTERF_VAR_LIST_ADD")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_INTERF_VAR_LIST_ADD

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
  SUBROUTINE SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(EQUATIONS_COL_TO_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE) :: EQUATIONS_COL_TO_SOLVER_COLS_MAP !<The equations col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_COL_TO_SOLVER_COLS_MAP%SOLVER_COLS)) &
      & DEALLOCATE(EQUATIONS_COL_TO_SOLVER_COLS_MAP%SOLVER_COLS)
    IF(ALLOCATED(EQUATIONS_COL_TO_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(EQUATIONS_COL_TO_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)
        
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
  SUBROUTINE SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE(EQUATIONS_COL_TO_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE) :: EQUATIONS_COL_TO_SOLVER_COLS_MAP !<The equations column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_COL_TO_SOLVER_COLS_MAP%NUMBER_OF_SOLVER_COLS=0
    
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
    INTEGER(INTG) :: equations_set_idx,matrix_idx,solver_matrix_idx,variable_idx,variable_type
    INTEGER(INTG), ALLOCATABLE :: OLD_DYNAMIC_VARIABLE_TYPE(:),OLD_MATRIX_VARIABLE_TYPES(:,:,:),OLD_RHS_VARIABLE_TYPE(:), &
      & OLD_RESIDUAL_VARIABLE_TYPE(:),OLD_SOURCE_VARIABLE_TYPE(:)
    LOGICAL :: MATRIX_DONE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_SET_PTR_TYPE), ALLOCATABLE :: OLD_EQUATIONS_SETS(:)
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(LIST_PTR_TYPE), POINTER :: NEW_INTERFACE_INDICES(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_INTERFACE_INDICES)
    
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
                      ALLOCATE(NEW_INTERFACE_INDICES(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interface indices.",ERR,ERROR,*999)
                      OLD_DYNAMIC_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE
                      OLD_MATRIX_VARIABLE_TYPES=SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
                      OLD_RESIDUAL_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE
                      OLD_RHS_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE
                      OLD_SOURCE_VARIABLE_TYPE=SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE
                      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        NEW_INTERFACE_INDICES(equations_set_idx)%PTR=>SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES( &
                          & equations_set_idx)%PTR
                      ENDDO !equations_sets
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
                      IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES)
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
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES=>NEW_INTERFACE_INDICES
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1,:)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE(SOLVER_MAPPING% &
                        & NUMBER_OF_EQUATIONS_SETS+1)=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES=>NEW_INTERFACE_INDICES
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
                      IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES)) &
                        & DEALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES)
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
                      ALLOCATE(NEW_INTERFACE_INDICES(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interface indices.",ERR,ERROR,*999)
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=0
                      SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES=>NEW_INTERFACE_INDICES
                    ELSE
                      CALL FLAG_ERROR("The number of equations sets is < 0.",ERR,ERROR,*999)
                    ENDIF
                    NULLIFY(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+1)%PTR)
                    CALL LIST_CREATE_START(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(SOLVER_MAPPING% &
                      & NUMBER_OF_EQUATIONS_SETS+1)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(SOLVER_MAPPING% &
                      & NUMBER_OF_EQUATIONS_SETS+1)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_DATA_DIMENSION_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(SOLVER_MAPPING% &
                      & NUMBER_OF_EQUATIONS_SETS+1)%PTR,2,ERR,ERROR,*999)
                    CALL LIST_KEY_DIMENSION_SET(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(SOLVER_MAPPING% &
                      & NUMBER_OF_EQUATIONS_SETS+1)%PTR,1,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(SOLVER_MAPPING% &
                      & NUMBER_OF_EQUATIONS_SETS+1)%PTR,ERR,ERROR,*999)
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
                    
                    !Add the variables to the list of variables
                    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(EQUATIONS_SET_INDEX)
                    CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD(SOLVER_MAPPING,1,EQUATIONS_SET_INDEX,variable_type, &
                      & ERR,ERROR,*999)
                    DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                      DO matrix_idx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                        & MATRIX_VARIABLE_TYPES(0,EQUATIONS_SET_INDEX,solver_matrix_idx)
                        variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                          & MATRIX_VARIABLE_TYPES(matrix_idx,EQUATIONS_SET_INDEX,solver_matrix_idx)
                        CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD(SOLVER_MAPPING,1,EQUATIONS_SET_INDEX, &
                          & variable_type,ERR,ERROR,*999)
                      ENDDO !matrix_idx
                    ENDDO !solver_matrix_idx
                    variable_type=SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(EQUATIONS_SET_INDEX)
                    CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_EQN_VAR_LIST_ADD(SOLVER_MAPPING,1,EQUATIONS_SET_INDEX,variable_type, &
                      & ERR,ERROR,*999)
                    
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
    INTEGER(INTG) :: equations_matrix_idx,interface_condition_idx,row_idx,solver_matrix_idx
    
    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE)) THEN
      DO interface_condition_idx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE,1)
        CALL SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_FINALISE(EQUATIONS_SET_TO_SOLVER_MAP% &
          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(interface_condition_idx),ERR,ERROR,*999)
      ENDDO !interface_condition_idx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE)
    ENDIF
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
    EQUATIONS_SET_TO_SOLVER_MAP%NUMBER_OF_INTERFACE_CONDITIONS=0
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
      IF(ALLOCATED(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP)) THEN
        DO column_idx=1,SIZE(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP,1)
          CALL SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(EQUATIONS_TO_SOLVER_MAP% &
            & EQUATIONS_COL_TO_SOLVER_COLS_MAP(column_idx),ERR,ERROR,*999)
        ENDDO !column_idx
        DEALLOCATE(EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_TO_SOLVER_COLS_MAP)
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

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_FINALISE(EQUATIONS_TO_SOLVER_INTERFACE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE_TYPE) :: EQUATIONS_TO_SOLVER_INTERFACE_MAP !<The equations set to solver map interface to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_FINALISE",ERR,ERROR,*999)

    EQUATIONS_TO_SOLVER_INTERFACE_MAP%INTERFACE_CONDITION_INDEX=0
    NULLIFY(EQUATIONS_TO_SOLVER_INTERFACE_MAP%INTERFACE_CONDITION)
    EQUATIONS_TO_SOLVER_INTERFACE_MAP%INTERFACE_MATRIX_NUMBER=0
        
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver matrix interface map.
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_INITIALISE(EQUATIONS_TO_SOLVER_INTERFACE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE_TYPE) :: EQUATIONS_TO_SOLVER_INTERFACE_MAP !<The equations set to solver map interface to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_TO_SOLVER_INTERFACE_MAP%INTERFACE_CONDITION_INDEX=0
    NULLIFY(EQUATIONS_TO_SOLVER_INTERFACE_MAP%INTERFACE_CONDITION)
    EQUATIONS_TO_SOLVER_INTERFACE_MAP%INTERFACE_MATRIX_NUMBER=0
        
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_TO_SOLVER_INTERFACE_INITIALISE

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
    INTEGER(INTG) :: equations_set_idx,interface_condition_idx,row_idx,solver_matrix_idx

    CALL ENTERS("SOLVER_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ALLOCATED(SOLVER_MAPPING%VARIABLES_LIST)) THEN
        DO solver_matrix_idx=1,SIZE(SOLVER_MAPPING%VARIABLES_LIST,1)
          CALL SOLVER_MAPPING_VARIABLES_FINALISE(SOLVER_MAPPING%VARIABLES_LIST(solver_matrix_idx),ERR,ERROR,*999)
        ENDDO ! solver_matrix_idx
        DEALLOCATE(SOLVER_MAPPING%VARIABLES_LIST)
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%EQUATIONS_SETS)) DEALLOCATE(SOLVER_MAPPING%EQUATIONS_SETS)        
      IF(ALLOCATED(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)) THEN
        DO equations_set_idx=1,SIZE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP,1)
          CALL SOLVER_MAPPING_EQUATIONS_SET_TO_SOLVER_MAP_FINALISE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
            & equations_set_idx),ERR,ERROR,*999)
        ENDDO !equations_set_idx
        DEALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%INTERFACE_CONDITIONS)) DEALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITIONS)
      IF(ALLOCATED(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP)) THEN
        DO interface_condition_idx=1,SIZE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP,1)
          CALL SOLVER_MAPPING_INTERF_CONDITION_TO_SOLVER_MAP_FINALISE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
            & interface_condition_idx),ERR,ERROR,*999)
        ENDDO !interface_condition_idx
        DEALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP)
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP)) THEN
        DO solver_matrix_idx=1,SIZE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP,1)
          CALL SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_FINALISE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP( &
            & solver_matrix_idx),ERR,ERROR,*999)
        ENDDO !solver_matrix_idx
        DEALLOCATE(SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_COLS_MAP)
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP)) THEN
        DO row_idx=1,SIZE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP,1)
          CALL SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_FINALISE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP( &
            & row_idx),ERR,ERROR,*999)
        ENDDO !row_idx
        DEALLOCATE(SOLVER_MAPPING%SOLVER_ROW_TO_EQUATIONS_ROWS_MAP)
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
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS=0
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

  !>Adds an interface condition to a solver mapping
  SUBROUTINE SOLVER_MAPPING_INTERFACE_CONDITION_ADD(SOLVER_MAPPING,INTERFACE_CONDITION,INTERFACE_CONDITION_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to add the interface condition to
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to add
    INTEGER(INTG), INTENT(OUT) :: INTERFACE_CONDITION_INDEX !<On exit, the index of the interface condition in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_matrix_idx,equations_set_idx,interface_condition_idx,interface_matrix_idx,LIST_ITEM(2)
    LOGICAL :: EQUATIONS_SET_FOUND,VARIABLE_FOUND
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_CONDITION_PTR_TYPE), ALLOCATABLE :: OLD_INTERFACE_CONDITIONS(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MAPPING_INTERFACE_CONDITION_ADD",ERR,ERROR,*999)

    INTERFACE_CONDITION_INDEX=0
    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solver mapping has been finished.",ERR,ERROR,*999)
      ELSE
        SOLVER_EQUATIONS=>SOLVER_MAPPING%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
            IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                IF(ASSOCIATED(INTERFACE_MAPPING)) THEN                
                  IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
                    !Check that the interface variables are already part of an added equations set.
                    SELECT CASE(INTERFACE_CONDITION%METHOD)
                    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
                      IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                        DO interface_matrix_idx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                          EQUATIONS_SET=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%EQUATIONS_SET
                          IF(ASSOCIATED(EQUATIONS_SET)) THEN
                            EQUATIONS_SET_FOUND=.FALSE.
                            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                              IF(ASSOCIATED(EQUATIONS_SET,SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR)) THEN
                                EQUATIONS_SET_FOUND=.TRUE.
                                EXIT
                              ENDIF
                            ENDDO !equations_set_idx
                            IF(EQUATIONS_SET_FOUND) THEN
                              !See if the variable is in the equations set.
                              DEPENDENT_VARIABLE=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE
                              IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                                VARIABLE_FOUND=.FALSE.
                                !Check dynamic variables
                                IF(SOLVER_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE(equations_set_idx)== &
                                  & DEPENDENT_VARIABLE%VARIABLE_TYPE) THEN
                                  VARIABLE_FOUND=.TRUE.
                                ELSE
                                  !Check linear matrices. Just check for solver matrix 1 and the moment
                                  DO equations_matrix_idx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                                    & MATRIX_VARIABLE_TYPES(0,equations_set_idx,1)
                                    IF(SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                                      & MATRIX_VARIABLE_TYPES(equations_matrix_idx,equations_set_idx,1)== &
                                      & DEPENDENT_VARIABLE%VARIABLE_TYPE) THEN
                                      VARIABLE_FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !equations matrix_idx
                                  IF(.NOT.VARIABLE_FOUND) THEN
                                    !Check residual variable type
                                    IF(SOLVER_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE(equations_set_idx)== &
                                      & DEPENDENT_VARIABLE%VARIABLE_TYPE) THEN
                                      VARIABLE_FOUND=.TRUE.
                                    ENDIF
                                  ENDIF
                                ENDIF
                                IF(VARIABLE_FOUND) THEN
                                  !Add in interface condition to equations set (just for solver matrix 1 at the moment)
                                  LIST_ITEM(1)=SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS+1
                                  LIST_ITEM(2)=interface_matrix_idx
                                  CALL LIST_ITEM_ADD(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_INDICES(equations_set_idx)% &
                                    & PTR,LIST_ITEM,ERR,ERROR,*999)                                  
                                ELSE
                                  LOCAL_ERROR="The dependent variable associated with interface matrix number "// &
                                    & TRIM(NUMBER_TO_VSTRING(interface_matrix_idx,"*",ERR,ERROR))// &
                                    & " is not mapped to the solver equations."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                LOCAL_ERROR="The dependent variable associated with interface matrix number "// &
                                  & TRIM(NUMBER_TO_VSTRING(interface_matrix_idx,"*",ERR,ERROR))//" is not associated."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              LOCAL_ERROR="The equations set for the dependent variable associated with interface "// &
                                & "matrix number "//TRIM(NUMBER_TO_VSTRING(interface_matrix_idx,"*",ERR,ERROR))// &
                                & " has not been added to the solver equations."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="Equations set is not associated for interface matrix number "// &
                              & TRIM(NUMBER_TO_VSTRING(interface_matrix_idx,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                        ENDDO !interface_matrix_idx
                      ELSE
                        CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interface condition method of "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT                    
                    IF(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS>0) THEN
                      ALLOCATE(OLD_INTERFACE_CONDITIONS(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interface conditions.",ERR,ERROR,*999)
                      DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
                        OLD_INTERFACE_CONDITIONS(interface_condition_idx)%PTR=>SOLVER_MAPPING% &
                          & INTERFACE_CONDITIONS(interface_condition_idx)%PTR
                      ENDDO !interface_condition_idx
                    ELSE IF(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS==0) THEN
                      !Do nothing
                    ELSE
                      CALL FLAG_ERROR("The number of interface conditions is < 0.",ERR,ERROR,*999)
                    ENDIF
                    ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITIONS(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS+1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface conditions.",ERR,ERROR,*999)
                    DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
                      SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR=> &
                        & OLD_INTERFACE_CONDITIONS(interface_condition_idx)%PTR
                    ENDDO !interface_condition_idx
                    SOLVER_MAPPING%INTERFACE_CONDITIONS(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS+1)%PTR=>INTERFACE_CONDITION
                    SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS=SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS+1
                    INTERFACE_CONDITION_INDEX=SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS

!!TODO: SORT OUT LAGRANGE FIELD VARIABLE
                    CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_INTERF_VAR_LIST_ADD(SOLVER_MAPPING,1,INTERFACE_CONDITION_INDEX, &
                      & 1,ERR,ERROR,*999)
                    
                    IF(ALLOCATED(OLD_INTERFACE_CONDITIONS)) DEALLOCATE(OLD_INTERFACE_CONDITIONS)
                  ELSE
                    CALL FLAG_ERROR("Solvers mapping create values cache is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Interface equations mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface condition has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver mapping solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_INTERFACE_CONDITION_ADD")
    RETURN
999 IF(ALLOCATED(OLD_INTERFACE_CONDITIONS)) DEALLOCATE(OLD_INTERFACE_CONDITIONS)
    CALL ERRORS("SOLVER_MAPPING_INTERFACE_CONDITION_ADD",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_CONDITION_ADD")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_CONDITION_ADD

  !
  !================================================================================================================================
  !

  !>Finalises an interface condition to solver map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERF_CONDITION_TO_SOLVER_MAP_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE) :: INTERFACE_CONDITION_TO_SOLVER_MAP !<The interface condition to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,equations_set_idx,interface_matrix_idx,solver_matrix_idx
    
    CALL ENTERS("SOLVER_MAPPING_INTERF_CONDITION_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_CONDITION_INDEX=0
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%SOLVER_MAPPING)
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_EQUATIONS)
    INTERFACE_CONDITION_TO_SOLVER_MAP%NUMBER_OF_EQUATIONS_SETS=0
    IF(ALLOCATED(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS)) THEN
      DO equations_set_idx=1,SIZE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS,1)
        CALL SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP% &
          & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_set_idx),ERR,ERROR,*999)
      ENDDO !equations_set_idx
      DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS)
    ENDIF
    IF(ALLOCATED(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM)) THEN
      DO solver_matrix_idx=1,SIZE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM,1)
        CALL SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP% &
          & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx),ERR,ERROR,*999)
      ENDDO !solver_matrix_idx
      DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_SM)
    ENDIF
    IF(ALLOCATED(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM)) THEN
      DO interface_matrix_idx=1,SIZE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM,1)
        CALL SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP% &
          & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx),ERR,ERROR,*999)
      ENDDO !interface_matrix_idx
      DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_IM)
    ENDIF
    IF(ALLOCATED(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS)) THEN
      DO column_idx=1,SIZE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS,1)
        CALL SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP% &
          & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(column_idx),ERR,ERROR,*999)
      ENDDO !column_idx
      DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_INTERF_CONDITION_TO_SOLVER_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_CONDITION_TO_SOLVER_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_CONDITION_TO_SOLVER_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_CONDITION_TO_SOLVER_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition to solver map.
  SUBROUTINE SOLVER_MAPPING_INTERF_CONDITON_TO_SOLVER_MAP_INITIALISE(INTERFACE_CONDITION_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE) :: INTERFACE_CONDITION_TO_SOLVER_MAP !<The interface condition to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_INTERF_CONDITON_TO_SOLVER_MAP_INITIALISE",ERR,ERROR,*999)

    INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_CONDITION_INDEX=0
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%SOLVER_MAPPING)
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_EQUATIONS)
    INTERFACE_CONDITION_TO_SOLVER_MAP%NUMBER_OF_EQUATIONS_SETS=0
    
    CALL EXITS("SOLVER_MAPPING_INTERF_CONDITON_TO_SOLVER_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_CONDITON_TO_SOLVER_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_CONDITON_TO_SOLVER_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_CONDITON_TO_SOLVER_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a interface condition to solver map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE), POINTER :: INTERFACE_CONDITION_TO_SOLVER_MAP !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,interface_column_idx
    
    CALL ENTERS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION_TO_SOLVER_MAP)) THEN
      IF(ALLOCATED(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS)) THEN
        DO equations_set_idx=1,SIZE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS,1)
          CALL SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP% &
            & INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_set_idx),ERR,ERROR,*999)
        ENDDO !equations_set_idx
        DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS)
      ENDIF
      IF(ALLOCATED(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS)) THEN
        DO interface_column_idx=1,SIZE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS,1)
          CALL SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE(INTERFACE_CONDITION_TO_SOLVER_MAP% &
            & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(interface_column_idx),ERR,ERROR,*999)
        ENDDO !interface_column_idx
        DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS)
      ENDIF
      DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP)
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition to solver map.
  SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_INITIALISE(INTERFACE_CONDITION_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE) :: INTERFACE_CONDITION_TO_SOLVER_MAP !<The interface condition to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_INITIALISE",ERR,ERROR,*999)

    INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_CONDITION_INDEX=0
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%SOLVER_MAPPING)
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_EQUATIONS)
    INTERFACE_CONDITION_TO_SOLVER_MAP%NUMBER_OF_EQUATIONS_SETS=0
    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAP_INITIALISE

 !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix maps and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE(INTERFACE_TO_SOLVER_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MAPS_TYPE), POINTER :: INTERFACE_TO_SOLVER_MAPS !<The interface to solver maps to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: row_idx
    
    CALL ENTERS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_TO_SOLVER_MAPS)) THEN
      IF(ALLOCATED(INTERFACE_TO_SOLVER_MAPS%INTERFACE_ROW_TO_SOLVER_COLS_MAP)) THEN
        DO row_idx=1,SIZE(INTERFACE_TO_SOLVER_MAPS%INTERFACE_ROW_TO_SOLVER_COLS_MAP,1)
          CALL SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(INTERFACE_TO_SOLVER_MAPS% &
            & INTERFACE_ROW_TO_SOLVER_COLS_MAP(row_idx),ERR,ERROR,*999)
        ENDDO !row_idx
        DEALLOCATE(INTERFACE_TO_SOLVER_MAPS%INTERFACE_ROW_TO_SOLVER_COLS_MAP)
      ENDIF
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver maps
  SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_INITIALISE(INTERFACE_TO_SOLVER_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MAPS_TYPE), POINTER :: INTERFACE_TO_SOLVER_MAPS !<The interface to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_TO_SOLVER_MAPS)) THEN
      INTERFACE_TO_SOLVER_MAPS%INTERFACE_MATRIX_NUMBER=0
      INTERFACE_TO_SOLVER_MAPS%SOLVER_MATRIX_NUMBER=0
      NULLIFY(INTERFACE_TO_SOLVER_MAPS%INTERFACE_MATRIX)
      NULLIFY(INTERFACE_TO_SOLVER_MAPS%SOLVER_MATRIX)
    ELSE
      CALL FLAG_ERROR("Interface to solver maps is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix equations map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE(INTERFACE_TO_SOLVER_EQUATIONS_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS_TYPE) :: INTERFACE_TO_SOLVER_EQUATIONS_MAPS !<The interface to solver equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE",ERR,ERROR,*999)

    INTERFACE_TO_SOLVER_EQUATIONS_MAPS%EQUATIONS_SET_INDEX=0
    NULLIFY(INTERFACE_TO_SOLVER_EQUATIONS_MAPS%EQUATIONS_SET)
    INTERFACE_TO_SOLVER_EQUATIONS_MAPS%INTERFACE_MATRIX_INDEX=0
         
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver matrix equations map.
  SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_INITIALISE(INTERFACE_TO_SOLVER_EQUATIONS_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS_TYPE) :: INTERFACE_TO_SOLVER_EQUATIONS_MAPS !<The interface to solver equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_INITIALISE",ERR,ERROR,*999)

    INTERFACE_TO_SOLVER_EQUATIONS_MAPS%EQUATIONS_SET_INDEX=0
    NULLIFY(INTERFACE_TO_SOLVER_EQUATIONS_MAPS%EQUATIONS_SET)
    INTERFACE_TO_SOLVER_EQUATIONS_MAPS%INTERFACE_MATRIX_INDEX=0
         
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_TO_SOLVER_EQUATIONS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a interface column to solver row map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE(INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP_TYPE) :: INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP !<The interface column to solver row map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR,*999)

    INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP%NUMBER_OF_SOLVER_ROWS=0
    INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP%SOLVER_ROW=0
    INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP%COUPLING_COEFFICIENT=0.0_DP
         
    CALL EXITS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises am interface column to solver row map.
  SUBROUTINE SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_INITIALISE(INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP_TYPE) :: INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP !<The interface column to solver row map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR,*999)

    INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP%NUMBER_OF_SOLVER_ROWS=0
    INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP%SOLVER_ROW=0
    INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP%COUPLING_COEFFICIENT=0.0_DP
        
    CALL EXITS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_COL_TO_SOL_ROWS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a interface row to solver row map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_FINALISE(INTERFACE_ROW_TO_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_ROW_TO_SOLVER_ROWS_MAP_TYPE) :: INTERFACE_ROW_TO_SOLVER_ROWS_MAP !<The interface row to solver row map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR,*999)

    INTERFACE_ROW_TO_SOLVER_ROWS_MAP%NUMBER_OF_SOLVER_ROWS=0
    INTERFACE_ROW_TO_SOLVER_ROWS_MAP%SOLVER_ROW=0
    INTERFACE_ROW_TO_SOLVER_ROWS_MAP%COUPLING_COEFFICIENT=0.0_DP
        
    CALL EXITS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises am interface row to solver row map.
  SUBROUTINE SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_INITIALISE(INTERFACE_ROW_TO_SOLVER_ROWS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_ROW_TO_SOLVER_ROWS_MAP_TYPE) :: INTERFACE_ROW_TO_SOLVER_ROWS_MAP !<The interface row to solver row map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR,*999)

    INTERFACE_ROW_TO_SOLVER_ROWS_MAP%NUMBER_OF_SOLVER_ROWS=0
    INTERFACE_ROW_TO_SOLVER_ROWS_MAP%SOLVER_ROW=0
    INTERFACE_ROW_TO_SOLVER_ROWS_MAP%COUPLING_COEFFICIENT=0.0_DP
        
    CALL EXITS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix map im and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM_TYPE) :: INTERFACE_TO_SOLVER_MATRIX_MAPS_IM !<The interface to solver matrix maps Im to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: row_idx,solver_matrix_idx
    
    CALL ENTERS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_FINALISE",ERR,ERROR,*999)

    INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_MATRIX_NUMBER=0
    INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%NUMBER_OF_SOLVER_MATRICES=0
    IF(ALLOCATED(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_TO_SOLVER_MATRIX_MAPS)) THEN
      DO solver_matrix_idx=1,SIZE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_TO_SOLVER_MATRIX_MAPS,1)
        CALL SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM% &
          & INTERFACE_TO_SOLVER_MATRIX_MAPS(solver_matrix_idx)%PTR,ERR,ERROR,*999)        
      ENDDO !solver_matrix_idx
      DEALLOCATE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_TO_SOLVER_MATRIX_MAPS)
    ENDIF
    IF(ALLOCATED(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_ROW_TO_SOLVER_ROWS_MAP)) THEN
      DO row_idx=1,SIZE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_ROW_TO_SOLVER_ROWS_MAP,1)
        CALL SOLVER_MAPPING_INTERF_ROW_TO_SOL_ROWS_MAP_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM% &
          INTERFACE_ROW_TO_SOLVER_ROWS_MAP(row_idx),ERR,ERROR,*999)
      ENDDO !row_idx
      DEALLOCATE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_ROW_TO_SOLVER_ROWS_MAP)
    ENDIF
    
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver matrix maps im.
  SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_INITIALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM_TYPE) :: INTERFACE_TO_SOLVER_MATRIX_MAPS_IM !<The interface to solver matrix maps Im to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_INITIALISE",ERR,ERROR,*999)

    INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%INTERFACE_MATRIX_NUMBER=0
    INTERFACE_TO_SOLVER_MATRIX_MAPS_IM%NUMBER_OF_SOLVER_MATRICES=0
    
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_IM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix map sm and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM_TYPE) :: INTERFACE_TO_SOLVER_MATRIX_MAPS_SM !<The interface to solver matrix maps sm to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,interface_matrix_idx
    
    CALL ENTERS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_FINALISE",ERR,ERROR,*999)

    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%SOLVER_MATRIX_NUMBER=0
    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%LAGRANGE_VARIABLE_TYPE=0
    NULLIFY(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%LAGRANGE_VARIABLE)
    CALL SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM% &
      & LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP,ERR,ERROR,*999)
    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_DEPENDENT_VARIABLES=0
    IF(ALLOCATED(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%DEPENDENT_VARIABLE_TYPES)) &
      & DEALLOCATE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%DEPENDENT_VARIABLE_TYPES)
    IF(ALLOCATED(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%DEPENDENT_VARIABLES)) &
      & DEALLOCATE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%DEPENDENT_VARIABLES)
    IF(ALLOCATED(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS)) THEN
      DO interface_matrix_idx=1,SIZE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS,1)
        CALL SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM% &
          & DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(interface_matrix_idx),ERR,ERROR,*999)
      ENDDO !interface_matrix_idx
      DEALLOCATE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS)
    ENDIF
    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_INTERFACE_MATRICES=0
    IF(ALLOCATED(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS)) THEN
      DO interface_matrix_idx=1,SIZE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS,1)
        CALL SOLVER_MAPPING_INTERFACE_TO_SOLVER_MAPS_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM% &
          INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx)%PTR,ERR,ERROR,*999)
      ENDDO !interface_matrix_idx
      DEALLOCATE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS)
    ENDIF
    IF(ALLOCATED(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%INTERFACE_COL_TO_SOLVER_COLS_MAP)) THEN
      DO column_idx=1,SIZE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%INTERFACE_COL_TO_SOLVER_COLS_MAP,1)
        CALL SOLVER_MAPPING_EQUATS_COL_TO_SOL_COLS_MAP_FINALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM% &
          & INTERFACE_COL_TO_SOLVER_COLS_MAP(column_idx),ERR,ERROR,*999)
      ENDDO !column_idx
      DEALLOCATE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%INTERFACE_COL_TO_SOLVER_COLS_MAP)
    ENDIF
   
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver matrix maps sm.
  SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM_TYPE) :: INTERFACE_TO_SOLVER_MATRIX_MAPS_SM !<The interface to solver matrix maps sm to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE",ERR,ERROR,*999)

    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%SOLVER_MATRIX_NUMBER=0
    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%LAGRANGE_VARIABLE_TYPE=0
    NULLIFY(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%LAGRANGE_VARIABLE)
    CALL SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_INITIALISE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM% &
      & LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP,ERR,ERROR,*999)
    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_DEPENDENT_VARIABLES=0
    INTERFACE_TO_SOLVER_MATRIX_MAPS_SM%NUMBER_OF_INTERFACE_MATRICES=0
        
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERF_TO_SOL_MAT_MAPS_SM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a Jacobian column to solver columns map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE(JACOBIAN_COL_TO_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE) :: JACOBIAN_COL_TO_SOLVER_COLS_MAP !<The Jacobian col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(JACOBIAN_COL_TO_SOLVER_COLS_MAP%SOLVER_COLS)) &
      & DEALLOCATE(JACOBIAN_COL_TO_SOLVER_COLS_MAP%SOLVER_COLS)
    IF(ALLOCATED(JACOBIAN_COL_TO_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(JACOBIAN_COL_TO_SOLVER_COLS_MAP%COUPLING_COEFFICIENTS)
        
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
  SUBROUTINE SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE(JACOBIAN_COL_TO_SOLVER_COLS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE) :: JACOBIAN_COL_TO_SOLVER_COLS_MAP !<The Jacobian column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_INITIALISE",ERR,ERROR,*999)

    JACOBIAN_COL_TO_SOLVER_COLS_MAP%NUMBER_OF_SOLVER_COLS=0
    
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
      IF(ALLOCATED(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP)) THEN
        DO column_idx=1,SIZE(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP,1)
          CALL SOLVER_MAPPING_JAC_COL_TO_SOL_COLS_MAP_FINALISE(JACOBIAN_TO_SOLVER_MAP% &
            & JACOBIAN_COL_TO_SOLVER_COLS_MAP(column_idx),ERR,ERROR,*999)
        ENDDO !column_idx
        DEALLOCATE(JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_TO_SOLVER_COLS_MAP)
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
      CALL FLAG_ERROR("Jacobian to solver matrix map is not associated.",ERR,ERROR,*999)
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
      SOLVER_COL_TO_EQUATIONS_SET_MAP%HAVE_DYNAMIC=.FALSE.
    ENDIF
    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS)) THEN
      DO col_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS,1)
        CALL SOLVER_MAPPING_SOLVER_COL_TO_S_EQUATIONS_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_SET_MAP% &
          & SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(col_idx),ERR,ERROR,*999)
      ENDDO !col_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_SET_MAP%SOLVER_COL_TO_STATIC_EQUATIONS_MAPS)
      SOLVER_COL_TO_EQUATIONS_SET_MAP%HAVE_STATIC=.FALSE.
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
    SOLVER_COL_TO_EQUATIONS_SET_MAP%HAVE_DYNAMIC=.FALSE.
    SOLVER_COL_TO_EQUATIONS_SET_MAP%HAVE_STATIC=.FALSE.
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_FINALISE(SOLVER_COL_TO_EQUATIONS_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_MAPS_TYPE) :: SOLVER_COL_TO_EQUATIONS_MAPS !<The solver column to equations sets map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: col_idx,equations_set_idx
    
    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_COL_TO_EQUATIONS_SET_MAPS)) THEN
      DO equations_set_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_COL_TO_EQUATIONS_SET_MAPS,1)
        CALL SOLVER_MAPPING_SOL_COL_TO_EQUATS_SET_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_MAPS% &
          SOLVER_COL_TO_EQUATIONS_SET_MAPS(equations_set_idx),ERR,ERROR,*999)
      ENDDO !equations_set_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_COL_TO_EQUATIONS_SET_MAPS)
    ENDIF
    IF(ALLOCATED(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_DOF_TO_VARIABLE_MAPS)) THEN
      DO col_idx=1,SIZE(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_DOF_TO_VARIABLE_MAPS,1)
        CALL SOLVER_MAPPING_SOLVER_DOF_TO_VARIABLE_MAP_FINALISE(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_DOF_TO_VARIABLE_MAPS( &
          & col_idx),ERR,ERROR,*999)
      ENDDO !col_idx
      DEALLOCATE(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_DOF_TO_VARIABLE_MAPS)
    ENDIF
    CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(SOLVER_COL_TO_EQUATIONS_MAPS%COLUMN_DOFS_MAPPING,ERR,ERROR,*999)
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_INITIALISE(SOLVER_COL_TO_EQUATIONS_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_EQUATIONS_MAPS_TYPE) :: SOLVER_COL_TO_EQUATIONS_MAPS !<The solver column to equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_INITIALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_MATRIX_NUMBER=0
    NULLIFY(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_MATRIX)
    NULLIFY(SOLVER_COL_TO_EQUATIONS_MAPS%SOLVER_MAPPING)
    SOLVER_COL_TO_EQUATIONS_MAPS%NUMBER_OF_COLUMNS=0
    SOLVER_COL_TO_EQUATIONS_MAPS%NUMBER_OF_DOFS=0
    SOLVER_COL_TO_EQUATIONS_MAPS%TOTAL_NUMBER_OF_DOFS=0
    SOLVER_COL_TO_EQUATIONS_MAPS%NUMBER_OF_GLOBAL_DOFS=0
    NULLIFY(SOLVER_COL_TO_EQUATIONS_MAPS%COLUMN_DOFS_MAPPING)
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_EQUATIONS_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to interface map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_FINALISE(SOLVER_COL_TO_INTERFACE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_INTERFACE_MAP_TYPE) :: SOLVER_COL_TO_INTERFACE_MAP !<The solver column to interface map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx
    
    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_FINALISE",ERR,ERROR,*999)

    NULLIFY(SOLVER_COL_TO_INTERFACE_MAP%INTERFACE_EQUATIONS)
    IF(ALLOCATED(SOLVER_COL_TO_INTERFACE_MAP%SOLVER_COL_TO_INTERFACE_EQUATIONS_MAPS)) THEN
      DO column_idx=1,SIZE(SOLVER_COL_TO_INTERFACE_MAP%SOLVER_COL_TO_INTERFACE_EQUATIONS_MAPS,1)
        CALL SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_FINALISE(SOLVER_COL_TO_INTERFACE_MAP% &
          & SOLVER_COL_TO_INTERFACE_EQUATIONS_MAPS(column_idx),ERR,ERROR,*999)
      ENDDO !column_idx
      DEALLOCATE(SOLVER_COL_TO_INTERFACE_MAP%SOLVER_COL_TO_INTERFACE_EQUATIONS_MAPS)
    ENDIF
     
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to interface mapping.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_INITIALISE(SOLVER_COL_TO_INTERFACE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_INTERFACE_MAP_TYPE) :: SOLVER_COL_TO_INTERFACE_MAP !<The solver column to interface map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_INITIALISE",ERR,ERROR,*999)

    NULLIFY(SOLVER_COL_TO_INTERFACE_MAP%INTERFACE_EQUATIONS)
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_INITIALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to interface equations map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_FINALISE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP !<The solver column to interface equatiosn map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_FINALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%NUMBER_OF_INTERFACE_MATRICES=0
    IF(ALLOCATED(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%INTERFACE_MATRIX_NUMBERS))  &
      & DEALLOCATE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%INTERFACE_MATRIX_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%INTERFACE_COL_NUMBERS))  &
      & DEALLOCATE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%INTERFACE_COL_NUMBERS)
    IF(ALLOCATED(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%COUPLING_COEFFICIENTS))  &
      & DEALLOCATE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%COUPLING_COEFFICIENTS)
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to interface equations mapping.
  SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_INITIALISE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP_TYPE) :: SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP !<The solver column to interface equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP%NUMBER_OF_INTERFACE_MATRICES=0
    
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_INITIALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_MAPPING_SOL_COL_TO_INTERF_EQUATS_MAP_INITIALISE

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

    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_TYPES)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_TYPES)
    IF(ALLOCATED(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_INDICES)) DEALLOCATE(SOLVER_DOF_TO_VARIABLE_MAP%EQUATIONS_INDICES)
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

    SOLVER_DOF_TO_VARIABLE_MAP%NUMBER_OF_EQUATIONS=0
    
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

  !>Finalises a solver row to equations map and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_FINALISE(SOLVER_ROW_TO_EQUATIONS_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_ROW_TO_EQUATIONS_MAPS_TYPE) :: SOLVER_ROW_TO_EQUATIONS_MAPS !<The solver row to equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_FINALISE",ERR,ERROR,*999)

    SOLVER_ROW_TO_EQUATIONS_MAPS%NUMBER_OF_EQUATIONS_SETS=0
    SOLVER_ROW_TO_EQUATIONS_MAPS%INTERFACE_CONDITION_INDEX=0
    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_MAPS%EQUATIONS_INDEX)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_MAPS%EQUATIONS_INDEX)
    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_MAPS%ROWCOL_NUMBER)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_MAPS%ROWCOL_NUMBER)
    IF(ALLOCATED(SOLVER_ROW_TO_EQUATIONS_MAPS%COUPLING_COEFFICIENTS)) &
      & DEALLOCATE(SOLVER_ROW_TO_EQUATIONS_MAPS%COUPLING_COEFFICIENTS)
    
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_FINALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solver row to equations map.
  SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE(SOLVER_ROW_TO_EQUATIONS_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_ROW_TO_EQUATIONS_MAPS_TYPE) :: SOLVER_ROW_TO_EQUATIONS_MAPS !<The solver row to equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE",ERR,ERROR,*999)

    SOLVER_ROW_TO_EQUATIONS_MAPS%NUMBER_OF_EQUATIONS_SETS=0
    SOLVER_ROW_TO_EQUATIONS_MAPS%INTERFACE_CONDITION_INDEX=0
        
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_MAPPING_SOL_ROW_TO_EQUATIONS_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variable and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_VARIABLE_FINALISE(SOLVER_MAPPING_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE) :: SOLVER_MAPPING_VARIABLE !<The solver mapping variable to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MAPPING_VARIABLE_FINALISE",ERR,ERROR,*999)

    NULLIFY(SOLVER_MAPPING_VARIABLE%VARIABLE)
    SOLVER_MAPPING_VARIABLE%VARIABLE_TYPE=0
    SOLVER_MAPPING_VARIABLE%NUMBER_OF_EQUATIONS=0
    IF(ALLOCATED(SOLVER_MAPPING_VARIABLE%EQUATION_TYPES)) DEALLOCATE(SOLVER_MAPPING_VARIABLE%EQUATION_TYPES)
    IF(ALLOCATED(SOLVER_MAPPING_VARIABLE%EQUATION_INDICES)) DEALLOCATE(SOLVER_MAPPING_VARIABLE%EQUATION_INDICES)
       
    CALL EXITS("SOLVER_MAPPING_VARIABLE_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_VARIABLE_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_VARIABLE_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_VARIABLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_VARIABLE_INITIALISE(SOLVER_MAPPING_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE) :: SOLVER_MAPPING_VARIABLE !<The solver mapping variable to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("SOLVER_MAPPING_VARIABLE_INITIALISE",ERR,ERROR,*999)

    NULLIFY(SOLVER_MAPPING_VARIABLE%VARIABLE)
    SOLVER_MAPPING_VARIABLE%VARIABLE_TYPE=0
    SOLVER_MAPPING_VARIABLE%NUMBER_OF_EQUATIONS=0
    
    CALL EXITS("SOLVER_MAPPING_VARIABLE_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_VARIABLE_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_VARIABLE_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_VARIABLE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variables and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_VARIABLES_FINALISE(SOLVER_MAPPING_VARIABLES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_VARIABLES_TYPE) :: SOLVER_MAPPING_VARIABLES !<The solver mapping variables to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx

    CALL ENTERS("SOLVER_MAPPING_VARIABLES_FINALISE",ERR,ERROR,*999)

    SOLVER_MAPPING_VARIABLES%NUMBER_OF_VARIABLES=0
    DO variable_idx=1,SIZE(SOLVER_MAPPING_VARIABLES%VARIABLES,1)
      CALL SOLVER_MAPPING_VARIABLE_FINALISE(SOLVER_MAPPING_VARIABLES%VARIABLES(variable_idx),ERR,ERROR,*999)
    ENDDO !variable_idx
     
    CALL EXITS("SOLVER_MAPPING_VARIABLES_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_VARIABLES_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_VARIABLES_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_VARIABLES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping variables and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_VARIABLES_INITIALISE(SOLVER_MAPPING_VARIABLES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_VARIABLES_TYPE) :: SOLVER_MAPPING_VARIABLES !<The solver mapping variables to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("SOLVER_MAPPING_VARIABLES_INITIALISE",ERR,ERROR,*999)

    SOLVER_MAPPING_VARIABLES%NUMBER_OF_VARIABLES=0
    
    CALL EXITS("SOLVER_MAPPING_VARIABLES_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_VARIABLES_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVER_MAPPING_VARIABLES_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_VARIABLES_INITIALISE

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

  !>Initialises a variable to solver column map.
  SUBROUTINE SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_INITIALISE(VARIABLE_TO_SOLVER_COL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VARIABLE_TO_SOLVER_COL_MAP_TYPE) :: VARIABLE_TO_SOLVER_COL_MAP !<The variable to solver column map to initalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_INITIALISE",ERR,ERROR,*999)

       
    CALL EXITS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_VARIABLE_TO_SOLVER_COL_MAP_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE SOLVER_MAPPING_ROUTINES

