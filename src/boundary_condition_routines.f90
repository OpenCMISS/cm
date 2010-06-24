!> \file
!> $Id$
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
  USE BASIS_ROUTINES
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE COORDINATE_ROUTINES
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
  USE LISTS
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions BOUNDARY_CONDITIONS_ROUTINES::BoundaryConditions
  !> \brief Boundary conditions type parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NOT_FIXED=0 !<The dof is not fixed. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED=1 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INLET=2 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_OUTLET=3 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_WALL=4 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL=5 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE_WALL=6 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MIXED=7 !<The dof is set as a mixed boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN=8 !<The dof is set to a Neumann boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT=9 !<The dof is set to a Neumann point boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DIRICHLET=10 !<The dof is set to a Dirichlet boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CAUCHY=11 !<The dof is set to a Cauchy boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_ROBIN=12 !<The dof is set to a Robin boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INCREMENTED_TYPE=13 !<The dof is a fixed boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Adds to the value of the specified local DOF and sets this as a boundary condition on the specified local DOF.
  INTERFACE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF
    MODULE PROCEDURE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1
    MODULE PROCEDURE BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS
  END INTERFACE !BOUNDARY_CONDITIONS_ADD_LOCAL_DOF

  !>Sets a boundary condition on the specified local DOF. 
  INTERFACE BOUNDARY_CONDITIONS_SET_LOCAL_DOF
    MODULE PROCEDURE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1
    MODULE PROCEDURE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS
  END INTERFACE !BOUNDARY_CONDITIONS_SET_LOCAL_DOF

  !>Add a specified dof(s) to a specified boundary condition
  INTERFACE BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF
    MODULE PROCEDURE BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF1
    MODULE PROCEDURE BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS
  END INTERFACE !BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF

  PUBLIC BOUNDARY_CONDITION_NOT_FIXED,BOUNDARY_CONDITION_FIXED,BOUNDARY_CONDITION_MIXED,BOUNDARY_CONDITION_FIXED_INLET, & 
    & BOUNDARY_CONDITION_FIXED_OUTLET,BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL,BOUNDARY_CONDITION_FREE_WALL, &
    & BOUNDARY_CONDITION_NEUMANN,BOUNDARY_CONDITION_DIRICHLET,BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_CAUCHY, &
    & BOUNDARY_CONDITION_ROBIN, BOUNDARY_CONDITION_FIXED_INCREMENTED_TYPE

  PUBLIC BOUNDARY_CONDITIONS_CREATE_FINISH,BOUNDARY_CONDITIONS_CREATE_START,BOUNDARY_CONDITIONS_DESTROY
  
  PUBLIC BOUNDARY_CONDITIONS_ADD_CONSTANT,BOUNDARY_CONDITIONS_ADD_LOCAL_DOF,BOUNDARY_CONDITIONS_ADD_ELEMENT, &
    & BOUNDARY_CONDITIONS_ADD_NODE

  PUBLIC BOUNDARY_CONDITIONS_SET_CONSTANT,BOUNDARY_CONDITIONS_SET_LOCAL_DOF,BOUNDARY_CONDITIONS_SET_ELEMENT, &
    & BOUNDARY_CONDITIONS_SET_NODE, BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF, BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE, &
    & BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE, BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE, &
    & BOUNDARY_CONDITIONS_BOUNDARY_SET_NUMBER_OF_BOUNDARIES

  PUBLIC EQUATIONS_SET_BOUNDARY_CONDITIONS_GET
  
CONTAINS  

  !
  !================================================================================================================================
  !

  !>Finish the creation of boundary conditions.
  SUBROUTINE BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MPI_IERROR,SEND_COUNT,NUMBER_OF_DIRICHLET_CONDITIONS, STORAGE_TYPE, NUMBER_OF_NON_ZEROS, NUMBER_OF_ROWS, COUNT
    INTEGER(INTG) :: variable_type_idx,dof_idx, equ_matrix_idx, dirichlet_idx, sparse_idx, row_idx, DUMMY, LAST, DIRICHLET_DOF
    INTEGER(INTG), POINTER :: ROW_INDICES(:), COLUMN_INDICES(:)
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITION_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: VARIABLE_DOMAIN_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: BOUNDARY_CONDITIONS_DIRICHLET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATION_MATRIX
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES
    TYPE(LIST_TYPE), POINTER :: SPARSE_INDICES

    CALL ENTERS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have already been finished.",ERR,ERROR,*999)        
      ELSE
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP)) THEN
          IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>0) THEN
            !Transfer all the boundary conditions to all the computational nodes.
 !!TODO \todo Look at this. ?????
            DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
              BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type_idx)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
                FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(VARIABLE_DOMAIN_MAPPING)) THEN
                    SEND_COUNT=VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL
!!This operation is a little expensive as we are doing an unnecessary sum across all the ranks in order to combin
!!the data from each rank into all ranks. We will see how this goes for now.
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS, &
                      & SEND_COUNT,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="Field variable domain mapping is not associated for variable type "// &
                      & TRIM(NUMBER_TO_VSTRING(variable_type_idx,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  ! Find out how many dirichlet conditions in problem  \todo make below just fixed
                  NUMBER_OF_DIRICHLET_CONDITIONS=0
                  DO dof_idx=1,FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                    IF(BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED.OR. &
                      & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED_INLET.OR. &
                      & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED_OUTLET.OR. &
                      & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED_WALL.OR. &
                      & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_MOVED_WALL.OR. &
                      & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx) &
                      & ==BOUNDARY_CONDITION_FIXED_INCREMENTED_TYPE) THEN
                      NUMBER_OF_DIRICHLET_CONDITIONS=NUMBER_OF_DIRICHLET_CONDITIONS+1
                    ENDIF
                  ENDDO
                  BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS=NUMBER_OF_DIRICHLET_CONDITIONS
                  ! Check that there is at least one dirichlet condition
                  IF(NUMBER_OF_DIRICHLET_CONDITIONS>0) THEN
                    CALL BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                    BOUNDARY_CONDITIONS_DIRICHLET=>BOUNDARY_CONDITION_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_DIRICHLET)) THEN
                      ! Find dirichlet conditions \todo make below just fixed
                      dirichlet_idx=1
                      DO dof_idx=1,FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                        IF(BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED.OR. &
                          & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED_INLET.OR. &
                          & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED_OUTLET.OR. &
                          & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_FIXED_WALL.OR. &
                          & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx)==BOUNDARY_CONDITION_MOVED_WALL.OR. &
                          & BOUNDARY_CONDITION_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(dof_idx) &
                          & ==BOUNDARY_CONDITION_FIXED_INCREMENTED_TYPE) THEN
                          BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)=dof_idx
                          dirichlet_idx=dirichlet_idx+1
                        ENDIF
                      ENDDO

                      ! Store Dirichlet dof indices
                      EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        EQUATIONS=>EQUATIONS_SET%EQUATIONS
                        IF(ASSOCIATED(EQUATIONS)) THEN
                          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
                          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
                            IF(ASSOCIATED(LINEAR_MATRICES)) THEN
                              ! Iterate through equations matrices
                              DO equ_matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                                EQUATION_MATRIX=>LINEAR_MATRICES%MATRICES(equ_matrix_idx)%PTR
                                CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATION_MATRIX%MATRIX, STORAGE_TYPE, ERR, ERROR,*999)
                                IF(ASSOCIATED(EQUATION_MATRIX)) THEN
                                  SELECT CASE(STORAGE_TYPE)
                                  CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for block storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for diagonal storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for column major storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for row major storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                    ! Get Sparsity pattern, number of non zeros, number of rows
                                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATION_MATRIX%MATRIX,ROW_INDICES, &
                                      & COLUMN_INDICES,ERR,ERROR,*999)
                                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET(EQUATION_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
                                      & ERR,ERROR,*999)
                                    NUMBER_OF_ROWS=EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS
                                    ! Initialise sparsity indices arrays
                                    CALL BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE(BOUNDARY_CONDITIONS_DIRICHLET% &
                                      & LINEAR_SPARSITY_INDICES(equ_matrix_idx)%PTR,BOUNDARY_CONDITION_VARIABLE% &
                                      & NUMBER_OF_DIRICHLET_CONDITIONS,ERR,ERROR,*999)

                                    ! Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                                    NULLIFY(SPARSITY_INDICES)
                                    SPARSITY_INDICES=>BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES(equ_matrix_idx)%PTR
                                    IF(ASSOCIATED(SPARSITY_INDICES)) THEN
                                      ! Setup list for storing dirichlet non zero indices
                                      NULLIFY(SPARSE_INDICES)
                                      CALL LIST_CREATE_START(SPARSE_INDICES,ERR,ERROR,*999)
                                      CALL LIST_DATA_TYPE_SET(SPARSE_INDICES,LIST_INTG_TYPE,ERR,ERROR,*999)
                                      CALL LIST_INITIAL_SIZE_SET(SPARSE_INDICES, &
                                        & NUMBER_OF_DIRICHLET_CONDITIONS*(NUMBER_OF_NON_ZEROS/NUMBER_OF_ROWS),ERR,ERROR,*999)
                                      CALL LIST_CREATE_FINISH(SPARSE_INDICES,ERR,ERROR,*999)

                                      COUNT=0
                                      SPARSITY_INDICES%SPARSE_COLUMN_INDICES(1)=1
                                      LAST=1
                                      DO dirichlet_idx=1,BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
                                        DIRICHLET_DOF=BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)
                                        DO row_idx=1,NUMBER_OF_ROWS
                                          IF(DIRICHLET_DOF>=COLUMN_INDICES(ROW_INDICES(row_idx)).AND. &
                                            & DIRICHLET_DOF<=COLUMN_INDICES(ROW_INDICES(row_idx+1)-1)) THEN
                                            DO sparse_idx=ROW_INDICES(row_idx),ROW_INDICES(row_idx+1)-1
                                              IF(DIRICHLET_DOF==COLUMN_INDICES(sparse_idx)) THEN
                                                CALL LIST_ITEM_ADD(SPARSE_INDICES,row_idx,ERR,ERROR,*999)
                                                COUNT=COUNT+1
                                                LAST=row_idx+1
                                                EXIT
                                              ENDIF
                                            ENDDO
                                          ENDIF
                                        ENDDO
                                        SPARSITY_INDICES%SPARSE_COLUMN_INDICES(dirichlet_idx+1)=COUNT+1
                                      ENDDO
                                      CALL LIST_DETACH_AND_DESTROY(SPARSE_INDICES,DUMMY,SPARSITY_INDICES%SPARSE_ROW_INDICES, &
                                        & ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Sparsity indices arrays are not associated for this equations matrix"
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for compressed column storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for row column storage",ERR,ERROR,*999)
                                  CASE DEFAULT
                                    LOCAL_ERROR="The storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE,"*",ERR,ERROR)) &
                                      //" is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                ELSE
                                  CALL FLAG_ERROR("The equation matrix is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ENDDO
                            ENDIF

                            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
                            IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                              ! Iterate through equations matrices
                              DO equ_matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
                                EQUATION_MATRIX=>DYNAMIC_MATRICES%MATRICES(equ_matrix_idx)%PTR
                                CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATION_MATRIX%MATRIX, STORAGE_TYPE, ERR, ERROR,*999)
                                IF(ASSOCIATED(EQUATION_MATRIX)) THEN
                                  SELECT CASE(STORAGE_TYPE)
                                  CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for block storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for diagonal storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for column major storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for row major storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                    ! Get Sparsity pattern, number of non zeros, number of rows
                                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATION_MATRIX%MATRIX,ROW_INDICES, &
                                      & COLUMN_INDICES,ERR,ERROR,*999)
                                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET(EQUATION_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
                                      & ERR,ERROR,*999)
                                    NUMBER_OF_ROWS=EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS
                                    ! Intialise sparsity indices arrays
                                    CALL BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE(BOUNDARY_CONDITIONS_DIRICHLET% &
                                      & DYNAMIC_SPARSITY_INDICES(equ_matrix_idx)%PTR,BOUNDARY_CONDITION_VARIABLE% &
                                      & NUMBER_OF_DIRICHLET_CONDITIONS,ERR,ERROR,*999)

                                    ! Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                                    NULLIFY(SPARSITY_INDICES)
                                    SPARSITY_INDICES=>BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES(equ_matrix_idx)%PTR
                                    IF(ASSOCIATED(SPARSITY_INDICES)) THEN
                                      ! Setup list for storing dirichlet non zero indices
                                      NULLIFY(SPARSE_INDICES)
                                      CALL LIST_CREATE_START(SPARSE_INDICES,ERR,ERROR,*999)
                                      CALL LIST_DATA_TYPE_SET(SPARSE_INDICES,LIST_INTG_TYPE,ERR,ERROR,*999)
                                      CALL LIST_INITIAL_SIZE_SET(SPARSE_INDICES, &
                                        & NUMBER_OF_DIRICHLET_CONDITIONS*(NUMBER_OF_NON_ZEROS/NUMBER_OF_ROWS),ERR,ERROR,*999)
                                      CALL LIST_CREATE_FINISH(SPARSE_INDICES,ERR,ERROR,*999)

                                      COUNT=0
                                      SPARSITY_INDICES%SPARSE_COLUMN_INDICES(1)=1
                                      LAST=1
                                      DO dirichlet_idx=1,BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
                                        DIRICHLET_DOF=BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)
                                        DO row_idx=1,NUMBER_OF_ROWS
                                          IF(DIRICHLET_DOF>=COLUMN_INDICES(ROW_INDICES(row_idx)).AND. &
                                            & DIRICHLET_DOF<=COLUMN_INDICES(ROW_INDICES(row_idx+1)-1)) THEN
                                            DO sparse_idx=ROW_INDICES(row_idx),ROW_INDICES(row_idx+1)-1
                                              IF(DIRICHLET_DOF==COLUMN_INDICES(sparse_idx)) THEN
                                                CALL LIST_ITEM_ADD(SPARSE_INDICES,row_idx,ERR,ERROR,*999)
                                                COUNT=COUNT+1
                                                LAST=row_idx+1
                                                EXIT
                                              ENDIF
                                            ENDDO
                                          ENDIF
                                        ENDDO
                                        SPARSITY_INDICES%SPARSE_COLUMN_INDICES(dirichlet_idx+1)=COUNT+1
                                      ENDDO
                                      CALL LIST_DETACH_AND_DESTROY(SPARSE_INDICES,DUMMY,SPARSITY_INDICES%SPARSE_ROW_INDICES, &
                                        & ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Sparsity indices arrays are not associated for this equations matrix"
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for compressed column storage",ERR,ERROR,*999)
                                  CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                    CALL FLAG_ERROR("Not implemented for row column storage",ERR,ERROR,*999)
                                  CASE DEFAULT
                                    LOCAL_ERROR="The storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE,"*",ERR,ERROR)) &
                                      //" is invalid."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  END SELECT
                                ELSE
                                  CALL FLAG_ERROR("The equation matrix is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ENDDO
                            ENDIF
                          ELSE
                            LOCAL_ERROR="Equations Matrices is not associated for these Equations "
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          LOCAL_ERROR="Equations is not associated for this Equations Set "
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Equations Set is not associated for boundary conditions variable "// &
                          & TRIM(NUMBER_TO_VSTRING(variable_type_idx,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Dirichlet Boundary Conditions type is not associated for boundary condition variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(variable_type_idx,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  LOCAL_ERROR="Field variable is not associated for variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(variable_type_idx,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ENDDO ! variable_type_idx
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
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Variable type = ",variable_type_idx,ERR,ERROR,*999)
        IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
          FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
          VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
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
999 CALL ERRORS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR)
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
          SELECT CASE(EQUATIONS_SET%TYPE)
          CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
            !Initialise the boundary conditions for load increment loop
            CALL BOUNDARY_CONDITIONS_INITIALISE_LOAD_INCREMENT(EQUATIONS_SET,ERR,ERROR,*999)
          CASE DEFAULT
            !Initialise the boundary conditions
            CALL BOUNDARY_CONDITIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
          END SELECT
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
        ENDDO !variable_type_idx
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

  !>Initialises the boundary conditions for an equations set to be used with load increment loop.
  SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE_LOAD_INCREMENT(EQUATIONS_SET,ERR,ERROR,*)

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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING    
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ! The difference between normal BOUNDARY_CONDITIONS_INITIALISE and this one is that this one
    ! allocates a BOUNDARY_CONDITIONS_SET_TYPE for dependent field to store the conditions in the intrim.

    CALL ENTERS("BOUNDARY_CONDITIONS_INITIALISE_LOAD_INCREMENT",ERR,ERROR,*998)

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
                          VARIABLE => LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(linear_variable_type)%VARIABLE !
                          CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,VARIABLE,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_CREATE(VARIABLE%FIELD,linear_variable_type, &
                            & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)  !
                        ENDIF
                      ENDDO !linear_variable_idx
                    ELSE
                      CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,RHS_MAPPING% &
                        & RHS_VARIABLE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(RHS_MAPPING%RHS_VARIABLE%FIELD,RHS_MAPPING%RHS_VARIABLE_TYPE, &
                        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
                    ENDIF
                  CASE(EQUATIONS_NONLINEAR)
                    NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,NONLINEAR_MAPPING% &
                        & RESIDUAL_VARIABLE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%VARIABLE%FIELD, &
                        & NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999) !
                    ELSE
                      CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,RHS_MAPPING% &
                        & RHS_VARIABLE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(RHS_MAPPING%RHS_VARIABLE%FIELD,RHS_MAPPING%RHS_VARIABLE_TYPE, &
                        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999) !
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
                      CALL FIELD_PARAMETER_SET_CREATE(DYNAMIC_MAPPING%DYNAMIC_VARIABLE%FIELD, &
                        & DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999) !
                    ELSE
                      CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,RHS_MAPPING% &
                        & RHS_VARIABLE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(RHS_MAPPING%RHS_VARIABLE%FIELD,RHS_MAPPING%RHS_VARIABLE_TYPE, &
                        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999) !
                    ELSE
                      CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE(EQUATIONS_NONLINEAR)
                    DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                    IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,DYNAMIC_MAPPING% &
                        & DYNAMIC_VARIABLE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(DYNAMIC_MAPPING%DYNAMIC_VARIABLE%FIELD, &
                        & DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999) !
                    ELSE
                      CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,RHS_MAPPING% &
                        & RHS_VARIABLE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_CREATE(RHS_MAPPING%RHS_VARIABLE%FIELD,RHS_MAPPING%RHS_VARIABLE_TYPE, &
                        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999) !
                    ELSE
                      CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
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
       
    CALL EXITS("BOUNDARY_CONDITIONS_INITIALISE_LOAD_INCREMENT")
    RETURN
999 CALL BOUNDARY_CONDITIONS_FINALISE(EQUATIONS_SET%BOUNDARY_CONDITIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("BOUNDARY_CONDITIONS_INITIALISE_LOAD_INCREMENT",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_INITIALISE_LOAD_INCREMENT")
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE_LOAD_INCREMENT

  !
  !================================================================================================================================
  !
 
  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant. \see OPENCMISS::CMISSBoundaryConditionAddConstant
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_CONSTANT(BOUNDARY_CONDITIONS,VARIABLE_TYPE,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_CONSTANT",ERR,ERROR,*999)

    !Note: This routine is for constant interpolation
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
                CASE(BOUNDARY_CONDITION_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                CASE(BOUNDARY_CONDITION_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                  CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                    & VALUE,ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_MIXED)
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
       
    CALL EXITS("BOUNDARY_CONDITION_ADD_CONSTANT")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_ADD_CONSTANT",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_ADD_CONSTANT")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_CONSTANT
  
 !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified constant. \see OPENCMISS::CMISSBoundaryConditionsSetConstant
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

    !Note: This routine is for constant interpolation
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
                CASE(BOUNDARY_CONDITION_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                CASE(BOUNDARY_CONDITION_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                    & VALUE,ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_MIXED)
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
 
  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOF.
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1(BOUNDARY_CONDITIONS,VARIABLE_TYPE,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS(BOUNDARY_CONDITIONS,VARIABLE_TYPE,(/DOF_INDEX/),(/CONDITION/),(/VALUE/),ERR,ERROR,*999)
           
    CALL EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1
  
  !
  !================================================================================================================================
  !
 
  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOFs.
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS(BOUNDARY_CONDITIONS,VARIABLE_TYPE,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDICES(:) !<DOF_INDICES(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITIONS(:) !<CONDITIONS(:). The boundary condition type to set for the i'th dof \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(:). The value of the boundary condition for the i'th dof to add
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
 
    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS",ERR,ERROR,*999)

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
                          CASE(BOUNDARY_CONDITION_NOT_FIXED)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                          CASE(BOUNDARY_CONDITION_FIXED)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                            CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(BOUNDARY_CONDITION_MIXED)
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
       
    CALL EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS
  
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
           
    CALL EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1")
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
                          CASE(BOUNDARY_CONDITION_NOT_FIXED)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                          CASE(BOUNDARY_CONDITION_FIXED)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(BOUNDARY_CONDITION_FIXED_INLET)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED_INLET
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED_OUTLET
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(BOUNDARY_CONDITION_FIXED_WALL)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED_WALL
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(BOUNDARY_CONDITION_MOVED_WALL)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_MOVED_WALL
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(BOUNDARY_CONDITION_FREE_WALL)
                            BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FREE_WALL
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & local_ny,VALUES(i),ERR,ERROR,*999)
                          CASE(BOUNDARY_CONDITION_MIXED)
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
       
    CALL EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS
  
  !
  !================================================================================================================================
  !
 
  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsAddElement
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_ELEMENT(BOUNDARY_CONDITIONS,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
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
 
    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_ELEMENT",ERR,ERROR,*999)

    !Note: this routine is for element based interpolation
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
                CASE(BOUNDARY_CONDITION_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                CASE(BOUNDARY_CONDITION_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                  CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                    & VALUE,ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_MIXED)
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
       
    CALL EXITS("BOUNDARY_CONDITION_ADD_ELEMENT")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_ADD_ELEMENT",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_ADD_ELEMENT")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_ELEMENT
  
  !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsSetElement
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

    !Note: this routine is for element based interpolation
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
                CASE(BOUNDARY_CONDITION_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                CASE(BOUNDARY_CONDITION_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                    & VALUE,ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_MIXED)
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
 
  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user node. \see OPENCMISS_CMISSBoundaryConditionsAddNode
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_NODE(BOUNDARY_CONDITIONS,VARIABLE_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER, &
    & CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

 
    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_NODE",ERR,ERROR,*999)

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
                CASE(BOUNDARY_CONDITION_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                CASE(BOUNDARY_CONDITION_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                  CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUE, &
                    & ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_MIXED)
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
       
    CALL EXITS("BOUNDARY_CONDITION_ADD_NODE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITION_ADD_NODE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITION_ADD_NODE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_NODE

  !
  !================================================================================================================================
  !

  !>Set the total number of Neumann boundaries 
  SUBROUTINE BOUNDARY_CONDITIONS_BOUNDARY_SET_NUMBER_OF_BOUNDARIES(BOUNDARY_CONDITIONS,VARIABLE_TYPE,NUMBER_OF_NEUMANN,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NEUMANN !<The total number of user-set Neumann boundaries
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    INTEGER(INTG) :: i
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_BOUNDARY_SET_NUMBER_OF_BOUNDARIES",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN

          BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_NEUMANN_BOUNDARIES=NUMBER_OF_NEUMANN

          ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann Boundary Conditions",ERR,ERROR,*999)

          ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%NEUMANN_BOUNDARY_IDENTIFIER(NUMBER_OF_NEUMANN) &
                 & ,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann identifier array",ERR,ERROR,*999)

          DO i=1,NUMBER_OF_NEUMANN
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%NEUMANN_BOUNDARY_IDENTIFIER(i)%PTR &
                 & ,STAT=ERR)
          ENDDO
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann identifier ptr array",ERR,ERROR,*999)

        ELSE
          LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " has not been created."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_BOUNDARY_SET_NUMBER_OF_BOUNDARIES")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_BOUNDARY_SET_NUMBER_OF_BOUNDARIES",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_BOUNDARY_SET_NUMBER_OF_BOUNDARIES")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_BOUNDARY_SET_NUMBER_OF_BOUNDARIES

  !
  !================================================================================================================================
  !

  !>Adds a single specified dof to the specified boundary condition. 
  SUBROUTINE BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF1(BOUNDARY_CONDITIONS,VARIABLE_TYPE,USER_GLOBAL_DOF_NUMBER,VALUE, &
                                                                                             & CONDITION,ID,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: USER_GLOBAL_DOF_NUMBER !<The global dof index to set the boundary condition at
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    INTEGER(INTG), INTENT(IN) :: ID !<Identifier for the Neumann boundary condition
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS(BOUNDARY_CONDITIONS,VARIABLE_TYPE,(/USER_GLOBAL_DOF_NUMBER/),(/VALUE/), &
                                                                                          & (/CONDITION/),ID,ERR,ERROR,*999)
 
    CALL EXITS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF1")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF1",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF1")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOF1

  !
  !================================================================================================================================
  !

  !>Adds a specified dofs to the specified boundary condition. 
  SUBROUTINE BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS(BOUNDARY_CONDITIONS,VARIABLE_TYPE,USER_GLOBAL_DOF_NUMBER,VALUE, &
                                                                                        & CONDITION,ID,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: USER_GLOBAL_DOF_NUMBER(:) !<USER_GLOBAL_DOF_NUMBER(:). The global dof numbers to set the boundary condition at
    REAL(DP), INTENT(IN) :: VALUE(:) !<VALUE(:). The value of the boundary condition for the added nodes
    INTEGER(INTG), INTENT(IN) :: CONDITION(:) !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: ID !<Identifier for the Neumann boundary condition
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: i,NUMBER_OF_VALUES
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_VALUES_TYPE), POINTER :: NEUMANN_VALUES
 
    CALL ENTERS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD           
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
                IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN

                  NUMBER_OF_VALUES=SIZE(USER_GLOBAL_DOF_NUMBER,1)

                  NEUMANN_VALUES=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%NEUMANN_BOUNDARY_IDENTIFIER(ID)%PTR

                  IF (.NOT.ALLOCATED(NEUMANN_VALUES%SET_DOF)) THEN !This is only performed once

                    ALLOCATE(NEUMANN_VALUES%SET_DOF(NUMBER_OF_VALUES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Boundary Conditions dof array.",ERR,ERROR,*999)
                    ALLOCATE(NEUMANN_VALUES%SET_DOF_VALUES(NUMBER_OF_VALUES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Boundary Conditions dof values array",ERR,ERROR,*999)

                    NEUMANN_VALUES%SET_DOF=0
                    NEUMANN_VALUES%SET_DOF_VALUES=0.0_DP

                  ENDIF

                  DO i=1,SIZE(CONDITION,1)
                    SELECT CASE(CONDITION(i))
                    CASE(BOUNDARY_CONDITION_NOT_FIXED)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BOUNDARY_CONDITION_FIXED)
                      BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(USER_GLOBAL_DOF_NUMBER(i))= &
                                                                                         & BOUNDARY_CONDITION_FIXED
                    CASE(BOUNDARY_CONDITION_FIXED_WALL)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BOUNDARY_CONDITION_FIXED_INLET)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BOUNDARY_CONDITION_MOVED_WALL)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BOUNDARY_CONDITION_FREE_WALL)
                      BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(USER_GLOBAL_DOF_NUMBER(i))= &
                                                                                         & BOUNDARY_CONDITION_FREE_WALL
                    CASE(BOUNDARY_CONDITION_NEUMANN)
                      BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(USER_GLOBAL_DOF_NUMBER(i))= &
                                                                                         & BOUNDARY_CONDITION_NEUMANN
                    CASE(BOUNDARY_CONDITION_MIXED)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The specified boundary condition type of "&
                      &//TRIM(NUMBER_TO_VSTRING(CONDITION(i),"*",ERR,ERROR))//" is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT 
                  ENDDO

                  !Add dof number and value to an array on the boundary condition object
                  IF(SIZE(USER_GLOBAL_DOF_NUMBER,1)==SIZE(VALUE,1)) THEN
                    NEUMANN_VALUES%SET_DOF(:)=USER_GLOBAL_DOF_NUMBER(:)
                    NEUMANN_VALUES%SET_DOF_VALUES(:)=VALUE(:)
                  ELSE
                    LOCAL_ERROR="The size of the dof array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(USER_GLOBAL_DOF_NUMBER,1),"*",ERR,ERROR))// &
                      & ") does not match the size of the values array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(VALUE,1),"*",ERR,ERROR))//")."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF

                ELSE
                  CALL FLAG_ERROR("The boundary condition Neumann is not associated",ERR,ERROR,*999)
                ENDIF
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

    CALL EXITS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_BOUNDARY_ADD_DOFS

  !
  !================================================================================================================================
  !

  !>Given listing of faces on boundary calculate value at individual nodes and set boundary condition 
  SUBROUTINE BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE(BOUNDARY_CONDITIONS,VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<Pointer to Topology
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), POINTER :: TEMP_SET_NODES(:),TEMP_FACES(:),TEMP_LINES(:),TEMP_DOFS(:) !Temporary lists of faces, line, nodes and dofs on domain of interest
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_VALUES_TYPE), POINTER :: NEUMANN_VALUES
    INTEGER(INTG) :: TEMP_NUMBER_OF_SET_NODES,TEMP_NUMBER_OF_FACES,TEMP_NUMBER_OF_LINES
    INTEGER(INTG) :: NUMBER_OF_SET_DOF,local_ny,component_idx,nd,nn,nf,nl,nf1,nl1,M,N,i,j,ms,ns
    INTEGER(INTG) :: derivative,id,NUMBER_DOFS_IN_FACE,NUMBER_DOFS_IN_LINE,TOTAL_NUMBER_OF_FACE_DOF
    INTEGER(INTG) :: TOTAL_NUMBER_OF_LINE_DOF,global_ny,MAX_NUMBER_DOFS_IN_FACE,MAX_NUMBER_DOFS_IN_LINE
    INTEGER(INTG) :: FACE_NUMBER,LINE_NUMBER,NODE_NUMBER,x_pos,NODE,face_local_node_index,line_local_node_index
    INTEGER(INTG) :: face_local_deriv_index,line_local_deriv_index
    LOGICAL :: DOF_LOCATED
    TYPE(LIST_TYPE), POINTER :: TEMP_SET_NODES_LIST,TEMP_FACES_LIST,TEMP_LINES_LIST,TEMP_DOFS_LIST
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE",ERR,ERROR,*999)

    NULLIFY(TEMP_SET_NODES)
    NULLIFY(TEMP_FACES)
    NULLIFY(TEMP_LINES)
    NULLIFY(TEMP_DOFS)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS% &
                                      & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
                SELECT CASE(EQUATIONS_SET%TYPE)
                CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE,EQUATIONS_SET_LAPLACE_EQUATION_TYPE)

                  IF(BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_NEUMANN_BOUNDARIES>0) THEN

!Addition needed here to iterate over whole VALUES vector and take these values to produce SET_DOF

                    DO id=1,BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_NEUMANN_BOUNDARIES   !!!Remove this
                      !Iterate over components u,v,w,p
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                        TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY
                        IF(ASSOCIATED(TOPOLOGY)) THEN                        
                          DOMAIN_NODES=>TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            SELECT CASE(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                            CASE(FIELD_NODE_BASED_INTERPOLATION)                      
                              IF(DEPENDENT_FIELD%DECOMPOSITION% &
                                & DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
                                & PTR%NUMBER_OF_DIMENSIONS==3) THEN

                                DOMAIN_FACES=>TOPOLOGY%FACES
                                IF(ASSOCIATED(DOMAIN_FACES)) THEN

                                  TEMP_NUMBER_OF_SET_NODES=0 !Counts the number of nodes that correspond to the set dofs
                                  TEMP_NUMBER_OF_FACES=0 !Counts the number of unique faces of interest
                                  TOTAL_NUMBER_OF_FACE_DOF=0 !Counts the number of unique dof in faces of interest

                                  NEUMANN_VALUES=>BOUNDARY_CONDITIONS_NEUMANN%NEUMANN_BOUNDARY_IDENTIFIER(ID)%PTR

                                  !Calculate number of rows in SET_DOF
                                  NUMBER_OF_SET_DOF=SIZE(NEUMANN_VALUES%SET_DOF,1) 
                                  !This includes all set dof for all components and derivatives

                                  !The following two lists are for calculating the size of the INTEGRATION_MATRIX
                                  !Create listing of nodes on which DOF have been set
                                  NULLIFY(TEMP_SET_NODES_LIST)
                                  CALL LIST_CREATE_START(TEMP_SET_NODES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DATA_TYPE_SET(TEMP_SET_NODES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                                  CALL LIST_INITIAL_SIZE_SET(TEMP_SET_NODES_LIST,NUMBER_OF_SET_DOF,ERR,ERROR,*999)
                                  CALL LIST_CREATE_FINISH(TEMP_SET_NODES_LIST,ERR,ERROR,*999)

                                  DO nd=1,NUMBER_OF_SET_DOF
                                    !Check if dof in current component (uvwp)
                                    IF(FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                                      & NODE_DOF2PARAM_MAP(3,NEUMANN_VALUES%SET_DOF(nd))==component_idx) THEN
                                      !Convert dof to node
                                      CALL LIST_ITEM_ADD(TEMP_SET_NODES_LIST,FIELD_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP &
                                        & (2,NEUMANN_VALUES%SET_DOF(nd)),ERR,ERROR,*999)
                                    ENDIF
                                  ENDDO !nd

                                  CALL LIST_REMOVE_DUPLICATES(TEMP_SET_NODES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DETACH_AND_DESTROY(TEMP_SET_NODES_LIST, &
                                                  & TEMP_NUMBER_OF_SET_NODES,TEMP_SET_NODES,ERR,ERROR,*999)

                                  !Create list of surrounding faces
                                  NULLIFY(TEMP_FACES_LIST)
                                  CALL LIST_CREATE_START(TEMP_FACES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DATA_TYPE_SET(TEMP_FACES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                                  CALL LIST_INITIAL_SIZE_SET(TEMP_FACES_LIST,DOMAIN_FACES%NUMBER_OF_FACES,ERR,ERROR,*999)
                                  CALL LIST_CREATE_FINISH(TEMP_FACES_LIST,ERR,ERROR,*999)

                                  MAX_NUMBER_DOFS_IN_FACE=0
                                  !Find the faces this node is a part of, by looping over set nodes that correspond to set dof
                                  DO nn=1,TEMP_NUMBER_OF_SET_NODES
                                    DOMAIN_NODE=>DOMAIN_NODES%NODES(TEMP_SET_NODES(nn))
                                    DO face_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_FACES

                                      DO nf1=1,DOMAIN_FACES%NUMBER_OF_FACES
                                        DOMAIN_FACE=>DOMAIN_FACES%FACES(nf1)
                                        IF(DOMAIN_FACE%NUMBER==DOMAIN_NODE%NODE_FACES(face_local_node_index)) THEN
                                          IF(DOMAIN_FACE%BOUNDARY_FACE) THEN
                                            CALL LIST_ITEM_ADD(TEMP_FACES_LIST,DOMAIN_NODE% &
                                                                        & NODE_FACES(face_local_node_index),ERR,ERROR,*999)
                                            MAX_NUMBER_DOFS_IN_FACE=MAX(MAX_NUMBER_DOFS_IN_FACE,DOMAIN_FACE% & 
                                                                                    & BASIS%NUMBER_OF_ELEMENT_PARAMETERS)
                                          ENDIF
                                        ENDIF
                                      ENDDO !nf1
                                    ENDDO !face_local_node_index
                                  ENDDO !nn

                                  CALL LIST_REMOVE_DUPLICATES(TEMP_FACES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DETACH_AND_DESTROY(TEMP_FACES_LIST,TEMP_NUMBER_OF_FACES,TEMP_FACES,ERR,ERROR,*999)

                                  !Use a list to collect all the dof, sort and substitute into INTEGRATION_MATRIX_MAPPING array 
                                  NULLIFY(TEMP_DOFS_LIST)
                                  CALL LIST_CREATE_START(TEMP_DOFS_LIST,ERR,ERROR,*999)
                                  CALL LIST_DATA_TYPE_SET(TEMP_DOFS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                                  CALL LIST_INITIAL_SIZE_SET(TEMP_DOFS_LIST,TEMP_NUMBER_OF_FACES* &
                                                                                      & MAX_NUMBER_DOFS_IN_FACE,ERR,ERROR,*999)
                                  CALL LIST_CREATE_FINISH(TEMP_DOFS_LIST,ERR,ERROR,*999)


                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF &
                                                              & (DOMAIN_FACES%NUMBER_OF_FACES,MAX_NUMBER_DOFS_IN_FACE))

                                  BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF = 0

                                  !Find the faces this node is a part of, by looping over set nodes that correspond to set dof (global dofs)
                                  DO nf=1,TEMP_NUMBER_OF_FACES
                                    DO nf1=1,DOMAIN_FACES%NUMBER_OF_FACES
                                      DOMAIN_FACE=>DOMAIN_FACES%FACES(nf1)
                                      IF(DOMAIN_FACE%NUMBER==TEMP_FACES(nf).AND.DOMAIN_FACE%BOUNDARY_FACE) THEN
                                       
                                        DO face_local_node_index=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES !nnf

                                          NODE=DOMAIN_FACE%NODES_IN_FACE(face_local_node_index)
                                          DO face_local_deriv_index=1,DOMAIN_FACE%BASIS% &
                                                                      & NUMBER_OF_DERIVATIVES(face_local_node_index) !nnd

                                            DERIVATIVE=DOMAIN_FACE%DERIVATIVES_IN_FACE(face_local_deriv_index, &
                                                                                     & face_local_node_index)
                                            local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)% &
                                                                   & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(DERIVATIVE,NODE)
                                            x_pos = DOMAIN_FACE%BASIS%ELEMENT_PARAMETER_INDEX(face_local_deriv_index, &
                                                                                            & face_local_node_index)

                                            CALL LIST_ITEM_ADD(TEMP_DOFS_LIST,local_ny,ERR,ERROR,*999)

                                            !The mapping of local_ny to element parameter number (ms or ns) is stored in a series of vectors
                                            !that represent the mapping for each face. FACES_ELEMENT_PARAM_2_LOCAL_DOF is indexed by face number
                                            !and element parameter number; corresponding local dof number is stored.

                                            BOUNDARY_CONDITIONS_NEUMANN% &
                                                          & FACES_ELEMENT_PARAM_2_LOCAL_DOF(DOMAIN_FACE%NUMBER,x_pos) = local_ny
                                          ENDDO !face_local_deriv_index
                                        ENDDO !face_local_node_index
                                      ENDIF
                                    ENDDO !nf1
                                  ENDDO !nf

                                  CALL LIST_REMOVE_DUPLICATES(TEMP_DOFS_LIST,ERR,ERROR,*999)
                                  CALL LIST_DETACH_AND_DESTROY(TEMP_DOFS_LIST,TOTAL_NUMBER_OF_FACE_DOF,TEMP_DOFS,ERR,ERROR,*999)

                                  !Number of DOF involved in boundary condition calculation
                                  M = TOTAL_NUMBER_OF_FACE_DOF
                                  !Number of DOF where boundary conditions have been fixed by user 
                                  N = NUMBER_OF_SET_DOF

                                  !Integration_matrix = 'A' matrix
                                  !Point_values_vector = 'x' vector
                                  !Integrated_values_vector = 'B' vector

                                  !Set size of matrix and vectors
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(MAX_NUMBER_DOFS_IN_FACE),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann face integration matrix array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING(MAX_NUMBER_DOFS_IN_FACE &
                                                                                                                & ),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann face integration matrix mapping array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(M,N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix array." &
                                                                                                            & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X(N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix mapping X array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y(M),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix mapping Y array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR(N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann point values vector.",ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING(N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann point values vector mapping."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(M),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integrated values vector." &
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(M),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integrated values vector mapping."&
                                                                                                             & ,ERR,ERROR,*999)

                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX=0.0_DP
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X=0
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y=0
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR=0.0_DP
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING=0
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR=0.0_DP
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING=0

                                  !Calculate the mapping of INTEGRATION_MATRIX, POINT_VALUES_VECTOR and INTEGRATED_VALUES_VECTOR
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X = NEUMANN_VALUES%SET_DOF
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING = NEUMANN_VALUES%SET_DOF
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y=TEMP_DOFS(1:TOTAL_NUMBER_OF_FACE_DOF)
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING= &
                                                                                         & TEMP_DOFS(1:TOTAL_NUMBER_OF_FACE_DOF)
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR = NEUMANN_VALUES%SET_DOF_VALUES


                                  !Iterate over set dofs
                                  DO nd=1,NUMBER_OF_SET_DOF
                                    !Check if dof in current component (uvwp)
                                    IF(FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                                      & NODE_DOF2PARAM_MAP(3,NEUMANN_VALUES%SET_DOF(nd))==component_idx) THEN

                                      !Convert dof to node number
                                      NODE_NUMBER=FIELD_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(2,NEUMANN_VALUES%SET_DOF(nd))
                                      DOMAIN_NODE=>DOMAIN_NODES%NODES(NODE_NUMBER)

                                      !Iterate over the faces on which that node is a part of
                                      DO face_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_FACES

                                        DO nf1=1,DOMAIN_FACES%NUMBER_OF_FACES
                                          DOMAIN_FACE=>DOMAIN_FACES%FACES(nf1)
                                          ns = 0

                                          !Locate face and ensure face is a boundary face
                                          IF(DOMAIN_FACE%NUMBER==DOMAIN_NODE%NODE_FACES(face_local_node_index) &
                                            & .AND.(DOMAIN_FACE%BOUNDARY_FACE)) THEN

                                            FACE_NUMBER = DOMAIN_NODE%NODE_FACES(face_local_node_index)

                                            !Calculate the ns corresponding to the dof
                                            DO ms=1,DOMAIN_FACE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                                              local_ny=BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF(FACE_NUMBER,ms)
                                              IF(local_ny==NEUMANN_VALUES%SET_DOF(nd)) THEN
                                                ns = ms
                                              ENDIF
                                            ENDDO !ms

                                            !Calculate basis for face
                                            CALL BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,VARIABLE_TYPE, &
                                              & component_idx,FACE_NUMBER,NUMBER_DOFS_IN_FACE,ns,ERR,ERROR,*999)

                                            !Iterate over the number of dofs in the face
                                            DO ms=1,NUMBER_DOFS_IN_FACE 
                                              DOF_LOCATED=.FALSE.

                                              !Locate in INTEGRATION_MATRIX and add into, one column for each SET_DOF
                                              DO j=1,TOTAL_NUMBER_OF_FACE_DOF
                                                IF((DOF_LOCATED.NEQV..TRUE.).AND. &
                                                  &(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y(j) == &
                                                  & BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING(ms))) THEN

                                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,nd)= &
                                                    & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,nd) + & 
                                                    & BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms)

                                                  DOF_LOCATED=.TRUE.
                                                ENDIF
                                              ENDDO !j
                                            ENDDO !ms

                                          ENDIF
                                        ENDDO !nf1
                                      ENDDO !face_local_node_index
                                    ENDIF
                                  ENDDO !nd

                                  !Perform Ax=B calculation
                                  DO j=1,TOTAL_NUMBER_OF_FACE_DOF
                                    DO i=1,NUMBER_OF_SET_DOF
                                      BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) = &
                                        & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) + &
                                        & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,i) * &
                                        & BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR(i)
                                    ENDDO !i
!!Check here that global_ny is the best thing to do here
                                    global_ny=BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)
                                    !Ensure the set dof are set as BOUNDARY_CONDITION_NEUMANN
                                    BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NEUMANN
                                  ENDDO !j

                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_SIZE=TOTAL_NUMBER_OF_FACE_DOF

                                  DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%NEUMANN_BOUNDARY_IDENTIFIER)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING)
                                  DEALLOCATE(TEMP_SET_NODES)
                                  DEALLOCATE(TEMP_FACES)

                                ELSE
                                  CALL FLAG_ERROR("Topology faces is not associated",ERR,ERROR,*999)
                                ENDIF

                              ELSEIF(DEPENDENT_FIELD%DECOMPOSITION% &
                                & DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
                                & PTR%NUMBER_OF_DIMENSIONS==2) THEN

                                DOMAIN_LINES=>TOPOLOGY%LINES
                                IF(ASSOCIATED(DOMAIN_LINES)) THEN

                                  TEMP_NUMBER_OF_SET_NODES=0 !Counts the number of nodes that correspond to the set dofs
                                  TEMP_NUMBER_OF_LINES=0 !Counts the number of unique lines of interest
                                  TOTAL_NUMBER_OF_LINE_DOF=0 !Counts the number of unique dof in lines of interest

                                  NEUMANN_VALUES=>BOUNDARY_CONDITIONS_NEUMANN%NEUMANN_BOUNDARY_IDENTIFIER(ID)%PTR

                                  !Calculate number of rows in SET_DOF
                                  NUMBER_OF_SET_DOF=SIZE(NEUMANN_VALUES%SET_DOF,1) 
                                  !This includes all set dof for all components and derivatives

                                  !The following two lists are for calculating the size of the INTEGRATION_MATRIX
                                  !Create listing of nodes on which DOF have been set
                                  NULLIFY(TEMP_SET_NODES_LIST)
                                  CALL LIST_CREATE_START(TEMP_SET_NODES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DATA_TYPE_SET(TEMP_SET_NODES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                                  CALL LIST_INITIAL_SIZE_SET(TEMP_SET_NODES_LIST,NUMBER_OF_SET_DOF,ERR,ERROR,*999)
                                  CALL LIST_CREATE_FINISH(TEMP_SET_NODES_LIST,ERR,ERROR,*999)

                                  DO nd=1,NUMBER_OF_SET_DOF
                                    !Check if dof in current component (uvwp)
                                    IF(FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                                      & NODE_DOF2PARAM_MAP(3,NEUMANN_VALUES%SET_DOF(nd))==component_idx) THEN
                                      !Convert dof to node
                                      CALL LIST_ITEM_ADD(TEMP_SET_NODES_LIST,FIELD_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP &
                                        & (2,NEUMANN_VALUES%SET_DOF(nd)),ERR,ERROR,*999)
                                    ENDIF
                                  ENDDO !nd

                                  CALL LIST_REMOVE_DUPLICATES(TEMP_SET_NODES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DETACH_AND_DESTROY(TEMP_SET_NODES_LIST, &
                                                  & TEMP_NUMBER_OF_SET_NODES,TEMP_SET_NODES,ERR,ERROR,*999)

                                  !Create list of surrounding lines
                                  NULLIFY(TEMP_LINES_LIST)
                                  CALL LIST_CREATE_START(TEMP_LINES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DATA_TYPE_SET(TEMP_LINES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                                  CALL LIST_INITIAL_SIZE_SET(TEMP_LINES_LIST,DOMAIN_LINES%NUMBER_OF_LINES,ERR,ERROR,*999)
                                  CALL LIST_CREATE_FINISH(TEMP_LINES_LIST,ERR,ERROR,*999)

                                  MAX_NUMBER_DOFS_IN_LINE=0
                                  !Find the lines this node is a part of, by looping over set nodes that correspond to set dof
                                  DO nn=1,TEMP_NUMBER_OF_SET_NODES
                                    DOMAIN_NODE=>DOMAIN_NODES%NODES(TEMP_SET_NODES(nn))
                                    DO line_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_LINES

                                      DO nl1=1,DOMAIN_LINES%NUMBER_OF_LINES
                                        DOMAIN_LINE=>DOMAIN_LINES%LINES(nl1)
                                        IF(DOMAIN_LINE%NUMBER==DOMAIN_NODE%NODE_LINES(line_local_node_index)) THEN
                                          IF(DOMAIN_LINE%BOUNDARY_LINE) THEN
                                            CALL LIST_ITEM_ADD(TEMP_LINES_LIST,DOMAIN_NODE% &
                                                                        & NODE_LINES(line_local_node_index),ERR,ERROR,*999)
                                            MAX_NUMBER_DOFS_IN_LINE=MAX(MAX_NUMBER_DOFS_IN_LINE,DOMAIN_LINE% & 
                                                                                    & BASIS%NUMBER_OF_ELEMENT_PARAMETERS)
                                          ENDIF
                                        ENDIF
                                      ENDDO !nl1
                                    ENDDO !line_local_node_index
                                  ENDDO !nn

                                  CALL LIST_REMOVE_DUPLICATES(TEMP_LINES_LIST,ERR,ERROR,*999)
                                  CALL LIST_DETACH_AND_DESTROY(TEMP_LINES_LIST,TEMP_NUMBER_OF_LINES,TEMP_LINES,ERR,ERROR,*999)

                                  !Use a list to collect all the dof, sort and substitute into INTEGRATION_MATRIX_MAPPING array 
                                  NULLIFY(TEMP_DOFS_LIST)
                                  CALL LIST_CREATE_START(TEMP_DOFS_LIST,ERR,ERROR,*999)
                                  CALL LIST_DATA_TYPE_SET(TEMP_DOFS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                                  CALL LIST_INITIAL_SIZE_SET(TEMP_DOFS_LIST,TEMP_NUMBER_OF_LINES* &
                                                                                      & MAX_NUMBER_DOFS_IN_LINE,ERR,ERROR,*999)
                                  CALL LIST_CREATE_FINISH(TEMP_DOFS_LIST,ERR,ERROR,*999)


                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF &
                                                              & (DOMAIN_LINES%NUMBER_OF_LINES,MAX_NUMBER_DOFS_IN_LINE))

                                  BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF = 0

                                  !Find the lines this node is a part of, by looping over set nodes that correspond to set dof (global dofs)
                                  DO nl=1,TEMP_NUMBER_OF_LINES
                                    DO nl1=1,DOMAIN_LINES%NUMBER_OF_LINES
                                      DOMAIN_LINE=>DOMAIN_LINES%LINES(nl1)
                                      IF(DOMAIN_LINE%NUMBER==TEMP_LINES(nl).AND.DOMAIN_LINE%BOUNDARY_LINE) THEN
                                       
                                        DO line_local_node_index=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES !nnl

                                          NODE=DOMAIN_LINE%NODES_IN_LINE(line_local_node_index)
                                          DO line_local_deriv_index=1,DOMAIN_LINE%BASIS% &
                                                                      & NUMBER_OF_DERIVATIVES(line_local_node_index) !nnd

                                            DERIVATIVE=DOMAIN_LINE%DERIVATIVES_IN_LINE(line_local_deriv_index, &
                                                                                     & line_local_node_index)
                                            local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)% &
                                                                   & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(DERIVATIVE,NODE)
                                            x_pos = DOMAIN_LINE%BASIS%ELEMENT_PARAMETER_INDEX(line_local_deriv_index, &
                                                                                            & line_local_node_index)

                                            CALL LIST_ITEM_ADD(TEMP_DOFS_LIST,local_ny,ERR,ERROR,*999)

                                            !The mapping of local_ny to element parameter number (ms or ns) is stored in a series of vectors
                                            !that represent the mapping for each line. LINES_ELEMENT_PARAM_2_LOCAL_DOF is indexed by line number
                                            !and element parameter number; corresponding local dof number is stored.

                                            BOUNDARY_CONDITIONS_NEUMANN% &
                                                          & LINES_ELEMENT_PARAM_2_LOCAL_DOF(DOMAIN_LINE%NUMBER,x_pos) = local_ny
                                          ENDDO !line_local_deriv_index
                                        ENDDO !line_local_node_index
                                      ENDIF
                                    ENDDO !nl1
                                  ENDDO !nl

                                  CALL LIST_REMOVE_DUPLICATES(TEMP_DOFS_LIST,ERR,ERROR,*999)
                                  CALL LIST_DETACH_AND_DESTROY(TEMP_DOFS_LIST,TOTAL_NUMBER_OF_LINE_DOF,TEMP_DOFS,ERR,ERROR,*999)

                                  !Number of DOF involved in boundary condition calculation
                                  M = TOTAL_NUMBER_OF_LINE_DOF
                                  !Number of DOF where boundary conditions have been fixed by user 
                                  N = NUMBER_OF_SET_DOF

                                  !Integration_matrix = 'A' matrix
                                  !Point_values_vector = 'x' vector
                                  !Integrated_values_vector = 'B' vector

                                  !Set size of matrix and vectors
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX(MAX_NUMBER_DOFS_IN_LINE),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann line integration matrix array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX_MAPPING(MAX_NUMBER_DOFS_IN_LINE &
                                                                                                                & ),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann line integration matrix mapping array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(M,N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix array." &
                                                                                                            & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X(N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix mapping X array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y(M),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix mapping Y array."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR(N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann point values vector.",ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING(N),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann point values vector mapping."&
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(M),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integrated values vector." &
                                                                                                             & ,ERR,ERROR,*999)
                                  ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(M),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann integrated values vector mapping."&
                                                                                                             & ,ERR,ERROR,*999)

                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX=0.0_DP
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X=0
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y=0
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR=0.0_DP
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING=0
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR=0.0_DP
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING=0

                                  !Calculate the mapping of INTEGRATION_MATRIX, POINT_VALUES_VECTOR and INTEGRATED_VALUES_VECTOR
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X = NEUMANN_VALUES%SET_DOF
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING = NEUMANN_VALUES%SET_DOF
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y=TEMP_DOFS(1:TOTAL_NUMBER_OF_LINE_DOF)
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING= &
                                                                                         & TEMP_DOFS(1:TOTAL_NUMBER_OF_LINE_DOF)
                                  BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR = NEUMANN_VALUES%SET_DOF_VALUES


                                  !Iterate over set dofs
                                  DO nd=1,NUMBER_OF_SET_DOF
                                    !Check if dof in current component (uvwp)
                                    IF(FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                                      & NODE_DOF2PARAM_MAP(3,NEUMANN_VALUES%SET_DOF(nd))==component_idx) THEN

                                      !Convert dof to node number
                                      NODE_NUMBER=FIELD_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(2,NEUMANN_VALUES%SET_DOF(nd))
                                      DOMAIN_NODE=>DOMAIN_NODES%NODES(NODE_NUMBER)

                                      !Iterate over the lines on which that node is a part of
                                      DO line_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_LINES

                                        DO nl1=1,DOMAIN_LINES%NUMBER_OF_LINES
                                          DOMAIN_LINE=>DOMAIN_LINES%LINES(nl1)
                                          ns = 0

                                          !Locate line and ensure line is a boundary line
                                          IF(DOMAIN_LINE%NUMBER==DOMAIN_NODE%NODE_LINES(line_local_node_index) &
                                            & .AND.(DOMAIN_LINE%BOUNDARY_LINE)) THEN

                                            LINE_NUMBER = DOMAIN_NODE%NODE_LINES(line_local_node_index)

                                            !Calculate the ns corresponding to the dof
                                            DO ms=1,DOMAIN_LINE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                                              local_ny=BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF(LINE_NUMBER,ms)
                                              IF(local_ny==NEUMANN_VALUES%SET_DOF(nd)) THEN
                                                ns = ms
                                              ENDIF
                                            ENDDO !ms


                                            !Calculate basis for line
                                            CALL BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,VARIABLE_TYPE, &
                                              & component_idx,LINE_NUMBER,NUMBER_DOFS_IN_LINE,ns,ERR,ERROR,*999)

                                            !Iterate over the number of dofs in the line
                                            DO ms=1,NUMBER_DOFS_IN_LINE 
                                              DOF_LOCATED=.FALSE.

                                              !Locate in INTEGRATION_MATRIX and add into, one column for each SET_DOF
                                              DO j=1,TOTAL_NUMBER_OF_LINE_DOF
                                                IF((DOF_LOCATED.NEQV..TRUE.).AND. &
                                                  &(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y(j) == &
                                                  & BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX_MAPPING(ms))) THEN

                                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,nd)= &
                                                    & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,nd) + & 
                                                    & BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX(ms)

                                                  DOF_LOCATED=.TRUE.
                                                ENDIF
                                              ENDDO !j
                                            ENDDO !ms

                                          ENDIF
                                        ENDDO !nl1
                                      ENDDO !line_local_node_index
                                    ENDIF
                                  ENDDO !nd

                                  !Perform Ax=B calculation
                                  DO j=1,TOTAL_NUMBER_OF_LINE_DOF
                                    DO i=1,NUMBER_OF_SET_DOF
                                      BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) = &
                                        & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) + &
                                        & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,i) * &
                                        & BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR(i)
                                    ENDDO !i
!!Check here that global_ny is the best thing to do here
                                    global_ny=BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)
                                    !Ensure the set dof are set as BOUNDARY_CONDITION_NEUMANN
                                    BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NEUMANN
                                  ENDDO !j

                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_SIZE=TOTAL_NUMBER_OF_LINE_DOF

                                  DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%NEUMANN_BOUNDARY_IDENTIFIER)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX_MAPPING)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR)
                                  DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING)
                                  DEALLOCATE(TEMP_SET_NODES)
                                  DEALLOCATE(TEMP_LINES)
                                ELSE
                                  CALL FLAG_ERROR("Topology lines is not associated",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Invalid number of dimensions for Neumann boundary conditions.",ERR,ERROR,*999)
                              ENDIF
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE DEFAULT
                              LOCAL_ERROR="The interpolation type of "// &
                                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(component_idx)&
                                & %INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                                & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                                & " of dependent variable number "// &
                                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))//"."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                            END SELECT                 
                          ELSE
                            CALL FLAG_ERROR("Topology nodes is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
                        ENDIF
                      ENDDO !component_idx
                    ENDDO !id
                  ENDIF
                CASE DEFAULT
                  LOCAL_ERROR="Equations set equation type "//&
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
                    & " is not a valid equations set type for the integrated boundary calculation."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)  
                END SELECT
              ELSE
                CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The boundary condition neumann is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type "&
                      &//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " has not been created."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE

  !
  !================================================================================================================================
  !

  !>For a given face, calculates the contributions of the basis functions and integrated flux to RHS
  SUBROUTINE BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,VARIABLE_TYPE,component_idx,FACE_NUMBER, &
                                                                                       & NUMBER_DOFS_IN_FACE,ns,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: component_idx !<The component number, i.e. one of u,v,w
    INTEGER(INTG), INTENT(IN) :: FACE_NUMBER !<The face number (the OpenCMISS face number)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_DOFS_IN_FACE !<The number of dofs in the face
    INTEGER(INTG), INTENT(IN) :: ns !<The element parameter number of the set dof
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the calculations on
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ELEMENT_NUMBER,ng,ms,local_ny
    REAL(DP) :: RWG,PHIMS,PHINS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS% &
                         & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
              DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
              TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY
              IF(ASSOCIATED(TOPOLOGY)) THEN
                DOMAIN_FACES=>TOPOLOGY%FACES
                IF(ASSOCIATED(DOMAIN_FACES)) THEN
                  DOMAIN_FACE=>DOMAIN_FACES%FACES(FACE_NUMBER)
                  IF(ASSOCIATED(DOMAIN_FACE)) THEN
 
                    ELEMENT_NUMBER=DOMAIN_FACE%ELEMENT_NUMBER 

                    QUADRATURE_SCHEME=>DOMAIN_FACE%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                    IF(ASSOCIATED(QUADRATURE_SCHEME)) THEN
 
                      SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_IDX)%INTERPOLATION_TYPE)
                      CASE(FIELD_NODE_BASED_INTERPOLATION)

                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
                          & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                        BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX=0.0_DP
                        BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING=0

                        !Loop over gauss points
                        DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS

                          CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                          CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
                            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                          !Calculate RWG
                          RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN &
                                & *QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

                          !ns is singular because we are doing an integration of the ns with all ms in the face
                          PHINS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          !Loop over element rows
                          DO ms=1,DOMAIN_FACE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS

                            PHIMS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms) = &
                              & BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms) + PHIMS*PHINS*RWG

                            local_ny = BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF(FACE_NUMBER,ms)
                            BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING(ms) = local_ny

                          ENDDO !ms
                        ENDDO !ng
 
                        NUMBER_DOFS_IN_FACE=DOMAIN_FACE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS

                      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(FIELD_CONSTANT_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The interpolation type of "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                          & " is invalid for component number "// &
                          & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                          & " of dependent variable number "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Quadrature scheme is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain topology face is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain topology faces is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain topology is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The boundary condition neumann is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type " &
                    & //TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " has not been created."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE

  !
  !================================================================================================================================
  !

  !>For a given line, calculates the contributions of the basis functions and integrated flux to RHS
  SUBROUTINE BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,VARIABLE_TYPE,component_idx,LINE_NUMBER, &
                                                                                       & NUMBER_DOFS_IN_LINE,ns,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: component_idx !<The component number, i.e. one of u,v,w
    INTEGER(INTG), INTENT(IN) :: LINE_NUMBER !<The line number (the OpenCMISS line number)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_DOFS_IN_LINE !<The number of dofs in the line
    INTEGER(INTG), INTENT(IN) :: ns !<The element parameter number of the set dof
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the calculations on
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ELEMENT_NUMBER,ng,ms,local_ny
    REAL(DP) :: RWG,PHIMS,PHINS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      EQUATIONS_SET=>BOUNDARY_CONDITIONS%EQUATIONS_SET
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS% &
                         & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
              DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
              TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY
              IF(ASSOCIATED(TOPOLOGY)) THEN
                DOMAIN_LINES=>TOPOLOGY%LINES
                IF(ASSOCIATED(DOMAIN_LINES)) THEN
                  DOMAIN_LINE=>DOMAIN_LINES%LINES(LINE_NUMBER)
                  IF(ASSOCIATED(DOMAIN_LINE)) THEN
 
                    ELEMENT_NUMBER=DOMAIN_LINE%ELEMENT_NUMBER 

                    QUADRATURE_SCHEME=>DOMAIN_LINE%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

                    IF(ASSOCIATED(QUADRATURE_SCHEME)) THEN

                      SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_IDX)%INTERPOLATION_TYPE)
                      CASE(FIELD_NODE_BASED_INTERPOLATION)

                        CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,LINE_NUMBER, &
                          & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                        BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX=0.0_DP
                        BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX_MAPPING=0

                        !Loop over gauss points
                        DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS

                          CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                          CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
                            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                          !Calculate RWG
                          RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN &
                                & *QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

                          !ns is singular because we are doing an integration of the ns with all ms in the line
                          PHINS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          !Loop over element rows
                          DO ms=1,DOMAIN_LINE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS

                            PHIMS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX(ms) = &
                              & BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX(ms) + PHIMS*PHINS*RWG

                            local_ny = BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF(LINE_NUMBER,ms)
                            BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX_MAPPING(ms) = local_ny

                          ENDDO !ms
                        ENDDO !ng

                        NUMBER_DOFS_IN_LINE=DOMAIN_LINE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS

                      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(FIELD_CONSTANT_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The interpolation type of "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                          & " is invalid for component number "// &
                          & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                          & " of dependent variable number "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Quadrature scheme is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain topology line is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain topology lines is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain topology is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The boundary condition Neumann is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type " &
                    & //TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " has not been created."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user node. \see OPENCMISS_CMISSBoundaryConditionsSetNode
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

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(DEPENDENT_FIELD)

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
                CASE(BOUNDARY_CONDITION_NOT_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_NOT_FIXED
                CASE(BOUNDARY_CONDITION_FIXED)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUE, &
                    & ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_FIXED_WALL)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED_WALL
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUE, &
                    & ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_FIXED_INLET)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED_INLET
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUE, &
                    & ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_MOVED_WALL)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_MOVED_WALL
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUE, &
                    & ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_FREE_WALL)
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FREE_WALL
                CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED_TYPE) !For load increment loops
                  BOUNDARY_CONDITIONS_VARIABLE%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=BOUNDARY_CONDITION_FIXED_INCREMENTED_TYPE
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                    & local_ny,VALUE,ERR,ERROR,*999)
                CASE(BOUNDARY_CONDITION_MIXED)
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
!!      IF(ALLOCATED(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%NEUMANN_BOUNDARY_IDENTIFIER)) &
!!        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%NEUMANN_BOUNDARY_IDENTIFIER) !NOTE: SET_DOF not deallocated
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
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to initialise a variable type for.
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
              NULLIFY(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%DIRICHLET_BOUNDARY_CONDITIONS)
              BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%NUMBER_OF_DIRICHLET_CONDITIONS=0
              NULLIFY(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%NEUMANN_BOUNDARY_CONDITIONS)
              BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(variable_type)%PTR%NUMBER_OF_NEUMANN_BOUNDARIES=0
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

  !>Gets the boundary conditions for an equations set. \see OPENCMISS_CMISSEquationsSetBoundaryConditionsGet
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

  !>Initialise dirichlet boundary conditions for a boundary conditions.
  SUBROUTINE BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE !<A pointer to the boundary conditions variable to initialise a boundary conditions dirichlet type for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_DIRICHLET_CONDITIONS, NUMBER_OF_LINEAR_MATRICES, NUMBER_OF_DYNAMIC_MATRICES, matrix_idx
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: BOUNDARY_CONDITIONS_DIRICHLET
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING

    CALL ENTERS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS)) THEN
        CALL FLAG_ERROR("Dirichlet boundary conditions are already associated for this boundary conditions variable." &
           & ,ERR,ERROR,*998)
      ELSE
        ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Dirichlet Boundary Conditions",ERR,ERROR,*999)
        BOUNDARY_CONDITIONS_DIRICHLET=>BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS
        NUMBER_OF_DIRICHLET_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
        ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(NUMBER_OF_DIRICHLET_CONDITIONS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Dirichlet DOF indices array",ERR,ERROR,*999)

        EQUATIONS_SET=>BOUNDARY_CONDITIONS_VARIABLE%BOUNDARY_CONDITIONS%EQUATIONS_SET
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          EQUATIONS=>EQUATIONS_SET%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
            IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
              IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                NUMBER_OF_LINEAR_MATRICES=LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES(NUMBER_OF_LINEAR_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Dirichlet linear sparsity indices array",ERR,ERROR,*999)
                DO matrix_idx=1,NUMBER_OF_LINEAR_MATRICES
                  NULLIFY(BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES(matrix_idx)%PTR)
                ENDDO
              ENDIF
              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                NUMBER_OF_DYNAMIC_MATRICES=DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES(NUMBER_OF_DYNAMIC_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Dirichlet dynamic sparsity indices array",ERR,ERROR,*999)
                DO matrix_idx=1,NUMBER_OF_DYNAMIC_MATRICES
                  NULLIFY(BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES(matrix_idx)%PTR)
                ENDDO
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions variable is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE")
    RETURN

999 CALL ERRORS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE",ERR,ERROR)
!!TODO \todo write BOUNDARY_CONDITIONS_DIRICHLET_FINALISE
998 CALL ERRORS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialise Sparsity Indices type
  SUBROUTINE BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE(SPARSITY_INDICES,NUMBER_OF_DIRICHLET,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES !<A pointer to the Sparsity Indices type tp initialise
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIRICHLET !<The number of dirichlet conditions this sparsity indices type will hold
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SPARSITY_INDICES)) THEN
     CALL FLAG_ERROR("Sparsity Indices are already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(SPARSITY_INDICES,STAT=ERR)
      NULLIFY(SPARSITY_INDICES%SPARSE_ROW_INDICES)
      ALLOCATE(SPARSITY_INDICES%SPARSE_COLUMN_INDICES(NUMBER_OF_DIRICHLET+1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sparsity column indices array",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE",ERR,ERROR)
!!TODO \todo write BOUNDARY_CONDITIONS_SPARSITY_INDICES_FINALISE
998 CALL ERRORS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE BOUNDARY_CONDITIONS_ROUTINES
