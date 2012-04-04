!> \file 
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
  USE LINKEDLIST_ROUTINES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_DOFTypes BOUNDARY_CONDITIONS_ROUTINES::DOFTypes
  !> \brief DOF type for boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FREE=0 !<The dof is free. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FIXED=1 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_MIXED=2 !<The dof is set as a mixed boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  !>@}
  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions BOUNDARY_CONDITIONS_ROUTINES::BoundaryConditions
  !> \brief Boundary conditions types. These may be specific to a particular equation type and the solver routines should not need to use these.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE=0 !<The dof is free. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED=1 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INLET=2 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_OUTLET=3 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_WALL=4 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL=5 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE_WALL=6 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT=8 !<The dof is set to a Neumann point boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED=9 !<The dof is set to a Neumann integrated boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DIRICHLET=10 !<The dof is set to a Dirichlet boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CAUCHY=11 !<The dof is set to a Cauchy boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_ROBIN=12 !<The dof is set to a Robin boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INCREMENTED=13 !<The dof is a fixed boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE=14 !<The dof is a surface pressure boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE_INCREMENTED=15 !<The dof is a surface pressure boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED=17 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE=18 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_IMPERMEABLE_WALL=19 !<The dof is set such that (via penalty formulation): velocity * normal = 0. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  !>@}

  INTEGER(INTG), PARAMETER :: MAX_BOUNDARY_CONDITION_NUMBER=19 !The maximum boundary condition type identifier, used for allocating an array with an entry for each type
  
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

  PUBLIC BOUNDARY_CONDITION_DOF_FREE,BOUNDARY_CONDITION_DOF_FIXED,BOUNDARY_CONDITION_DOF_MIXED

  PUBLIC BOUNDARY_CONDITION_FREE,BOUNDARY_CONDITION_FIXED,BOUNDARY_CONDITION_FIXED_INLET,&
    & BOUNDARY_CONDITION_FIXED_OUTLET,BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL,BOUNDARY_CONDITION_FREE_WALL,&
    & BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_DIRICHLET,BOUNDARY_CONDITION_NEUMANN_POINT, &
    & BOUNDARY_CONDITION_CAUCHY,BOUNDARY_CONDITION_ROBIN,BOUNDARY_CONDITION_FIXED_INCREMENTED,BOUNDARY_CONDITION_PRESSURE,&
    & BOUNDARY_CONDITION_PRESSURE_INCREMENTED,BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED, &
    & BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE,BOUNDARY_CONDITION_IMPERMEABLE_WALL

  PUBLIC BOUNDARY_CONDITIONS_CREATE_FINISH,BOUNDARY_CONDITIONS_CREATE_START,BOUNDARY_CONDITIONS_DESTROY
  
  PUBLIC BOUNDARY_CONDITIONS_ADD_CONSTANT,BOUNDARY_CONDITIONS_ADD_LOCAL_DOF,BOUNDARY_CONDITIONS_ADD_ELEMENT, &
    & BOUNDARY_CONDITIONS_ADD_NODE,BOUNDARY_CONDITIONS_VARIABLE_GET

  PUBLIC BOUNDARY_CONDITIONS_SET_CONSTANT,BOUNDARY_CONDITIONS_SET_LOCAL_DOF,BOUNDARY_CONDITIONS_SET_ELEMENT, &
    & BOUNDARY_CONDITIONS_SET_NODE, BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE, &
    & BOUNDARY_CONDITIONS_FACE_BASIS_PRESSURE_POISSON_CALCULATE!,BOUNDARY_CONDITIONS_LINE_BASIS_PRESSURE_POISSON_CALCULATE

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
    INTEGER(INTG) :: MPI_IERROR,SEND_COUNT,STORAGE_TYPE, NUMBER_OF_NON_ZEROS, NUMBER_OF_ROWS,COUNT
    INTEGER(INTG) :: variable_idx,dof_idx, equ_matrix_idx, dirichlet_idx, row_idx, DUMMY, LAST, DIRICHLET_DOF
    INTEGER(INTG) :: col_idx,equations_set_idx,bc_idx,parameterSetIdx
    INTEGER(INTG), POINTER :: ROW_INDICES(:), COLUMN_INDICES(:)
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITION_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: VARIABLE_DOMAIN_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: BOUNDARY_CONDITIONS_DIRICHLET
    TYPE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE), POINTER :: BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATION_MATRIX
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES
    TYPE(LIST_TYPE), POINTER :: SPARSE_INDICES
    TYPE(LinkedList),POINTER :: LIST(:)
    INTEGER(INTG), ALLOCATABLE:: COLUMN_ARRAY(:)

    CALL ENTERS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
          IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>0) THEN
            !Transfer all the boundary conditions to all the computational nodes.
 !!TODO \todo Look at this. ?????
            DO variable_idx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
              BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
                FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(VARIABLE_DOMAIN_MAPPING)) THEN
                    SEND_COUNT=VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL
!!This operation is a little expensive as we are doing an unnecessary sum across all the ranks in order to combin
!!the data from each rank into all ranks. We will see how this goes for now.
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%DOF_TYPES, &
                      & SEND_COUNT,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES, &
                      & SEND_COUNT,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="Field variable domain mapping is not associated for variable type "// &
                      & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF

                  ! Update the total number of boundary condition types by summing across all nodes
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS, &
                    & MAX_BOUNDARY_CONDITION_NUMBER,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS, &
                    & 1,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)

                  ! Check that the boundary conditions set are appropriate for equations sets
                  CALL BoundaryConditions_CheckEquations(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)

                  !Make sure the required parameter sets are created on all computational nodes and begin updating them
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%parameterSetRequired, &
                    & FIELD_NUMBER_OF_SET_TYPES,MPI_LOGICAL,MPI_LOR,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
                    IF(BOUNDARY_CONDITION_VARIABLE%parameterSetRequired(parameterSetIdx)) THEN
                      CALL Field_ParameterSetEnsureCreated(FIELD_VARIABLE%FIELD,FIELD_VARIABLE%VARIABLE_TYPE, &
                        & parameterSetIdx,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD_VARIABLE%FIELD,FIELD_VARIABLE%VARIABLE_TYPE, &
                        & parameterSetIdx,ERR,ERROR,*999)
                    END IF
                  END DO

                  ! Set up pressure incremented condition, if it exists
                  IF(BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)>0) THEN
                    CALL BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                    BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED=>BOUNDARY_CONDITION_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED)) THEN
                      ! Find pressure incremented conditions
                      bc_idx=1
                      DO dof_idx=1,FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                        IF(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES(dof_idx) &
                          & ==BOUNDARY_CONDITION_PRESSURE_INCREMENTED) THEN
                          BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED%PRESSURE_INCREMENTED_DOF_INDICES(bc_idx)= &
                            & dof_idx
                          bc_idx=bc_idx+1
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDIF
                  ! Set up Neumann condition information if there are any conditions
                  IF(BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)>0) THEN
                    CALL BOUNDARY_CONDITIONS_NEUMANN_INITIALISE(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                    IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS)) THEN
                      ! Find pressure incremented conditions
                      bc_idx=1
                      DO dof_idx=1,FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                        IF(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES(dof_idx)==BOUNDARY_CONDITION_NEUMANN_POINT) THEN
                          BOUNDARY_CONDITION_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS%SET_DOF(bc_idx)=dof_idx
                          bc_idx=bc_idx+1
                        END IF
                      END DO
                    ELSE
                      CALL FLAG_ERROR("Neumann conditions not associated.",ERR,ERROR,*999)
                    END IF
                  END IF
                  ! Check that there is at least one dirichlet condition
                  IF(BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS>0) THEN
                    CALL BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                    BOUNDARY_CONDITIONS_DIRICHLET=>BOUNDARY_CONDITION_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_DIRICHLET)) THEN
                      ! Find dirichlet conditions
                      dirichlet_idx=1
                      DO dof_idx=1,FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                        IF(BOUNDARY_CONDITION_VARIABLE%DOF_TYPES(dof_idx)==BOUNDARY_CONDITION_DOF_FIXED) THEN
                          BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)=dof_idx
                          dirichlet_idx=dirichlet_idx+1
                        ENDIF
                      ENDDO

                      ! Store Dirichlet dof indices
                      SOLVER_EQUATIONS=>BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
                      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                        IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
                          DO equations_set_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                            EQUATIONS_SET=>SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
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
                                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATION_MATRIX%MATRIX,STORAGE_TYPE,ERR,ERROR,*999)
                                      IF(ASSOCIATED(EQUATION_MATRIX)) THEN
                                        SELECT CASE(STORAGE_TYPE)
                                        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                          !Do nothing
                                        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                          !Do nothing
                                        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for column major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for row major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                          ! Get Sparsity pattern, number of non zeros, number of rows
                                          CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATION_MATRIX%MATRIX,ROW_INDICES, &
                                            & COLUMN_INDICES,ERR,ERROR,*999)
                                          CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET(EQUATION_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
                                            & ERR,ERROR,*999)
                                          ! Get the matrix stored as a linked list
                                          CALL DISTRIBUTED_MATRIX_LINKLIST_GET(EQUATION_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                                          NUMBER_OF_ROWS=EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS
                                          ! Initialise sparsity indices arrays
                                          CALL BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE(BOUNDARY_CONDITIONS_DIRICHLET% &
                                            & LINEAR_SPARSITY_INDICES(equations_set_idx,equ_matrix_idx)%PTR, &
                                            & BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS,ERR,ERROR,*999)

                                          ! Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                                          NULLIFY(SPARSITY_INDICES)
                                          SPARSITY_INDICES=>BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES( &
                                              & equations_set_idx,equ_matrix_idx)%PTR
                                          IF(ASSOCIATED(SPARSITY_INDICES)) THEN
                                            ! Setup list for storing dirichlet non zero indices
                                            NULLIFY(SPARSE_INDICES)
                                            CALL LIST_CREATE_START(SPARSE_INDICES,ERR,ERROR,*999)
                                            CALL LIST_DATA_TYPE_SET(SPARSE_INDICES,LIST_INTG_TYPE,ERR,ERROR,*999)
                                            CALL LIST_INITIAL_SIZE_SET(SPARSE_INDICES, &
                                              & BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS*( &
                                              & NUMBER_OF_NON_ZEROS/NUMBER_OF_ROWS),ERR,ERROR,*999)
                                            CALL LIST_CREATE_FINISH(SPARSE_INDICES,ERR,ERROR,*999)

                                            COUNT=0
                                            SPARSITY_INDICES%SPARSE_COLUMN_INDICES(1)=1
                                            LAST=1
                                            DO dirichlet_idx=1,BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
                                              DIRICHLET_DOF=BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)
                                              CALL LinkedList_to_Array(list(DIRICHLET_DOF),column_array)
                                                DO row_idx=1,size(column_array)
                                                  CALL LIST_ITEM_ADD(SPARSE_INDICES,column_array(row_idx),ERR,ERROR,*999)
                                                  COUNT=COUNT+1
                                                  LAST=row_idx+1
                                                ENDDO
                                              SPARSITY_INDICES%SPARSE_COLUMN_INDICES(dirichlet_idx+1)=COUNT+1
                                            ENDDO
                                            CALL LIST_DETACH_AND_DESTROY(SPARSE_INDICES,DUMMY,SPARSITY_INDICES%SPARSE_ROW_INDICES, &
                                              & ERR,ERROR,*999)
                                            DO col_idx =1,NUMBER_OF_ROWS
                                              CALL LINKEDLIST_DESTROY(list(col_idx))
                                            ENDDO
                                          ELSE
                                            LOCAL_ERROR="Sparsity indices arrays are not associated for this equations matrix."
                                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                          ENDIF
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for compressed column storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for row column storage.",ERR,ERROR,*999)
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
                                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATION_MATRIX%MATRIX,STORAGE_TYPE,ERR,ERROR,*999)
                                      IF(ASSOCIATED(EQUATION_MATRIX)) THEN
                                        SELECT CASE(STORAGE_TYPE)
                                        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                          !Do nothing
                                        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                          !Do nothing
                                        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for column major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for row major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                          ! Get Sparsity pattern, number of non zeros, number of rows
                                          CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATION_MATRIX%MATRIX,ROW_INDICES, &
                                            & COLUMN_INDICES,ERR,ERROR,*999)
                                          CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET(EQUATION_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
                                            & ERR,ERROR,*999)
                                          ! Sparse matrix in a list
                                          CALL DISTRIBUTED_MATRIX_LINKLIST_GET(EQUATION_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                                          NUMBER_OF_ROWS=EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS
                                          ! Intialise sparsity indices arrays
                                          CALL BOUNDARY_CONDITIONS_SPARSITY_INDICES_INITIALISE(BOUNDARY_CONDITIONS_DIRICHLET% &
                                            & DYNAMIC_SPARSITY_INDICES(equations_set_idx,equ_matrix_idx)%PTR, &
                                            & BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS,ERR,ERROR,*999)

                                          ! Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                                          NULLIFY(SPARSITY_INDICES)
                                          SPARSITY_INDICES=>BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES( &
                                              & equations_set_idx,equ_matrix_idx)%PTR
                                          IF(ASSOCIATED(SPARSITY_INDICES)) THEN
                                            ! Setup list for storing dirichlet non zero indices
                                            NULLIFY(SPARSE_INDICES)
                                            CALL LIST_CREATE_START(SPARSE_INDICES,ERR,ERROR,*999)
                                            CALL LIST_DATA_TYPE_SET(SPARSE_INDICES,LIST_INTG_TYPE,ERR,ERROR,*999)
                                            CALL LIST_INITIAL_SIZE_SET(SPARSE_INDICES, &
                                              & BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS*( &
                                              & NUMBER_OF_NON_ZEROS/NUMBER_OF_ROWS),ERR,ERROR,*999)
                                            CALL LIST_CREATE_FINISH(SPARSE_INDICES,ERR,ERROR,*999)

                                            COUNT=0
                                            SPARSITY_INDICES%SPARSE_COLUMN_INDICES(1)=1
                                            LAST=1
                                            DO dirichlet_idx=1,BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
                                              ! Dirichlet columns
                                              DIRICHLET_DOF=BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)
                                              CALL LinkedList_to_Array(list(DIRICHLET_DOF),column_array)
                                              ! The row indices
                                              DO row_idx=1,size(column_array)
                                                CALL LIST_ITEM_ADD(SPARSE_INDICES,column_array(row_idx),ERR,ERROR,*999)
                                                 COUNT=COUNT+1
                                                 LAST=row_idx+1
                                              ENDDO
                                              SPARSITY_INDICES%SPARSE_COLUMN_INDICES(dirichlet_idx+1)=COUNT+1
                                            ENDDO
                                            CALL LIST_DETACH_AND_DESTROY(SPARSE_INDICES,DUMMY,SPARSITY_INDICES%SPARSE_ROW_INDICES, &
                                              & ERR,ERROR,*999)
                                            DO col_idx =1,NUMBER_OF_ROWS
                                              CALL LINKEDLIST_DESTROY(list(col_idx))
                                            ENDDO
                                          ELSE
                                            LOCAL_ERROR="Sparsity indices arrays are not associated for this equations matrix."
                                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                          ENDIF
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for compressed column storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                          CALL FLAG_ERROR("Not implemented for row column storage.",ERR,ERROR,*999)
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
                                  LOCAL_ERROR="Equations Matrices is not associated for these Equations."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                LOCAL_ERROR="Equations is not associated for this Equations Set."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              LOCAL_ERROR="Equations Set is not associated for boundary conditions variable "// &
                                & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !equations_set_idx
                        ELSE
                          LOCAL_ERROR="Solver equations solver mapping is not associated."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Solver equations is not associated."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Dirichlet Boundary Conditions type is not associated for boundary condition variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                  ! Finish field update
                  DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
                    IF(BOUNDARY_CONDITION_VARIABLE%parameterSetRequired(parameterSetIdx)) THEN
                      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD_VARIABLE%FIELD,FIELD_VARIABLE%VARIABLE_TYPE, &
                        & parameterSetIdx,ERR,ERROR,*999)
                    END IF
                  END DO
                ELSE
                  LOCAL_ERROR="Field variable is not associated for variable index "// &
                    & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Boundary conditions variable is not associated for variable index "// &
                    & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//".",ERR,ERROR,*999)
              ENDIF
            ENDDO ! variable_idx
          ENDIF
          !Set the finished flag
          BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.TRUE.
        ELSE
          CALL FLAG_ERROR("Boundary conditions variables array is not allocated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Boundary conditions:",ERR,ERROR,*999)
      DO variable_idx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
        BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Variable type = ",BOUNDARY_CONDITION_VARIABLE%VARIABLE_TYPE, &
            & ERR,ERROR,*999)
        IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
          FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
          VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global dofs = ",VARIABLE_DOMAIN_MAPPING% &
            & NUMBER_OF_GLOBAL,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,8,8, &
            & BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES,'("    Global BCs:",8(X,I8))','(15X,8(X,I8))', &
            & ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Boundary condition variable is not associated",ERR,ERROR,*999)
        ENDIF
      ENDDO !variable_idx
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
  SUBROUTINE BOUNDARY_CONDITIONS_CREATE_START(SOLVER_EQUATIONS,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to create boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<On exit, a pointer to the created boundary conditions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS)) THEN
        CALL FLAG_ERROR("Boundary conditions are already associated for the solver equations.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
          CALL FLAG_ERROR("Boundary conditions is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
            !Initialise the boundary conditions
            CALL BOUNDARY_CONDITIONS_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Solver equations solver mapping is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
          !Return the pointer
          BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
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
    INTEGER(INTG) :: variable_idx

    CALL ENTERS("BOUNDARY_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
        DO variable_idx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
          IF(ASSOCIATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR, &
                & ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Boundary conditions variable number "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                  & " is not associated",ERR,ERROR,*999)
          ENDIF
        ENDDO !variable_idx
        DEALLOCATE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)
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
  SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to initialise the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_idx,VARIABLE_TYPE,equations_set_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS)) THEN
        CALL FLAG_ERROR("Boundary conditions is already associated for these solver equations.",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
          ALLOCATE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate boundary conditions.",ERR,ERROR,*999)
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.FALSE.
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES=0
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%SOLVER_EQUATIONS=>SOLVER_EQUATIONS
          DO equations_set_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                  EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                  IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                    IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
                      EQUATIONS_SET%BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
                      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                        SELECT CASE(EQUATIONS%LINEARITY)
                        CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                            DO variable_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                              VARIABLE_TYPE=LINEAR_MAPPING%LINEAR_MATRIX_VARIABLE_TYPES(variable_idx)
                              IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(VARIABLE_TYPE)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                                CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                  & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(VARIABLE_TYPE)%VARIABLE,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !variable_idx
                          ELSE
                            CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ENDIF
                        CASE(EQUATIONS_NONLINEAR)
                          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                            DO variable_idx=1,NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
                              CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                  & NONLINEAR_MAPPING%RESIDUAL_VARIABLES(variable_idx)%PTR,ERR,ERROR,*999)
                            ENDDO
                          ELSE
                            CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*", &
                                & ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                        SELECT CASE(EQUATIONS%LINEARITY)
                        CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & DYNAMIC_MAPPING%DYNAMIC_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE(EQUATIONS_NONLINEAR)
                          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & DYNAMIC_MAPPING%DYNAMIC_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*", &
                                & ERR,ERROR))//" is invalid."
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
            ELSE
              CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
            ENDIF
          ENDDO !equations_set_idx
        ELSE
          CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_FINALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_INITIALISE")
    RETURN 1

  END SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant. \see OPENCMISS::CMISSBoundaryConditionAddConstant
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_CONSTANT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_CONSTANT",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(DEPENDENT_FIELD_VARIABLE)

    !Note: This routine is for constant interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,local_ny,global_ny, &
            & ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,DEPENDENT_FIELD_VARIABLE,ERR,ERROR,*999)
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
            & ERR,ERROR,*999)
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
              & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated.",ERR,ERROR,*999)
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
  SUBROUTINE BOUNDARY_CONDITIONS_SET_CONSTANT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_CONSTANT",ERR,ERROR,*999)

    !Note: This routine is for constant interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,local_ny,global_ny, &
            & ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
              & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated.",ERR,ERROR,*999)
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
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,(/DOF_INDEX/),(/CONDITION/),(/VALUE/), &
        & ERR,ERROR,*999)

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
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDICES(:) !<DOF_INDICES(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITIONS(:) !<CONDITIONS(:). The boundary condition type to set for the i'th dof \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(:). The value of the boundary condition for the i'th dof to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,global_ny,local_ny
    REAL(DP) :: INITIAL_VALUE
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS",ERR,ERROR,*999)
    NULLIFY(dependent_variable)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          NULLIFY(DEPENDENT_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,DEPENDENT_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
            DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
            IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                & ERR,ERROR,*999)
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                IF(SIZE(DOF_INDICES,1)==SIZE(CONDITIONS,1)) THEN
                  IF(SIZE(DOF_INDICES,1)==SIZE(VALUES,1)) THEN
                    DO i=1,SIZE(DOF_INDICES,1)
                      local_ny=DOF_INDICES(i)
                      IF(local_ny>=1.AND.local_ny<=DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                        global_ny=DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                        ! Set boundary condition and dof type, and make sure parameter sets are created
                        CALL BoundaryConditions_SetConditionType(BOUNDARY_CONDITIONS_VARIABLE,global_ny,CONDITIONS(i), &
                          & ERR,ERROR,*999)
                        ! Update field sets by adding boundary condition values
                        SELECT CASE(CONDITIONS(i))
                        CASE(BOUNDARY_CONDITION_FREE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_FIXED)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INLET)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FREE_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUES(i), &
                            & ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
                          ! For increment loops, we need to set the full BC parameter set value by
                          ! getting the current value from the values parameter set
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,INITIAL_VALUE,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,INITIAL_VALUE+VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                          ! For pressure incremented, adding to the values_set parameter value doesn't make sense,
                          ! so just increment the value in the pressure values parameter set
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
                          ! Point value is stored in boundary conditions field set, and is then integrated to
                          ! get the RHS variable value
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED)
                          ! For integrated Neumann condition, integration is already done, so set the RHS
                          ! dof value directly
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
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
            CALL FLAG_ERROR("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated..",ERR,ERROR,*999)
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
  SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,(/DOF_INDEX/),(/CONDITION/),(/VALUE/), &
      & ERR,ERROR,*999)

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
  SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          NULLIFY(DEPENDENT_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,DEPENDENT_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
            DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
            IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                  & ERR,ERROR,*999)
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                IF(SIZE(DOF_INDICES,1)==SIZE(CONDITIONS,1)) THEN
                  IF(SIZE(DOF_INDICES,1)==SIZE(VALUES,1)) THEN
                    DO i=1,SIZE(DOF_INDICES,1)
                      local_ny=DOF_INDICES(i)
                      IF(local_ny>=1.AND.local_ny<=DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                        global_ny=DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                        ! Set boundary condition and dof type
                        CALL BoundaryConditions_SetConditionType(BOUNDARY_CONDITIONS_VARIABLE,global_ny,CONDITIONS(i), &
                          & ERR,ERROR,*999)
                        ! Update field sets with boundary condition value
                        SELECT CASE(CONDITIONS(i))
                        CASE(BOUNDARY_CONDITION_FREE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_FIXED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INLET)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FREE_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUES(i), &
                            & ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
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
            CALL FLAG_ERROR("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated.",ERR,ERROR,*999)
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

  !> Checks the boundary condition type and sets the boundary condition type and dof type for the boundary conditions.
  !> Makes sure any field parameter sets required are created, and sets the parameter set required array value.
  SUBROUTINE BoundaryConditions_SetConditionType(boundaryConditionsVariable,dofIndex,condition,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: dofIndex !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dofType, previousCondition, previousDof

    CALL ENTERS("BoundaryConditions_SetConditionType",err,error,*999)

    ! We won't do much checking here as this is only used internally and everything has been checked for
    ! association already
    ! Don't need to make sure field values set type is available as this will always be there, but need
    ! to make sure any other parameter sets required are.
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_FREE)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_FIXED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_INLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_MOVED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FREE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_PRESSURE)
      ! Pressure boundary conditions leave the RHS dof as free, as the Neumann terms
      ! are calculated in finite elasticity routines when calculating the element residual
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PRESSURE_VALUES_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PRESSURE_VALUES_SET_TYPE)=.TRUE.
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_PREVIOUS_PRESSURE_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PREVIOUS_PRESSURE_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE DEFAULT
      CALL FLAG_ERROR("The specified boundary condition type for dof number "// &
        & TRIM(NUMBER_TO_VSTRING(dofIndex,"*",err,error))//" of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT

    !We have a valid boundary condition type
    !Update condition type counts
    previousCondition=boundaryConditionsVariable%CONDITION_TYPES(dofIndex)
    IF(previousCondition/=condition) THEN
      ! DOF_COUNTS array doesn't include a count for BOUNDARY_CONDITION_FREE, which equals zero
      IF(previousCondition/=BOUNDARY_CONDITION_FREE) THEN
        boundaryConditionsVariable%DOF_COUNTS(previousCondition)= &
          & boundaryConditionsVariable%DOF_COUNTS(previousCondition)-1
      END IF
      IF(condition/=BOUNDARY_CONDITION_FREE) THEN
        boundaryConditionsVariable%DOF_COUNTS(condition)= &
          & boundaryConditionsVariable%DOF_COUNTS(condition)+1
      END IF
    END IF
    !Update Dirichlet DOF count
    previousDof=boundaryConditionsVariable%DOF_TYPES(dofIndex)
    IF(dofType==BOUNDARY_CONDITION_DOF_FIXED.AND.previousDof==BOUNDARY_CONDITION_DOF_FREE) THEN
      boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS= &
        & boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS+1
    ELSE IF(dofType==BOUNDARY_CONDITION_DOF_FREE.AND.previousDof==BOUNDARY_CONDITION_DOF_FIXED) THEN
      boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS= &
        & boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS-1
    ELSE IF(dofType==BOUNDARY_CONDITION_DOF_MIXED.OR.previousDof==BOUNDARY_CONDITION_DOF_MIXED) THEN
      CALL FLAG_ERROR("Mixed condition DOFs are not implemented.",err,error,*999)
    END IF

    !Set the boundary condition type and DOF type
    boundaryConditionsVariable%CONDITION_TYPES(dofIndex)=condition
    boundaryConditionsVariable%DOF_TYPES(dofIndex)=dofType

    CALL EXITS("BoundaryConditions_SetConditionType")
    RETURN
999 CALL ERRORS("BoundaryConditions_SetConditionType",err,error)
    CALL EXITS("BoundaryConditions_SetConditionType")
    RETURN 1
  END SUBROUTINE BoundaryConditions_SetConditionType

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsAddElement
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_ELEMENT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_ELEMENT",ERR,ERROR,*999)

    !Note: this routine is for element based interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_ELEMENT(FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
            & local_ny,global_ny,ERR,ERROR,*999)
          NULLIFY(FIELD_VARIABLE)
          NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated.",ERR,ERROR,*999)
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
 
  !> Checks that the specified boundary condition is appropriate for the field variable interpolation type
  SUBROUTINE BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*)

    ! Argument variables
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type being set
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field the boundary condition is set on
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type the boundary condition is set on
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number the boundary condition is set on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: interpolationType
    LOGICAL :: validCondition

    CALL ENTERS("BoundaryConditions_CheckInterpolationType",err,error,*999)

    CALL FIELD_COMPONENT_INTERPOLATION_GET(field,variableType,componentNumber,interpolationType,err,error,*999)

    validCondition=.TRUE.
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_FREE, &
        & BOUNDARY_CONDITION_FIXED, &
        & BOUNDARY_CONDITION_FIXED_INCREMENTED)
      ! Valid for all interpolation types
    CASE(BOUNDARY_CONDITION_FIXED_INLET, &
        & BOUNDARY_CONDITION_FIXED_OUTLET)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_FIXED_WALL, &
        & BOUNDARY_CONDITION_MOVED_WALL, &
        & BOUNDARY_CONDITION_FREE_WALL, &
        & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_PRESSURE, &
        & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT, &
        & BOUNDARY_CONDITION_NEUMANN_INTEGRATED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE DEFAULT
      CALL FLAG_ERROR("The specified boundary condition type of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT
    IF(.NOT.validCondition) THEN
      CALL FLAG_ERROR("The specified boundary condition type of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is not valid for the field component "// &
        & "interpolation type of "//TRIM(NUMBER_TO_VSTRING(interpolationType,"*",err,error))//".", &
        & err,error,*999)
    END IF

    CALL EXITS("BoundaryConditions_CheckInterpolationType")
    RETURN
999 CALL ERRORS("BoundaryConditions_CheckInterpolationType",err,error)
    CALL EXITS("BoundaryConditions_CheckInterpolationType")
    RETURN 1
  END SUBROUTINE BoundaryConditions_CheckInterpolationType

  !
  !================================================================================================================================
  !

  !> Checks that the applied boundary conditions are supported by the equations sets in the solver equations
  SUBROUTINE BoundaryConditions_CheckEquations(boundaryConditionsVariable,err,error,*)

    ! Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    type(varying_string), intent(out) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: boundaryConditionType,equationsSetIdx
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    LOGICAL :: validEquationsSetFound

    CALL ENTERS("BoundaryConditions_CheckEquations",err,error,*999)

    !Get and check pointers we need
    solverEquations=>boundaryConditionsVariable%BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
    IF(.NOT.ASSOCIATED(solverEquations)) THEN
      CALL FLAG_ERROR("Boundary conditions solver equations are not associated.",err,error,*999)
    END IF
    solverMapping=>solverEquations%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(solverMapping)) THEN
      CALL FLAG_ERROR("Solver equations solver mapping is not associated.",err,error,*999)
    END IF

    DO boundaryConditionType=1,MAX_BOUNDARY_CONDITION_NUMBER
      !Check if any DOFs have been set to this BC type
      IF(boundaryConditionsVariable%DOF_COUNTS(boundaryConditionType)>0) THEN
        validEquationsSetFound=.FALSE.
        DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
          equationsSet=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%PTR
          IF(.NOT.ASSOCIATED(equationsSet)) THEN
            CALL FLAG_ERROR("Solver equations equations set is not associated.",err,error,*999)
          END IF

          SELECT CASE(boundaryConditionType)
          CASE(BOUNDARY_CONDITION_FREE)
            ! Valid for any equations set
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_INLET, &
              & BOUNDARY_CONDITION_FIXED_OUTLET)
            IF(equationsSet%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                & (equationsSet%TYPE==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                & equationsSet%TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL, &
              & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED,BOUNDARY_CONDITION_FREE_WALL)
            IF(equationsSet%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                & (equationsSet%TYPE==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                & equationsSet%TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE.OR. &
                & equationsSet%TYPE==EQUATIONS_SET_DARCY_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            ELSE IF(equationsSet%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                & equationsSet%TYPE==EQUATIONS_SET_LAPLACE_EQUATION_TYPE.AND. &
                & equationsSet%SUBTYPE==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_PRESSURE, &
              & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
            IF(equationsSet%CLASS==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                & equationsSet%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
            !Not actually used anywhere? So keep it as invalid, although maybe it should be removed?
            validEquationsSetFound=.FALSE.
          CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
            IF(equationsSet%CLASS==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                & equationsSet%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
                & (equationsSet%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE.OR. &
                & equationsSet%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
                & equationsSet%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED)
            validEquationsSetFound=.TRUE.
          CASE DEFAULT
            CALL FLAG_ERROR("The specified boundary condition type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
              & " is invalid.",err,error,*999)
          END SELECT
        END DO

        IF(.NOT.validEquationsSetFound) THEN
            CALL FLAG_ERROR("The specified boundary condition type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
              & " is invalid for the equations sets in the solver equations.",err,error,*999)
        END IF
      END IF
    END DO

    CALL EXITS("BoundaryConditions_CheckEquations")
    RETURN
999 CALL ERRORS("BoundaryConditions_CheckEquations",err,error)
    CALL EXITS("BoundaryConditions_CheckEquations")
    RETURN 1
  END SUBROUTINE BoundaryConditions_CheckEquations

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsSetElement
  SUBROUTINE BOUNDARY_CONDITIONS_SET_ELEMENT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("BOUNDARY_CONDITIONS_SET_ELEMENT",ERR,ERROR,*999)

    !Note: this routine is for element based interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_ELEMENT(FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
            & local_ny,global_ny,ERR,ERROR,*999)
          NULLIFY(FIELD_VARIABLE)
          NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated.",ERR,ERROR,*999)
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
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_NODE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
    & USER_NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The derivative version to set the boundary condition at
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_ADD_NODE",ERR,ERROR,*999)

    NULLIFY(FIELD_VARIABLE)
    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_NODE(FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
            & USER_NODE_NUMBER,COMPONENT_NUMBER,local_ny,global_ny,ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The dependent field variable is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_ADD_NODE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_ADD_NODE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_ADD_NODE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_NODE

  !
  !================================================================================================================================
  !

  !>Adds a specified dofs to the specified boundary condition. 
  SUBROUTINE BOUNDARY_CONDITIONS_NEUMANN_INITIALISE(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE !<A pointer to the boundary conditions variable to initialise Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    INTEGER(INTG) :: NUMBER_OF_VALUES

    CALL ENTERS("BOUNDARY_CONDITIONS_NEUMANN_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      NUMBER_OF_VALUES=BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)
      ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann Boundary Conditions",ERR,ERROR,*999)
      BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
        ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(NUMBER_OF_VALUES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Boundary Conditions dof array.",ERR,ERROR,*999)
        BOUNDARY_CONDITIONS_NEUMANN%SET_DOF=0
      ELSE
        CALL FLAG_ERROR("The boundary condition Neumann is not associated",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Boundary conditions variable is not associated.",ERR,ERROR,*999)
    END IF

    CALL EXITS("BOUNDARY_CONDITIONS_NEUMANN_INITIALISE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_NEUMANN_INITIALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_NEUMANN_INITIALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_NEUMANN_INITIALISE

  !
  !================================================================================================================================
  !

  !>Given listing of faces on boundary calculate value at individual nodes and set boundary condition 
  SUBROUTINE BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,component_idx,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The field to calculate the boundary conditions for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to calculate the boundary conditions for
    INTEGER(INTG), INTENT(IN) :: component_idx !<The component number, i.e. one of u,v,w
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<Pointer to Topology
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    INTEGER(INTG), ALLOCATABLE :: SET_NODES(:),FACES_IN_CALC(:),LINES_IN_CALC(:),DOFS_IN_CALC(:) !Lists of faces, line, nodes and dofs on domain of interest
    INTEGER(INTG), ALLOCATABLE :: FACES_IN_CALC_CORRECTED(:),LINES_IN_CALC_CORRECTED(:) !Lists of faces and lines on domain of interest corrected
    INTEGER(INTG) :: NUMBER_OF_SET_NODES,NUMBER_OF_FACES,NUMBER_OF_LINES,NUMBER_OF_FACES_CORRECTED,NUMBER_OF_LINES_CORRECTED
    INTEGER(INTG) :: NUMBER_OF_SET_DOF,local_ny,nd,nn,nf,nl,nf1,nl1,M,N,i,j,ms,ns,version
    INTEGER(INTG) :: derivative,NUMBER_DOFS_IN_FACE,NUMBER_DOFS_IN_LINE,TOTAL_NUMBER_OF_FACE_DOF
    INTEGER(INTG) :: TOTAL_NUMBER_OF_LINE_DOF,global_ny,MAX_NUMBER_DOFS_IN_FACE,MAX_NUMBER_DOFS_IN_LINE
    INTEGER(INTG) :: FACE_NUMBER,LINE_NUMBER,NODE_NUMBER,x_pos,NODE,face_local_node_index,line_local_node_index
    INTEGER(INTG) :: face_local_deriv_index,line_local_deriv_index,NUMBER_OF_SET_DOF_POINT,NUMBER_OF_SET_DOF_INTEGRATED
    INTEGER(INTG) :: NUMBER_OF_SET_DOF_FREE
    INTEGER(INTG), POINTER :: SET_DOF_POINT(:),SET_DOF_INTEGRATED(:),SET_DOF_FREE(:)
    REAL(DP), POINTER :: SET_DOF_POINT_VALUES(:),SET_DOF_INTEGRATED_VALUES(:)!,SET_DOF_POINT_VALUES_PREV(:)
    LOGICAL :: DOF_LOCATED,INCLUDE_FACE,INCLUDE_LINE
    TYPE(LIST_TYPE), POINTER :: SET_NODES_LIST,FACES_LIST,LINES_LIST,DOFS_LIST
    TYPE(LIST_TYPE), POINTER :: FACES_LIST_CORRECTED,LINES_LIST_CORRECTED
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_INTEGRATED_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ASSOCIATED(FIELD)) THEN
        FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
              TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY
              IF(ASSOCIATED(TOPOLOGY)) THEN
                DOMAIN_NODES=>TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  SELECT CASE(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    IF(FIELD%DECOMPOSITION% &
                      & DOMAIN(FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
                      & PTR%NUMBER_OF_DIMENSIONS==3) THEN

                      DOMAIN_FACES=>TOPOLOGY%FACES
                      IF(ASSOCIATED(DOMAIN_FACES)) THEN

                        NUMBER_OF_SET_NODES=0 !Counts the number of nodes that correspond to the set dofs
                        NUMBER_OF_FACES=0 !Counts the number of unique faces of interest
                        NUMBER_OF_FACES_CORRECTED=0 !The number of faces with all dof specified either FREE or POINT
                        TOTAL_NUMBER_OF_FACE_DOF=0 !Counts the number of unique dof in faces of interest

                        !Calculate number of rows in SET_DOF
                        NUMBER_OF_SET_DOF=SIZE(BOUNDARY_CONDITIONS_NEUMANN%SET_DOF,1)
                        !This includes all set dof for all components, derivatives, and BOUNDARY_CONDITION_NEUMANN_POINT
                        !and BOUNDARY_CONDITION_NEUMANN_INTEGRATED conditions

                        !SET_DOF_POINT and SET_DOF_POINT_VALUES are to contain a subset of all
                        !BOUNDARY_CONDITION_NEUMANN_POINT set DOFs in SET_DOF and SET_DOF_VALUES
                        ALLOCATE(SET_DOF_POINT(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof array.",ERR,ERROR,*999)
                        ALLOCATE(SET_DOF_POINT_VALUES(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof values array."&
                          & ,ERR,ERROR,*999)
                        !ALLOCATE(SET_DOF_POINT_VALUES_PREV(NUMBER_OF_SET_DOF),STAT=ERR)
                        !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate array set dof point values from previous time step."&
                        !  & ,ERR,ERROR,*999)

                        !SET_DOF_INTEGRATED and SET_DOF_INTEGRATED_VALUES are to contain a subset of all
                        !BOUNDARY_CONDITION_NEUMANN_INTEGRATED set DOFs in SET_DOF and SET_DOF_VALUES
                        ALLOCATE(SET_DOF_INTEGRATED(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof array.",ERR,ERROR,*999)
                        ALLOCATE(SET_DOF_INTEGRATED_VALUES(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof values array."&
                          & ,ERR,ERROR,*999)

                        !SET_DOF_FREE is to contain the dof not explicitly set, but included in the calculation
                        ALLOCATE(SET_DOF_FREE(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof array of free dof."&
                          & ,ERR,ERROR,*999)

                        SET_DOF_POINT=0
                        SET_DOF_POINT_VALUES=0.0_DP
                        !SET_DOF_POINT_VALUES_PREV=0.0_DP
                        SET_DOF_INTEGRATED=0
                        SET_DOF_INTEGRATED_VALUES=0.0_DP
                        SET_DOF_FREE=0


                        !Create listing of nodes on which DOF have been set
                        NULLIFY(SET_NODES_LIST)
                        CALL LIST_CREATE_START(SET_NODES_LIST,ERR,ERROR,*999)
                        CALL LIST_DATA_TYPE_SET(SET_NODES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                        CALL LIST_INITIAL_SIZE_SET(SET_NODES_LIST,NUMBER_OF_SET_DOF,ERR,ERROR,*999)
                        CALL LIST_CREATE_FINISH(SET_NODES_LIST,ERR,ERROR,*999)

                        NUMBER_OF_SET_DOF_POINT=0
                        NUMBER_OF_SET_DOF_INTEGRATED=0
                        NUMBER_OF_SET_DOF_FREE=0

                        !Separate out the POINT, INTEGRATED and FREE set dof
                        DO nd=1,NUMBER_OF_SET_DOF
                          !Check if dof in current component (uvwp)
                          IF(FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                            & NODE_DOF2PARAM_MAP(4,BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd))==component_idx) THEN
                            global_ny=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)
                            local_ny=FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_NUMBER(1)

                            IF(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(global_ny)==BOUNDARY_CONDITION_NEUMANN_POINT) THEN

                              NUMBER_OF_SET_DOF_POINT=NUMBER_OF_SET_DOF_POINT+1
                              SET_DOF_POINT(NUMBER_OF_SET_DOF_POINT)=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)
                              ! Get BC set value for point
                              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                                & local_ny,SET_DOF_POINT_VALUES(NUMBER_OF_SET_DOF_POINT),ERR,ERROR,*999)

                              !Convert dof to node and add to list
                              CALL LIST_ITEM_ADD(SET_NODES_LIST,FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                                & NODE_DOF2PARAM_MAP(3,BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)),ERR,ERROR,*999)

                            ELSE IF(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(global_ny)== &
                              & BOUNDARY_CONDITION_NEUMANN_INTEGRATED) THEN

                              NUMBER_OF_SET_DOF_INTEGRATED=NUMBER_OF_SET_DOF_INTEGRATED+1
                              SET_DOF_INTEGRATED(NUMBER_OF_SET_DOF_INTEGRATED)=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)
                              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                                & local_ny,SET_DOF_INTEGRATED_VALUES(NUMBER_OF_SET_DOF_INTEGRATED),ERR,ERROR,*999)

                            ELSE IF(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(global_ny)==BOUNDARY_CONDITION_FREE) THEN

                              NUMBER_OF_SET_DOF_FREE=NUMBER_OF_SET_DOF_FREE+1
                              SET_DOF_FREE(NUMBER_OF_SET_DOF_FREE)=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)

                            ENDIF
                          ENDIF
                        ENDDO !nd

                        CALL LIST_REMOVE_DUPLICATES(SET_NODES_LIST,ERR,ERROR,*999)
                        CALL LIST_DETACH_AND_DESTROY(SET_NODES_LIST,NUMBER_OF_SET_NODES,SET_NODES,ERR,ERROR,*999)

                        !If all Set DOF values set to Integrated, skip calculation and write values directly to RHS
                        IF(NUMBER_OF_SET_DOF_INTEGRATED>0.AND.NUMBER_OF_SET_DOF_POINT==0 &
                          & .AND.NUMBER_OF_SET_DOF_FREE==0) THEN
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN% &
                            & INTEGRATED_VALUES_VECTOR(NUMBER_OF_SET_DOF_INTEGRATED),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Integrated Values Vector."&
                            & ,ERR,ERROR,*999)
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN% &
                            & INTEGRATED_VALUES_VECTOR_MAPPING(NUMBER_OF_SET_DOF_INTEGRATED),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Integrated Values Vector Mapping."&
                            & ,ERR,ERROR,*999)
                          BOUNDARY_CONDITIONS_NEUMANN% &
                               & INTEGRATED_VALUES_VECTOR(1:NUMBER_OF_SET_DOF_INTEGRATED) = &
                               & SET_DOF_INTEGRATED_VALUES(1:NUMBER_OF_SET_DOF_INTEGRATED)
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(1:NUMBER_OF_SET_DOF_INTEGRATED)= &
                               & SET_DOF_INTEGRATED(1:NUMBER_OF_SET_DOF_INTEGRATED)
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_SIZE= &
                               & NUMBER_OF_SET_DOF_INTEGRATED
                          DEALLOCATE(SET_DOF_POINT)
                          DEALLOCATE(SET_DOF_POINT_VALUES)
                          !DEALLOCATE(SET_DOF_POINT_VALUES_PREV)
                          DEALLOCATE(SET_DOF_INTEGRATED)
                          DEALLOCATE(SET_DOF_INTEGRATED_VALUES)
                          DEALLOCATE(SET_DOF_FREE)
                          DEALLOCATE(SET_NODES)
                        ELSE

                          !Create list of surrounding faces
                          NULLIFY(FACES_LIST)
                          CALL LIST_CREATE_START(FACES_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(FACES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(FACES_LIST,DOMAIN_FACES%NUMBER_OF_FACES,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(FACES_LIST,ERR,ERROR,*999)

                          MAX_NUMBER_DOFS_IN_FACE=0
                          !Find the faces this node is a part of, by looping over set nodes that correspond to set dof
                          DO nn=1,NUMBER_OF_SET_NODES
                            DOMAIN_NODE=>DOMAIN_NODES%NODES(SET_NODES(nn))
                            DO face_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_FACES
                              DO nf1=1,DOMAIN_FACES%NUMBER_OF_FACES
                                DOMAIN_FACE=>DOMAIN_FACES%FACES(nf1)
                                IF(DOMAIN_FACE%NUMBER==DOMAIN_NODE%NODE_FACES(face_local_node_index)) THEN
                                  IF(DOMAIN_FACE%BOUNDARY_FACE) THEN
                                    CALL LIST_ITEM_ADD(FACES_LIST,DOMAIN_NODE% &
                                      & NODE_FACES(face_local_node_index),ERR,ERROR,*999)
                                    MAX_NUMBER_DOFS_IN_FACE=MAX(MAX_NUMBER_DOFS_IN_FACE,DOMAIN_FACE% &
                                      & BASIS%NUMBER_OF_ELEMENT_PARAMETERS)
                                  ENDIF
                                ENDIF
                              ENDDO !nf1
                            ENDDO !face_local_node_index
                          ENDDO !nn

                          CALL LIST_REMOVE_DUPLICATES(FACES_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(FACES_LIST,NUMBER_OF_FACES,FACES_IN_CALC,ERR,ERROR,*999)


                          !Create sub-set of faces for the calculation, by removing faces where all dof have not been set
                          NULLIFY(FACES_LIST_CORRECTED)
                          CALL LIST_CREATE_START(FACES_LIST_CORRECTED,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(FACES_LIST_CORRECTED,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(FACES_LIST_CORRECTED,&
                            & DOMAIN_FACES%NUMBER_OF_FACES,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(FACES_LIST_CORRECTED,ERR,ERROR,*999)

                          !Iterate over all faces in domain of interest
                          DO nf=1,NUMBER_OF_FACES
                            INCLUDE_FACE=.TRUE.
                            !Match face number to domain face number
                            DO nf1=1,DOMAIN_FACES%NUMBER_OF_FACES
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(nf1)
                              IF(DOMAIN_FACE%NUMBER==FACES_IN_CALC(nf).AND.DOMAIN_FACE%BOUNDARY_FACE) THEN
                                !Need whether dof have all been specified for face, if not disclude from calculation

                                !Iterate over number of nodes in face
                                DO face_local_node_index=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES !nnf
                                  NODE=DOMAIN_FACE%NODES_IN_FACE(face_local_node_index)
                                  !Iterate over number of derivatives in face
                                  DO face_local_deriv_index=1,DOMAIN_FACE%BASIS% &
                                                              & NUMBER_OF_DERIVATIVES(face_local_node_index) !nnd
                                    IF(INCLUDE_FACE) THEN
                                      DERIVATIVE=DOMAIN_FACE%DERIVATIVES_IN_FACE(1,face_local_deriv_index, &
                                        & face_local_node_index)
                                      version=DOMAIN_FACE%DERIVATIVES_IN_FACE(2,face_local_deriv_index, &
                                        & face_local_node_index)
                                      !Locate dof number on face
                                      local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP%NODES(node)%DERIVATIVES(derivative)% &
                                        & VERSIONS(version)
                                      !Check if that dof is specified as either point or free
                                      DOF_LOCATED=.FALSE.
                                      DO i=1,NUMBER_OF_SET_DOF_POINT
                                        IF(SET_DOF_POINT(i)==local_ny) THEN
                                          DOF_LOCATED=.TRUE.
                                        ENDIF
                                      ENDDO !i

                                      DO i=1,NUMBER_OF_SET_DOF_FREE
                                        IF(SET_DOF_FREE(i)==local_ny) THEN
                                          DOF_LOCATED=.TRUE.
                                        ENDIF
                                      ENDDO !i

                                      !If failed to locate dof then exclude face from list
                                      IF(DOF_LOCATED.EQV..FALSE.) THEN
                                        INCLUDE_FACE=.FALSE.
                                      ENDIF

                                    ENDIF
                                  ENDDO !face_local_deriv_index
                                ENDDO !face_local_node_index
                              ENDIF
                            ENDDO !nf1

                            IF(INCLUDE_FACE) THEN
                              CALL LIST_ITEM_ADD(FACES_LIST_CORRECTED,FACES_IN_CALC(nf),ERR,ERROR,*999)
                            ENDIF
                          ENDDO !nf

                          CALL LIST_REMOVE_DUPLICATES(FACES_LIST_CORRECTED,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(FACES_LIST_CORRECTED,NUMBER_OF_FACES_CORRECTED,&
                            & FACES_IN_CALC_CORRECTED,ERR,ERROR,*999)



                          !Use a list to collect all the dof, sort and substitute into INTEGRATION_MATRIX_MAPPING array
                          NULLIFY(DOFS_LIST)
                          CALL LIST_CREATE_START(DOFS_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(DOFS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(DOFS_LIST,NUMBER_OF_FACES_CORRECTED* &
                                                                            & MAX_NUMBER_DOFS_IN_FACE,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(DOFS_LIST,ERR,ERROR,*999)

                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF &
                                                      & (DOMAIN_FACES%NUMBER_OF_FACES,MAX_NUMBER_DOFS_IN_FACE),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Faces element param to &
                            & local DOF lookup array.",ERR,ERROR,*999)

                          BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF = 0

                          !Find the faces this node is a part of, by looping over set nodes that correspond to set dof (global dofs)
                          DO nf=1,NUMBER_OF_FACES_CORRECTED
                            DO nf1=1,DOMAIN_FACES%NUMBER_OF_FACES
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(nf1)
                              IF(DOMAIN_FACE%NUMBER==FACES_IN_CALC_CORRECTED(nf).AND.DOMAIN_FACE%BOUNDARY_FACE) THEN

                                DO face_local_node_index=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES !nnf

                                  NODE=DOMAIN_FACE%NODES_IN_FACE(face_local_node_index)
                                  DO face_local_deriv_index=1,DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(face_local_node_index) !nnd
                                    DERIVATIVE=DOMAIN_FACE%DERIVATIVES_IN_FACE(1,face_local_deriv_index,face_local_node_index)
                                    version=DOMAIN_FACE%DERIVATIVES_IN_FACE(2,face_local_deriv_index,face_local_node_index)
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                      & NODES(node)%DERIVATIVES(derivative)%VERSIONS(version)
                                    x_pos = DOMAIN_FACE%BASIS% &
                                            & ELEMENT_PARAMETER_INDEX(face_local_deriv_index,face_local_node_index)

                                    CALL LIST_ITEM_ADD(DOFS_LIST,local_ny,ERR,ERROR,*999)

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

                          CALL LIST_REMOVE_DUPLICATES(DOFS_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(DOFS_LIST,TOTAL_NUMBER_OF_FACE_DOF,&
                                                                                  & DOFS_IN_CALC,ERR,ERROR,*999)


                          !Number of DOF involved in boundary condition calculation
                          M = TOTAL_NUMBER_OF_FACE_DOF
                          !Number of DOF where boundary conditions have been fixed by user
                          N = NUMBER_OF_SET_DOF_POINT

                          !Integration_matrix = 'A' matrix
                          !Point_values_vector = 'x' vector
                          !Integrated_values_vector = 'B' vector

                          !Set size of matrix and vectors
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN% &
                            & FACE_INTEGRATION_MATRIX(MAX_NUMBER_DOFS_IN_FACE),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann face integration matrix array."&
                            & ,ERR,ERROR,*999)
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING(MAX_NUMBER_DOFS_IN_FACE &
                            & ),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann face &
                            & integration matrix mapping array.",ERR,ERROR,*999)
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
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X = &
                                                                            & SET_DOF_POINT(1:NUMBER_OF_SET_DOF_POINT)
                          BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING = &
                                                                            & SET_DOF_POINT(1:NUMBER_OF_SET_DOF_POINT)
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y = &
                                                                     & DOFS_IN_CALC(1:TOTAL_NUMBER_OF_FACE_DOF)
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING = &
                                                                     & DOFS_IN_CALC(1:TOTAL_NUMBER_OF_FACE_DOF)
                          BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR = &
                                                                     & SET_DOF_POINT_VALUES(1:NUMBER_OF_SET_DOF_POINT)

                          !Iterate over set dofs
                          DO nd=1,NUMBER_OF_SET_DOF_POINT

                            !Convert dof to node number
                            NODE_NUMBER=FIELD_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,SET_DOF_POINT(nd))
                            DOMAIN_NODE=>DOMAIN_NODES%NODES(NODE_NUMBER)

                            !Iterate over the faces on which that node is a part of
                            DO face_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_FACES

                              INCLUDE_FACE=.FALSE.
                              !Check face is included amongst faces involved in calculation
                              DO nf=1,NUMBER_OF_FACES_CORRECTED
                                IF(FACES_IN_CALC_CORRECTED(nf)==DOMAIN_NODE%NODE_FACES(face_local_node_index)) THEN
                                  INCLUDE_FACE=.TRUE.
                                ENDIF
                              ENDDO

                              IF(INCLUDE_FACE) THEN
                                DO nf1=1,DOMAIN_FACES%NUMBER_OF_FACES
                                  DOMAIN_FACE=>DOMAIN_FACES%FACES(nf1)
                                  ns = 0

                                  !Locate face and ensure face is a boundary face
                                  IF(DOMAIN_FACE%NUMBER==DOMAIN_NODE%NODE_FACES(face_local_node_index) &
                                    & .AND.(DOMAIN_FACE%BOUNDARY_FACE)) THEN

                                    FACE_NUMBER = DOMAIN_NODE%NODE_FACES(face_local_node_index)

                                    !Calculate the ns corresponding to the dof
                                    DO ms=1,DOMAIN_FACE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                                      local_ny=BOUNDARY_CONDITIONS_NEUMANN% &
                                                                   & FACES_ELEMENT_PARAM_2_LOCAL_DOF(FACE_NUMBER,ms)
                                      IF(local_ny==SET_DOF_POINT(nd)) THEN
                                        ns = ms
                                      ENDIF
                                    ENDDO !ms

                                    !Calculate basis for face
                                    CALL BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
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
                              ENDIF
                            ENDDO !face_local_node_index
                          ENDDO !nd

                          !Perform Ax=B calculation
                          DO j=1,TOTAL_NUMBER_OF_FACE_DOF
                            DO i=1,NUMBER_OF_SET_DOF_POINT
                              BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) = &
                                & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) + &
                                & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,i) * &
                                & BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR(i)
                            ENDDO !i
!!Check here that global_ny is the best thing to do here
                            global_ny=BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)
                          ENDDO !j

                          !For selected DOF, override calculated values with user defined BOUNDARY_CONDITION_NEUMANN_INTEGRATED values
                          IF(NUMBER_OF_SET_DOF_INTEGRATED>0) THEN
                            DO j=1,TOTAL_NUMBER_OF_FACE_DOF
                              DO nd=1,NUMBER_OF_SET_DOF_INTEGRATED
                                IF(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)== &
                                                                             & SET_DOF_INTEGRATED(nd)) THEN
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j)= &
                                                                             & SET_DOF_INTEGRATED_VALUES(nd)

                                  global_ny=BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)
                                ENDIF
                              ENDDO !nd
                            ENDDO !j
                          ENDIF

                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_SIZE=TOTAL_NUMBER_OF_FACE_DOF

                          DEALLOCATE(SET_DOF_POINT)
                          DEALLOCATE(SET_DOF_POINT_VALUES)
                          !DEALLOCATE(SET_DOF_POINT_VALUES_PREV)
                          DEALLOCATE(SET_DOF_INTEGRATED)
                          DEALLOCATE(SET_DOF_INTEGRATED_VALUES)
                          IF(NUMBER_OF_SET_DOF_FREE>0) DEALLOCATE(SET_DOF_FREE)
                          DEALLOCATE(SET_NODES)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING)
                          DEALLOCATE(FACES_IN_CALC_CORRECTED)
                          DEALLOCATE(FACES_IN_CALC)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Topology faces is not associated.",ERR,ERROR,*999)
                      ENDIF


                    ELSEIF(FIELD%DECOMPOSITION% &
                      & DOMAIN(FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
                      & PTR%NUMBER_OF_DIMENSIONS==2) THEN

                      DOMAIN_LINES=>TOPOLOGY%LINES
                      IF(ASSOCIATED(DOMAIN_LINES)) THEN

                        NUMBER_OF_SET_NODES=0 !Counts the number of nodes that correspond to the set dofs
                        NUMBER_OF_LINES=0 !Counts the number of unique lines of interest
                        NUMBER_OF_LINES_CORRECTED=0 !The number of lines with all dof specified either FREE or POINT
                        TOTAL_NUMBER_OF_LINE_DOF=0 !Counts the number of unique dof in lines of interest

                        !Calculate number of rows in SET_DOF
                        NUMBER_OF_SET_DOF=SIZE(BOUNDARY_CONDITIONS_NEUMANN%SET_DOF,1)
                        !This includes all set dof for all components, derivatives, and BOUNDARY_CONDITION_NEUMANN_POINT
                        !and BOUNDARY_CONDITION_NEUMANN_INTEGRATED conditions

                        !SET_DOF_POINT and SET_DOF_POINT_VALUES are to contain a subset of all
                        !BOUNDARY_CONDITION_NEUMANN_POINT set DOFs in SET_DOF and SET_DOF_VALUES
                        ALLOCATE(SET_DOF_POINT(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof array.",ERR,ERROR,*999)
                        ALLOCATE(SET_DOF_POINT_VALUES(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof values array."&
                          & ,ERR,ERROR,*999)
                        !ALLOCATE(SET_DOF_POINT_VALUES_PREV(NUMBER_OF_SET_DOF),STAT=ERR)
                        !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate array set dof point values from previous time step."&
                        !  & ,ERR,ERROR,*999)

                        !SET_DOF_INTEGRATED and SET_DOF_INTEGRATED_VALUES are to contain a subset of all
                        !BOUNDARY_CONDITION_NEUMANN_INTEGRATED set DOFs in SET_DOF and SET_DOF_VALUES
                        ALLOCATE(SET_DOF_INTEGRATED(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof array.",ERR,ERROR,*999)
                        ALLOCATE(SET_DOF_INTEGRATED_VALUES(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof values array."&
                          & ,ERR,ERROR,*999)

                        !SET_DOF_FREE is to contain the dof that are not set, but are to be included in the calculation
                        ALLOCATE(SET_DOF_FREE(NUMBER_OF_SET_DOF),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary set dof array of free dof."&
                          & ,ERR,ERROR,*999)

                        SET_DOF_POINT=0
                        SET_DOF_POINT_VALUES=0.0_DP
                        !SET_DOF_POINT_VALUES_PREV=0.0_DP
                        SET_DOF_INTEGRATED=0
                        SET_DOF_INTEGRATED_VALUES=0.0_DP
                        SET_DOF_FREE=0


                        !Create listing of nodes on which DOF have been set
                        NULLIFY(SET_NODES_LIST)
                        CALL LIST_CREATE_START(SET_NODES_LIST,ERR,ERROR,*999)
                        CALL LIST_DATA_TYPE_SET(SET_NODES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                        CALL LIST_INITIAL_SIZE_SET(SET_NODES_LIST,NUMBER_OF_SET_DOF,ERR,ERROR,*999)
                        CALL LIST_CREATE_FINISH(SET_NODES_LIST,ERR,ERROR,*999)

                        NUMBER_OF_SET_DOF_POINT=0
                        NUMBER_OF_SET_DOF_INTEGRATED=0
                        NUMBER_OF_SET_DOF_FREE=0

                        DO nd=1,NUMBER_OF_SET_DOF
                          !Check if dof in current component (uvwp)
                          IF(FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                            & NODE_DOF2PARAM_MAP(4,BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd))==component_idx) THEN
                            global_ny=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)
                            local_ny=FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_NUMBER(1)

                            IF(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(global_ny)==BOUNDARY_CONDITION_NEUMANN_POINT) THEN

                              NUMBER_OF_SET_DOF_POINT=NUMBER_OF_SET_DOF_POINT+1
                              SET_DOF_POINT(NUMBER_OF_SET_DOF_POINT)=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)
                              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                                & local_ny,SET_DOF_POINT_VALUES(NUMBER_OF_SET_DOF_POINT),ERR,ERROR,*999)

                              !Convert dof to node and add to list
                              CALL LIST_ITEM_ADD(SET_NODES_LIST,FIELD_VARIABLE%DOF_TO_PARAM_MAP% &
                                & NODE_DOF2PARAM_MAP(3,BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)),ERR,ERROR,*999)

                            ELSE IF(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(global_ny)== &
                              & BOUNDARY_CONDITION_NEUMANN_INTEGRATED) THEN

                              NUMBER_OF_SET_DOF_INTEGRATED=NUMBER_OF_SET_DOF_INTEGRATED+1
                              SET_DOF_INTEGRATED(NUMBER_OF_SET_DOF_INTEGRATED)=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)
                              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                                & local_ny,SET_DOF_INTEGRATED_VALUES(NUMBER_OF_SET_DOF_INTEGRATED),ERR,ERROR,*999)

                            ELSE IF(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(global_ny)== &
                              & BOUNDARY_CONDITION_FREE) THEN

                              NUMBER_OF_SET_DOF_FREE=NUMBER_OF_SET_DOF_FREE+1
                              SET_DOF_FREE(NUMBER_OF_SET_DOF_FREE)=BOUNDARY_CONDITIONS_NEUMANN%SET_DOF(nd)

                            ENDIF
                          ENDIF
                        ENDDO !nd

                        CALL LIST_REMOVE_DUPLICATES(SET_NODES_LIST,ERR,ERROR,*999)
                        CALL LIST_DETACH_AND_DESTROY(SET_NODES_LIST, &
                                        & NUMBER_OF_SET_NODES,SET_NODES,ERR,ERROR,*999)

                        !If all Set DOF values set to Integrated, skip calculation and write values directly to RHS
                        IF(NUMBER_OF_SET_DOF_INTEGRATED>0.AND.NUMBER_OF_SET_DOF_POINT==0 &
                          & .AND.NUMBER_OF_SET_DOF_FREE==0) THEN
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN% &
                            & INTEGRATED_VALUES_VECTOR(NUMBER_OF_SET_DOF_INTEGRATED),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Integrated Values Vector."&
                            & ,ERR,ERROR,*999)
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN% &
                            & INTEGRATED_VALUES_VECTOR_MAPPING(NUMBER_OF_SET_DOF_INTEGRATED),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Integrated Values Vector Mapping."&
                            & ,ERR,ERROR,*999)
                          BOUNDARY_CONDITIONS_NEUMANN% &
                               & INTEGRATED_VALUES_VECTOR(1:NUMBER_OF_SET_DOF_INTEGRATED)= &
                               & SET_DOF_INTEGRATED_VALUES(1:NUMBER_OF_SET_DOF_INTEGRATED)
                          BOUNDARY_CONDITIONS_NEUMANN% &
                               & INTEGRATED_VALUES_VECTOR_MAPPING(1:NUMBER_OF_SET_DOF_INTEGRATED)= &
                               & SET_DOF_INTEGRATED(1:NUMBER_OF_SET_DOF_INTEGRATED)
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_SIZE= &
                               & NUMBER_OF_SET_DOF_INTEGRATED
                          DEALLOCATE(SET_DOF_POINT)
                          DEALLOCATE(SET_DOF_POINT_VALUES)
                          !DEALLOCATE(SET_DOF_POINT_VALUES_PREV)
                          DEALLOCATE(SET_DOF_INTEGRATED)
                          DEALLOCATE(SET_DOF_INTEGRATED_VALUES)
                          DEALLOCATE(SET_DOF_FREE)
                          DEALLOCATE(SET_NODES)
                        ELSE

                          !Create list of surrounding lines
                          NULLIFY(LINES_LIST)
                          CALL LIST_CREATE_START(LINES_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(LINES_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(LINES_LIST,DOMAIN_LINES%NUMBER_OF_LINES,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(LINES_LIST,ERR,ERROR,*999)

                          MAX_NUMBER_DOFS_IN_LINE=0
                          !Find the lines this node is a part of, by looping over set nodes that correspond to set dof
                          DO nn=1,NUMBER_OF_SET_NODES
                            DOMAIN_NODE=>DOMAIN_NODES%NODES(SET_NODES(nn))
                            DO line_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_LINES

                              DO nl1=1,DOMAIN_LINES%NUMBER_OF_LINES
                                DOMAIN_LINE=>DOMAIN_LINES%LINES(nl1)
                                IF(DOMAIN_LINE%NUMBER==DOMAIN_NODE%NODE_LINES(line_local_node_index)) THEN
                                  IF(DOMAIN_LINE%BOUNDARY_LINE) THEN

                                    CALL LIST_ITEM_ADD(LINES_LIST,DOMAIN_NODE% &
                                                              & NODE_LINES(line_local_node_index),ERR,ERROR,*999)
                                    MAX_NUMBER_DOFS_IN_LINE=MAX(MAX_NUMBER_DOFS_IN_LINE,DOMAIN_LINE% &
                                                                          & BASIS%NUMBER_OF_ELEMENT_PARAMETERS)
                                  ENDIF
                                ENDIF
                              ENDDO !nl1
                            ENDDO !line_local_node_index
                          ENDDO !nn

                          CALL LIST_REMOVE_DUPLICATES(LINES_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(LINES_LIST,NUMBER_OF_LINES,LINES_IN_CALC,ERR,ERROR,*999)


                          !Create sub-set of lines for the calculation, by removing lines where all dof have not been set
                          NULLIFY(LINES_LIST_CORRECTED)
                          CALL LIST_CREATE_START(LINES_LIST_CORRECTED,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(LINES_LIST_CORRECTED,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(LINES_LIST_CORRECTED,&
                            & DOMAIN_LINES%NUMBER_OF_LINES,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(LINES_LIST_CORRECTED,ERR,ERROR,*999)

                          !Iterate over all lines in domain of interest
                          DO nl=1,NUMBER_OF_LINES
                            INCLUDE_LINE=.TRUE.
                            !Match line number to domain line number
                            DO nl1=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(nl1)
                              IF(DOMAIN_LINE%NUMBER==LINES_IN_CALC(nl).AND.DOMAIN_LINE%BOUNDARY_LINE) THEN
                                !Need whether dof have all been specified for line, if not disclude from calculation

                                !Iterate over number of nodes in line
                                DO line_local_node_index=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES !nnl
                                  NODE=DOMAIN_LINE%NODES_IN_LINE(line_local_node_index)
                                  !Iterate over number of derivatives in line
                                  DO line_local_deriv_index=1,DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(line_local_node_index) !nnd
                                    IF(INCLUDE_LINE) THEN
                                      DERIVATIVE=DOMAIN_LINE%DERIVATIVES_IN_LINE(1,line_local_deriv_index,line_local_node_index)
                                      version=DOMAIN_LINE%DERIVATIVES_IN_LINE(2,line_local_deriv_index,line_local_node_index)
                                      !Locate dof number on face
                                      local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                        & NODES(node)%DERIVATIVES(derivative)%VERSIONS(version)

                                      !Check if that dof is specified as either point or free
                                      DOF_LOCATED=.FALSE.
                                      DO i=1,NUMBER_OF_SET_DOF_POINT
                                        IF(SET_DOF_POINT(i)==local_ny) THEN
                                          DOF_LOCATED=.TRUE.
                                        ENDIF
                                      ENDDO !i

                                      DO i=1,NUMBER_OF_SET_DOF_FREE
                                        IF(SET_DOF_FREE(i)==local_ny) THEN
                                          DOF_LOCATED=.TRUE.
                                        ENDIF
                                      ENDDO !i

                                      !If failed to locate dof then exclude line from list
                                      IF(DOF_LOCATED.EQV..FALSE.) THEN
                                        INCLUDE_LINE=.FALSE.
                                      ENDIF

                                    ENDIF
                                  ENDDO !line_local_deriv_index
                                ENDDO !line_local_node_index
                              ENDIF
                            ENDDO !nl1

                            IF(INCLUDE_LINE) THEN
                              CALL LIST_ITEM_ADD(LINES_LIST_CORRECTED,LINES_IN_CALC(nl),ERR,ERROR,*999)
                            ENDIF
                          ENDDO !nl

                          CALL LIST_REMOVE_DUPLICATES(LINES_LIST_CORRECTED,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(LINES_LIST_CORRECTED,NUMBER_OF_LINES_CORRECTED,&
                            & LINES_IN_CALC_CORRECTED,ERR,ERROR,*999)



                          !Use a list to collect all the dof, sort and substitute into INTEGRATION_MATRIX_MAPPING array
                          NULLIFY(DOFS_LIST)
                          CALL LIST_CREATE_START(DOFS_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(DOFS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(DOFS_LIST,NUMBER_OF_LINES_CORRECTED* &
                                                                              & MAX_NUMBER_DOFS_IN_LINE,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(DOFS_LIST,ERR,ERROR,*999)

                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF &
                                                      & (DOMAIN_LINES%NUMBER_OF_LINES,MAX_NUMBER_DOFS_IN_LINE),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Lines element param to &
                            & local DOF lookup array.",ERR,ERROR,*999)

                          BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF = 0

                          !Find the lines this node is a part of, by looping over set nodes that correspond to set dof (global dofs)
                          DO nl=1,NUMBER_OF_LINES_CORRECTED
                            DO nl1=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(nl1)
                              IF(DOMAIN_LINE%NUMBER==LINES_IN_CALC_CORRECTED(nl).AND.DOMAIN_LINE%BOUNDARY_LINE) THEN

                                DO line_local_node_index=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES !nnl

                                  NODE=DOMAIN_LINE%NODES_IN_LINE(line_local_node_index)
                                  DO line_local_deriv_index=1,DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(line_local_node_index) !nnd
                                    DERIVATIVE=DOMAIN_LINE%DERIVATIVES_IN_LINE(1,line_local_deriv_index,line_local_node_index)
                                    version=DOMAIN_LINE%DERIVATIVES_IN_LINE(2,line_local_deriv_index,line_local_node_index)
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                      & NODES(node)%DERIVATIVES(derivative)%VERSIONS(version)
                                    x_pos = DOMAIN_LINE%BASIS% &
                                      & ELEMENT_PARAMETER_INDEX(line_local_deriv_index,line_local_node_index)

                                    CALL LIST_ITEM_ADD(DOFS_LIST,local_ny,ERR,ERROR,*999)

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

                          CALL LIST_REMOVE_DUPLICATES(DOFS_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(DOFS_LIST, &
                                                           & TOTAL_NUMBER_OF_LINE_DOF,DOFS_IN_CALC,ERR,ERROR,*999)

                          !Number of DOF involved in boundary condition calculation
                          M = TOTAL_NUMBER_OF_LINE_DOF
                          !Number of DOF where boundary conditions have been fixed by user
                          N = NUMBER_OF_SET_DOF_POINT

                          !Integration_matrix = 'A' matrix
                          !Point_values_vector = 'x' vector
                          !Integrated_values_vector = 'B' vector

                          !Set size of matrix and vectors
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN% &
                            & LINE_INTEGRATION_MATRIX(MAX_NUMBER_DOFS_IN_LINE),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann line integration matrix array."&
                            & ,ERR,ERROR,*999)
                          ALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX_MAPPING(MAX_NUMBER_DOFS_IN_LINE &
                            & ),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Neumann line &
                            & integration matrix mapping array.",ERR,ERROR,*999)
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
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X =  &
                                                                            & SET_DOF_POINT(1:NUMBER_OF_SET_DOF_POINT)
                          BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING = &
                                                                            & SET_DOF_POINT(1:NUMBER_OF_SET_DOF_POINT)
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y = &
                                                                     & DOFS_IN_CALC(1:TOTAL_NUMBER_OF_LINE_DOF)
                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING= &
                                                                     & DOFS_IN_CALC(1:TOTAL_NUMBER_OF_LINE_DOF)
                          BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR = &
                                                                     & SET_DOF_POINT_VALUES(1:NUMBER_OF_SET_DOF_POINT)

                          !Iterate over set dofs
                          DO nd=1,NUMBER_OF_SET_DOF_POINT

                            !Convert dof to node number
                            NODE_NUMBER=FIELD_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,SET_DOF_POINT(nd))
                            DOMAIN_NODE=>DOMAIN_NODES%NODES(NODE_NUMBER)

                            !Iterate over the lines on which that node is a part of
                            DO line_local_node_index=1,DOMAIN_NODE%NUMBER_OF_NODE_LINES

                              INCLUDE_LINE=.FALSE.
                              !Check line is included amongst lines involved in calculation
                              DO nl=1,NUMBER_OF_LINES_CORRECTED
                                IF(LINES_IN_CALC_CORRECTED(nl)==DOMAIN_NODE%NODE_LINES(line_local_node_index)) THEN
                                  INCLUDE_LINE=.TRUE.
                                ENDIF
                              ENDDO

                              IF(INCLUDE_LINE) THEN
                                DO nl1=1,DOMAIN_LINES%NUMBER_OF_LINES
                                  DOMAIN_LINE=>DOMAIN_LINES%LINES(nl1)
                                  ns = 0

                                  !Locate line and ensure line is a boundary line
                                  IF(DOMAIN_LINE%NUMBER==DOMAIN_NODE%NODE_LINES(line_local_node_index) &
                                    & .AND.(DOMAIN_LINE%BOUNDARY_LINE)) THEN

                                    LINE_NUMBER = DOMAIN_NODE%NODE_LINES(line_local_node_index)

                                    !Calculate the ns corresponding to the dof
                                    DO ms=1,DOMAIN_LINE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                                      local_ny=BOUNDARY_CONDITIONS_NEUMANN% &
                                                                   & LINES_ELEMENT_PARAM_2_LOCAL_DOF(LINE_NUMBER,ms)
                                      IF(local_ny==SET_DOF_POINT(nd)) THEN
                                        ns = ms
                                      ENDIF
                                    ENDDO !ms

                                    !Calculate basis for line
                                    CALL BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
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
                              ENDIF
                            ENDDO !line_local_node_index
                          ENDDO !nd

                          !Perform Ax=B calculation
                          DO j=1,TOTAL_NUMBER_OF_LINE_DOF
                            DO i=1,NUMBER_OF_SET_DOF_POINT
                              BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) = &
                                & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j) + &
                                & BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX(j,i) * &
                                & BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR(i)
                            ENDDO !i
!!Check here that global_ny is the best thing to do here
                            global_ny=BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)
                          ENDDO !j

                          !For selected DOF, override calculated values with user defined BOUNDARY_CONDITION_NEUMANN_INTEGRATED values
                          IF(NUMBER_OF_SET_DOF_INTEGRATED>0) THEN
                            DO j=1,TOTAL_NUMBER_OF_LINE_DOF
                              DO nd=1,NUMBER_OF_SET_DOF_INTEGRATED
                                IF(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)== &
                                                                               & SET_DOF_INTEGRATED(nd)) THEN
                                  BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR(j)= &
                                                                               & SET_DOF_INTEGRATED_VALUES(nd)

                                  global_ny=BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_MAPPING(j)
                                ENDIF
                              ENDDO !nd
                            ENDDO !j
                          ENDIF

                          BOUNDARY_CONDITIONS_NEUMANN%INTEGRATED_VALUES_VECTOR_SIZE=TOTAL_NUMBER_OF_LINE_DOF

                          DEALLOCATE(SET_DOF_POINT)
                          DEALLOCATE(SET_DOF_POINT_VALUES)
                          !DEALLOCATE(SET_DOF_POINT_VALUES_PREV)
                          DEALLOCATE(SET_DOF_INTEGRATED)
                          DEALLOCATE(SET_DOF_INTEGRATED_VALUES)
                          IF(NUMBER_OF_SET_DOF_FREE>0) DEALLOCATE(SET_DOF_FREE)
                          DEALLOCATE(SET_NODES)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINES_ELEMENT_PARAM_2_LOCAL_DOF)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX_MAPPING)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_X)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%INTEGRATION_MATRIX_MAPPING_Y)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR)
                          DEALLOCATE(BOUNDARY_CONDITIONS_NEUMANN%POINT_VALUES_VECTOR_MAPPING)
                          DEALLOCATE(LINES_IN_CALC_CORRECTED)
                          DEALLOCATE(LINES_IN_CALC)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Topology lines is not associated.",ERR,ERROR,*999)
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
                  CALL FLAG_ERROR("Topology nodes is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The boundary condition Neumann is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type "&
                      &//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " has not been created."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
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
  SUBROUTINE BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,component_idx,FACE_NUMBER, &
    & NUMBER_DOFS_IN_FACE,ns,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: component_idx !<The component number, i.e. one of u,v,w
    INTEGER(INTG), INTENT(IN) :: FACE_NUMBER !<The face number (the OpenCMISS face number)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_DOFS_IN_FACE !<The number of dofs in the face
    INTEGER(INTG), INTENT(IN) :: ns !<The element parameter number of the set dof
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG) :: ELEMENT_NUMBER,ng,ms,local_ny,equations_set_idx
    REAL(DP) :: RWG,PHIMS,PHINS
    LOGICAL :: EQUATIONS_SET_FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_FACE_BASIS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
          BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
            !Just take the first equations set from the solver equations that has the correct dependent field
            !for calculating the interpolation information.
            !If there is more than one equations set in the solver equations and they share a
            !dependent field variable, then the geometric field is also the same
            !\todo: Should interpolated_point etc stuff be associated with field variables rather than equations sets?
            SOLVER_EQUATIONS=>BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              equations_set_idx=0
              EQUATIONS_SET_FOUND = .FALSE.
              DO WHILE(.NOT.EQUATIONS_SET_FOUND .AND. &
                & equations_set_idx<SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                EQUATIONS_SET=>SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  !Check if equations set has the same dependent field we're after
                  IF(ASSOCIATED(FIELD_VARIABLE%FIELD,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)) THEN
                    EQUATIONS_SET_FOUND=.TRUE.
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
                equations_set_idx=equations_set_idx+1
              ENDDO
              IF(.NOT.EQUATIONS_SET_FOUND) THEN
                CALL FLAG_ERROR("Couldn't find an equations set with this dependent field.",ERR,ERROR,*999)
              ENDIF
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
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

                          !Scale factor adjustment required for handling with Cubic Hermite elements
                          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
                              & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                            DO ms=1,DOMAIN_FACE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms)= &
                                & BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms)* &
                                & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                                & PTR%SCALE_FACTORS(ms,component_idx)* &
                                & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                                & PTR%SCALE_FACTORS(ns,component_idx)
                            ENDDO !ms
                          ENDIF

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
                        CALL FLAG_ERROR("Quadrature scheme is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain topology face is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain topology faces is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The boundary condition Neumann is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The boundary conditions for variable type " &
                  & //TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " has not been created."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
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

  !>For a given face, calculates the contributions of the basis functions and integrated flux to RHS for the Pressure Poisson case
  SUBROUTINE BOUNDARY_CONDITIONS_FACE_BASIS_PRESSURE_POISSON_CALCULATE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,component_idx, &
    & FACE_NUMBER,NUMBER_DOFS_IN_FACE,ns,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: component_idx !<The component number, i.e. one of u,v,w
    INTEGER(INTG), INTENT(IN) :: FACE_NUMBER !<The face number (the OpenCMISS face number)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_DOFS_IN_FACE !<The number of dofs in the face
    INTEGER(INTG), INTENT(IN) :: ns !<The element parameter number of the set dof
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG) :: ELEMENT_NUMBER,ng,ms,local_ny,ni,mi,equations_set_idx
    REAL(DP) :: RWG,PHIMS,PHINS,DEL_PHINS,NONLINEAR_SUM,RHO_PARAM,MU_PARAM
    REAL(DP) :: U_VALUE(3),U_DERIV(3,3),DXI_DX(3,3)
    LOGICAL :: EQUATIONS_SET_FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BOUNDARY_CONDITIONS_FACE_BASIS_PRESSURE_POISSON_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
          BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
            SOLVER_EQUATIONS=>BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              equations_set_idx=0
              EQUATIONS_SET_FOUND = .FALSE.
              DO WHILE(.NOT.EQUATIONS_SET_FOUND .AND. &
                & equations_set_idx<SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                EQUATIONS_SET=>SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  !Check if equations set has the same dependent field we're after
                  IF(ASSOCIATED(FIELD_VARIABLE%FIELD,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)) THEN
                    EQUATIONS_SET_FOUND=.TRUE.
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
                equations_set_idx=equations_set_idx+1
              ENDDO
              IF(.NOT.EQUATIONS_SET_FOUND) THEN
                CALL FLAG_ERROR("Couldn't find an equations set with this dependent field.",ERR,ERROR,*999)
              ENDIF
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
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

                          !Define MU_PARAM, viscosity=1
                          MU_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                          !Define RHO_PARAM, density=2
                          RHO_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & VALUES(2,NO_PART_DERIV)

                          !Loop over gauss points
                          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS

                            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

                            !Calculate RWG
                            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN &
                                  & *QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

                            DO ni=1,DOMAIN_FACE%BASIS%NUMBER_OF_XI
                              DO mi=1,DOMAIN_FACE%BASIS%NUMBER_OF_XI
                                DXI_DX(mi,ni)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                                  & DXI_DX(mi,ni)
                              ENDDO !mi
                            ENDDO !ni

                            U_VALUE(1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(1,NO_PART_DERIV)
                            U_VALUE(2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(2,NO_PART_DERIV)
                            U_DERIV(1,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(1,PART_DERIV_S1)
                            U_DERIV(1,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(1,PART_DERIV_S2)
                            U_DERIV(2,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(2,PART_DERIV_S1)
                            U_DERIV(2,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(2,PART_DERIV_S2)
                            IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                              U_VALUE(3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(3,NO_PART_DERIV)
                              U_DERIV(3,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(3,PART_DERIV_S1)
                              U_DERIV(3,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(3,PART_DERIV_S2)
                              U_DERIV(3,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(3,PART_DERIV_S3)
                              U_DERIV(1,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(1,PART_DERIV_S3)
                              U_DERIV(2,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VARIABLE%VARIABLE_TYPE) &
                                         & %PTR%VALUES(2,PART_DERIV_S3)
                            ELSE
                              U_VALUE(3)=0.0_DP
                              U_DERIV(3,1)=0.0_DP
                              U_DERIV(3,2)=0.0_DP
                              U_DERIV(3,3)=0.0_DP
                              U_DERIV(1,3)=0.0_DP
                              U_DERIV(2,3)=0.0_DP
                            ENDIF

                            NONLINEAR_SUM=0.0_DP
                            !Calculate SUM
                            DO ni=1,DOMAIN_FACE%BASIS%NUMBER_OF_XI
                              NONLINEAR_SUM=NONLINEAR_SUM+RHO_PARAM*( &
                                & (U_VALUE(1))*(U_DERIV(component_idx,ni)*DXI_DX(ni,1))+ &
                                & (U_VALUE(2))*(U_DERIV(component_idx,ni)*DXI_DX(ni,2))+ &
                                & (U_VALUE(3))*(U_DERIV(component_idx,ni)*DXI_DX(ni,3)))
                            ENDDO !ni

                            !ns is singular because we are doing an integration of the ns with all ms in the face
                            PHINS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                            DEL_PHINS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,FIRST_PART_DERIV,ng)

                            !Loop over element rows
                            DO ms=1,DOMAIN_FACE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS

                              PHIMS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                              !Calculate mass matrix (M)
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms) = &
                                & BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms) + RHO_PARAM*PHIMS*PHINS*RWG

                              local_ny = BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF(FACE_NUMBER,ms)
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX_MAPPING(ms) = local_ny

                              !Calculate stiffness matrix (K)
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_STIFFNESS_MATRIX(ms) = &
                                & BOUNDARY_CONDITIONS_NEUMANN%FACE_STIFFNESS_MATRIX(ms) + MU_PARAM*PHIMS*DEL_PHINS*RWG

                              local_ny = BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF(FACE_NUMBER,ms)
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_STIFFNESS_MATRIX_MAPPING(ms) = local_ny

                              !Calculate Nonlinear matrix (NL)
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_NONLINEAR_MATRIX(ms) = &
                                & BOUNDARY_CONDITIONS_NEUMANN%FACE_NONLINEAR_MATRIX(ms) + PHIMS*NONLINEAR_SUM*RWG

                              local_ny = BOUNDARY_CONDITIONS_NEUMANN%FACES_ELEMENT_PARAM_2_LOCAL_DOF(FACE_NUMBER,ms)
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_NONLINEAR_MATRIX_MAPPING(ms) = local_ny

                            ENDDO !ms
                          ENDDO !ng

                          !Scale factor adjustment required for handling with Cubic Hermite elements
                          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
                              & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                            DO ms=1,DOMAIN_FACE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                              BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms)= &
                                & BOUNDARY_CONDITIONS_NEUMANN%FACE_INTEGRATION_MATRIX(ms)* &
                                & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                                & PTR%SCALE_FACTORS(ms,component_idx)* &
                                & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                                & PTR%SCALE_FACTORS(ns,component_idx)
                            ENDDO !ms
                          ENDIF

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
                        CALL FLAG_ERROR("Quadrature scheme is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain topology face is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain topology faces is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The boundary condition Neumann is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The boundary conditions for variable type " &
                  & //TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " has not been created."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_FACE_BASIS_PRESSURE_POISSON_CALCULATE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_FACE_BASIS_PRESSURE_POISSON_CALCULATE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_FACE_BASIS_PRESSURE_POISSON_CALCULATE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_FACE_BASIS_PRESSURE_POISSON_CALCULATE

  !
  !================================================================================================================================
  !

  !>For a given line, calculates the contributions of the basis functions and integrated flux to RHS
  SUBROUTINE BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,component_idx,LINE_NUMBER, &
    & NUMBER_DOFS_IN_LINE,ns,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: component_idx !<The component number, i.e. one of u,v,w
    INTEGER(INTG), INTENT(IN) :: LINE_NUMBER !<The line number (the OpenCMISS line number)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_DOFS_IN_LINE !<The number of dofs in the line
    INTEGER(INTG), INTENT(IN) :: ns !<The element parameter number of the set dof
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_NEUMANN_TYPE), POINTER :: BOUNDARY_CONDITIONS_NEUMANN
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG) :: ELEMENT_NUMBER,ng,ms,local_ny,equations_set_idx
    REAL(DP) :: RWG,PHIMS,PHINS
    LOGICAL :: EQUATIONS_SET_FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BOUNDARY_CONDITIONS_LINE_BASIS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
          BOUNDARY_CONDITIONS_NEUMANN=>BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_NEUMANN)) THEN
            SOLVER_EQUATIONS=>BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              equations_set_idx=0
              EQUATIONS_SET_FOUND = .FALSE.
              DO WHILE(.NOT.EQUATIONS_SET_FOUND .AND. &
                & equations_set_idx<SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                EQUATIONS_SET=>SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  !Check if equations set has the same dependent field we're after
                  IF(ASSOCIATED(FIELD_VARIABLE%FIELD,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)) THEN
                    EQUATIONS_SET_FOUND=.TRUE.
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
                equations_set_idx=equations_set_idx+1
              ENDDO
              IF(.NOT.EQUATIONS_SET_FOUND) THEN
                CALL FLAG_ERROR("Couldn't find an equations set with this dependent field.",ERR,ERROR,*999)
              ENDIF
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
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

                            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_LINE_TYPE, &
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

                          !Scale factor adjustment required for handling with Cubic Hermite elements
                          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
                              & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                            DO ms=1,DOMAIN_LINE%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                              BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX(ms)= &
                                & BOUNDARY_CONDITIONS_NEUMANN%LINE_INTEGRATION_MATRIX(ms)* &
                                & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                                & PTR%SCALE_FACTORS(ms,component_idx)* &
                                & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                                & PTR%SCALE_FACTORS(ns,component_idx)
                            ENDDO !ms
                          ENDIF

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
                CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The boundary condition Neumann is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The boundary conditions for variable type " &
                  & //TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " has not been created."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
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
  SUBROUTINE BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
    & USER_NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The derivative version to set the boundary condition at
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BOUNDARY_CONDITIONS_SET_NODE",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(FIELD_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FLAG_ERROR("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_NODE(FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
            & USER_NODE_NUMBER,COMPONENT_NUMBER,local_ny,global_ny,ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The dependent field variable is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The dependent field is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BOUNDARY_CONDITIONS_SET_NODE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_SET_NODE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_SET_NODE")
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
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: BOUNDARY_CONDITIONS_DIRICHLET

    CALL ENTERS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ALLOCATED(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES))  &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES)
      IF(ALLOCATED(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES))  &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS)) THEN
        BOUNDARY_CONDITIONS_DIRICHLET=>BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS
        CALL BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE(BOUNDARY_CONDITIONS_DIRICHLET% &
            & LINEAR_SPARSITY_INDICES,ERR,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE(BOUNDARY_CONDITIONS_DIRICHLET% &
            & DYNAMIC_SPARSITY_INDICES,ERR,ERROR,*999)
        IF(ALLOCATED(BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES)) THEN
          DEALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES)
        ENDIF
        DEALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET)
      ENDIF
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS)) &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)) &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)
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

  !>Finalise an array of sparcity indices and deallocate all memory.
  SUBROUTINE BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE(SPARSITY_INDICES_ARRAY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_PTR_TYPE), ALLOCATABLE :: SPARSITY_INDICES_ARRAY(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equ_set_idx, equ_matrix_idx
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES
    
    CALL ENTERS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE",ERR,ERROR,*999)
    
    IF (ALLOCATED(SPARSITY_INDICES_ARRAY)) THEN
      DO equ_set_idx=1,SIZE(SPARSITY_INDICES_ARRAY,1)
        DO equ_matrix_idx=1,SIZE(SPARSITY_INDICES_ARRAY,2)
          SPARSITY_INDICES=>SPARSITY_INDICES_ARRAY(equ_set_idx,equ_matrix_idx)%PTR
          IF(ASSOCIATED(SPARSITY_INDICES)) THEN
            IF(ALLOCATED(SPARSITY_INDICES%SPARSE_ROW_INDICES)) THEN
              DEALLOCATE(SPARSITY_INDICES%SPARSE_ROW_INDICES)
            ENDIF
            IF(ALLOCATED(SPARSITY_INDICES%SPARSE_COLUMN_INDICES)) THEN
              DEALLOCATE(SPARSITY_INDICES%SPARSE_COLUMN_INDICES)
            ENDIF
            DEALLOCATE(SPARSITY_INDICES)
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(SPARSITY_INDICES_ARRAY)
    ENDIF
    
    CALL EXITS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SPARSITY_INDICES_ARRAY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the boundary conditions variable for a variable type if that variable has not already been initialised, otherwise do nothing.
  SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(BOUNDARY_CONDITIONS,FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to initialise a variable type for.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable to initialise the boundary conditions variable for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_idx
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: VARIABLE_DOMAIN_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_PTR_TYPE), ALLOCATABLE :: NEW_BOUNDARY_CONDITIONS_VARIABLES(:)
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE

    CALL ENTERS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
        IF(ASSOCIATED(VARIABLE_DOMAIN_MAPPING)) THEN
          NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
          !Check if boundary conditions variable has already been added, if so then we don't do anything as different equations
          !sets can have the same dependent field variables and will both want to add the variable
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            ALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES(BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new boundary conditions variables array.",ERR,ERROR,*998)
            IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
              DO variable_idx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
                NEW_BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR=> &
                    & BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
              ENDDO
            ENDIF

            ALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES(BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate boundary condition variable.",ERR,ERROR,*998)
            BOUNDARY_CONDITIONS_VARIABLE=>NEW_BOUNDARY_CONDITIONS_VARIABLES( &
                & BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1)%PTR
            BOUNDARY_CONDITIONS_VARIABLE%BOUNDARY_CONDITIONS=>BOUNDARY_CONDITIONS
            BOUNDARY_CONDITIONS_VARIABLE%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            BOUNDARY_CONDITIONS_VARIABLE%VARIABLE=>FIELD_VARIABLE
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global boundary condition types.",ERR,ERROR,*999)
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES(VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global boundary condition dof types.",ERR,ERROR,*999)
            BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES=BOUNDARY_CONDITION_FREE
            BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES=BOUNDARY_CONDITION_DOF_FREE
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(MAX_BOUNDARY_CONDITION_NUMBER),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate boundary condition DOF counts array.",ERR,ERROR,*999)
            BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS=0
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS)
            BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS=0
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%NEUMANN_BOUNDARY_CONDITIONS)
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%parameterSetRequired(FIELD_NUMBER_OF_SET_TYPES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate boundary condition parameter set required array.",ERR,ERROR,*999)
            BOUNDARY_CONDITIONS_VARIABLE%parameterSetRequired=.FALSE.
            BOUNDARY_CONDITIONS_VARIABLE%parameterSetRequired(FIELD_VALUES_SET_TYPE)=.TRUE.

            CALL MOVE_ALLOC(NEW_BOUNDARY_CONDITIONS_VARIABLES,BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)
            BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES= &
                & BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1
          ENDIF
        ELSE
          CALL FLAG_ERROR("Field variable domain mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS_VARIABLE,DUMMY_ERR,DUMMY_ERROR,*998)
    DEALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES)
998 CALL ERRORS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions variable for a given field variable
  SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to initialise a variable type for.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable to initialise the boundary conditions variable for.
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER, INTENT(OUT) :: BOUNDARY_CONDITIONS_VARIABLE !<On return, a pointer to the boundary conditions variable, or NULL if it wasn't found
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    LOGICAL :: VARIABLE_FOUND

    CALL ENTERS("BOUNDARY_CONDITIONS_VARIABLE_GET",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        VARIABLE_FOUND=.FALSE.
        variable_idx=1
        DO WHILE(variable_idx<=BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES.AND..NOT.VARIABLE_FOUND)
          VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR%VARIABLE
          IF(ASSOCIATED(VARIABLE)) THEN
            IF(VARIABLE%VARIABLE_TYPE==FIELD_VARIABLE%VARIABLE_TYPE.AND. &
                & VARIABLE%FIELD%USER_NUMBER==FIELD_VARIABLE%FIELD%USER_NUMBER.AND. &
                & VARIABLE%FIELD%REGION%USER_NUMBER==FIELD_VARIABLE%FIELD%REGION%USER_NUMBER) THEN
              VARIABLE_FOUND=.TRUE.
              BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
            ENDIF
          ELSE
            CALL FLAG_ERROR("Boundary conditions variable field variable is not associated.",ERR,ERROR,*999)
          ENDIF
          variable_idx=variable_idx+1
        ENDDO
      ELSE
        CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_GET")
    RETURN
999 CALL ERRORS("BOUNDARY_CONDITIONS_VARIABLE_GET",ERR,ERROR)
    CALL EXITS("BOUNDARY_CONDITIONS_VARIABLE_GET")
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_GET

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
    INTEGER(INTG) :: NUMBER_OF_DIRICHLET_CONDITIONS,NUMBER_OF_LINEAR_MATRICES,NUMBER_OF_DYNAMIC_MATRICES,matrix_idx, &
      & MAX_NUMBER_LINEAR_MATRICES,MAX_NUMBER_DYNAMIC_MATRICES,equations_set_idx
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
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

        SOLVER_EQUATIONS=>BOUNDARY_CONDITIONS_VARIABLE%BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          MAX_NUMBER_LINEAR_MATRICES=0
          MAX_NUMBER_DYNAMIC_MATRICES=0
          DO equations_set_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                  LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                  DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                  IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                    NUMBER_OF_LINEAR_MATRICES=LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    IF(NUMBER_OF_LINEAR_MATRICES>MAX_NUMBER_LINEAR_MATRICES) &
                      & MAX_NUMBER_LINEAR_MATRICES=NUMBER_OF_LINEAR_MATRICES
                  ENDIF
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    NUMBER_OF_DYNAMIC_MATRICES=DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                    IF(NUMBER_OF_DYNAMIC_MATRICES>MAX_NUMBER_DYNAMIC_MATRICES) &
                      & MAX_NUMBER_DYNAMIC_MATRICES=NUMBER_OF_DYNAMIC_MATRICES
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
          ENDDO
          ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES(SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS, &
                & MAX_NUMBER_LINEAR_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Dirichlet linear sparsity indices array",ERR,ERROR,*999)
          ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES(SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,&
                & MAX_NUMBER_DYNAMIC_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Dirichlet dynamic sparsity indices array",ERR,ERROR,*999)
          DO equations_set_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            DO matrix_idx=1,MAX_NUMBER_LINEAR_MATRICES
              NULLIFY(BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES(equations_set_idx,matrix_idx)%PTR)
            ENDDO
            DO matrix_idx=1,MAX_NUMBER_DYNAMIC_MATRICES
              NULLIFY(BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES(equations_set_idx,matrix_idx)%PTR)
            ENDDO
          ENDDO
        ELSE
          CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*998)
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
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sparsity indicies.",ERR,ERROR,*999)
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

  !>Initialises the pressure incremented boundary condition.
  SUBROUTINE BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE !<A pointer to the boundary conditions variable to initialise a boundary conditions dirichlet type for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE), POINTER :: BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED
    INTEGER(INTG) :: NUMBER_OF_PRESSURE_INCREMENTED_CONDITIONS

    CALL ENTERS("BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)) THEN
        CALL FLAG_ERROR("Pressure incremented boundary conditions are already associated for this boundary conditions variable." &
           & ,ERR,ERROR,*998)
      ELSE
        ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Pressure incremented Boundary Conditions",ERR,ERROR,*999)
        BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED=>BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS
        NUMBER_OF_PRESSURE_INCREMENTED_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        ALLOCATE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED%PRESSURE_INCREMENTED_DOF_INDICES &
          & (NUMBER_OF_PRESSURE_INCREMENTED_CONDITIONS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Pressure incremented DOF indices array",ERR,ERROR,*999)
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
  END SUBROUTINE BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE BOUNDARY_CONDITIONS_ROUTINES
