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
    & BOUNDARY_CONDITIONS_SET_NODE,BoundaryConditions_NeumannIntegrate

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
    INTEGER(INTG) :: col_idx,equations_set_idx,parameterSetIdx
    INTEGER(INTG) :: pressureIdx,neumannIdx
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

    NULLIFY(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED)

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
                  END IF
                  ! Set up Neumann condition information if there are any Neumann conditions
                  IF(BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)>0) THEN
                    CALL BoundaryConditions_NeumannInitialise(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                  END IF

                  ! Loop over all global DOFs, keeping track of the dof indices of specific BC types where required
                  pressureIdx=1
                  neumannIdx=1
                  DO dof_idx=1,FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                    IF(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES(dof_idx)== BOUNDARY_CONDITION_PRESSURE_INCREMENTED) THEN
                      BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED%PRESSURE_INCREMENTED_DOF_INDICES(pressureIdx)=dof_idx
                      pressureIdx=pressureIdx+1
                    ELSE IF(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES(dof_idx)==BOUNDARY_CONDITION_NEUMANN_POINT) THEN
                      BOUNDARY_CONDITION_VARIABLE%neumannBoundaryConditions%setDofs(neumannIdx)=dof_idx
                      neumannIdx=neumannIdx+1
                    END IF
                  END DO

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

  !>Initialise the Neumann boundary conditions information
  SUBROUTINE BoundaryConditions_NeumannInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann
    INTEGER(INTG) :: numberOfValues,numberOfLocalDofs

    CALL ENTERS("BoundaryConditions_NeumannInitialise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      numberOfValues=boundaryConditionsVariable%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)
      ALLOCATE(boundaryConditionsVariable%neumannBoundaryConditions,stat=err)
      IF(err/=0) CALL FLAG_ERROR("Could not allocate Neumann Boundary Conditions",err,error,*999)
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        numberOfLocalDofs=boundaryConditionsVariable%VARIABLE%NUMBER_OF_DOFS
        ALLOCATE(boundaryConditionsNeumann%setDofs(numberOfValues),stat=err)
        boundaryConditionsNeumann%setDofs=0
        IF(err/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix.",err,error,*999)
        ALLOCATE(boundaryConditionsNeumann%integrationMatrix(numberOfLocalDofs,numberOfValues),stat=err)
        IF(err/=0) CALL FLAG_ERROR("Could not allocate Neumann integration matrix.",err,error,*999)
        ALLOCATE(boundaryConditionsNeumann%pointValues(numberOfValues),stat=err)
        IF(err/=0) CALL FLAG_ERROR("Could not allocate Neumann point values vector.",err,error,*999)
        ALLOCATE(boundaryConditionsNeumann%integratedValues(numberOfLocalDofs),stat=err)
        IF(err/=0) CALL FLAG_ERROR("Could not allocate Neumann integrated values vector.",err,error,*999)
        ALLOCATE(boundaryConditionsNeumann%localDofs(numberOfValues),stat=err)
        IF(err/=0) CALL FLAG_ERROR("Could not allocate Neumann local DOFs vector.",err,error,*999)
      ELSE
        CALL FLAG_ERROR("The boundary condition Neumann is not associated",err,error,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    CALL EXITS("BoundaryConditions_NeumannInitialise")
    RETURN
999 CALL ERRORS("BoundaryConditions_NeumannInitialise",err,error)
    CALL EXITS("BoundaryConditions_NeumannInitialise")
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannInitialise

  !
  !================================================================================================================================
  !

  !Finalise the Neumann condition information for a boundary conditions variable
  SUBROUTINE BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann

    CALL ENTERS("BoundaryConditions_NeumannFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        DEALLOCATE(boundaryConditionsNeumann%setDofs)
        DEALLOCATE(boundaryConditionsNeumann%integrationMatrix)
        DEALLOCATE(boundaryConditionsNeumann%pointValues)
        DEALLOCATE(boundaryConditionsNeumann%integratedValues)
        DEALLOCATE(boundaryConditionsNeumann%localDofs)
        DEALLOCATE(boundaryConditionsNeumann)
      END IF
    ELSE
      CALL FLAG_ERROR("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    CALL EXITS("BoundaryConditions_NeumannFinalise")
    RETURN
999 CALL ERRORS("BoundaryConditions_NeumannFinalise",err,error)
    CALL EXITS("BoundaryConditions_NeumannFinalise")
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannFinalise

  !
  !================================================================================================================================
  !

  !>Calculate the integrated Neumann conditions for a field variable and update the RHS field variable value
  SUBROUTINE BoundaryConditions_NeumannIntegrate(rhsBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER, INTENT(IN) :: rhsBoundaryConditions !<The boundary conditions for the right hand side field variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
    INTEGER(INTG) :: componentIdx,globalDof,localDof,numberOfNeumann
    INTEGER(INTG) :: faceIdx,lineIdx,nodeIdx,derivIdx,gaussIdx
    INTEGER(INTG) :: faceNumberOfNeumann,lineNumberOfNeumann
    INTEGER(INTG) :: neumannConditionNumber(16) !Maps from face/line Neumann condition index to global Neumann condition index
    INTEGER(INTG) :: neumannNodeAndDeriv(16,2) !Store node and derivative number of Neumann condition DOF in face/line
    INTEGER(INTG) :: neumannConditionIdx,elementNeumannDofIdx
    INTEGER(INTG) :: matrixNeumannDofIdx
    INTEGER(INTG) :: neumannNode,neumannDeriv
    INTEGER(INTG) :: ms,os,localNode,version
    INTEGER(INTG) :: mpiError
    REAL(DP) :: integratedValue,pointValue,phim,phio
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannConditions
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointMetrics(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:), scalingParameters(:)
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology
    TYPE(DOMAIN_FACES_TYPE), POINTER :: faces
    TYPE(DOMAIN_LINES_TYPE), POINTER :: lines
    TYPE(DOMAIN_FACE_TYPE), POINTER :: face
    TYPE(DOMAIN_LINE_TYPE), POINTER :: line
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme

    CALL ENTERS("BoundaryConditions_NeumannIntegrate",err,error,*999)

    NULLIFY(scalingParameters)
    NULLIFY(interpolationParameters)
    NULLIFY(interpolatedPoints)
    NULLIFY(interpolatedPointMetrics)

    neumannConditions=>rhsBoundaryConditions%neumannBoundaryConditions
    !Check that Neumann conditions are associated, otherwise do nothing
    IF(ASSOCIATED(neumannConditions)) THEN
      rhsVariable=>rhsBoundaryConditions%VARIABLE
      IF(.NOT.ASSOCIATED(rhsVariable)) THEN
        CALL FLAG_ERROR("Field variable for RHS boundary conditions is not associated.",err,error,*999)
      END IF
      geometricField=>rhsVariable%field%GEOMETRIC_FIELD
      IF(.NOT.ASSOCIATED(geometricField)) THEN
        CALL FLAG_ERROR("Geometric field for the rhs field variable is not associated.",err,error,*999)
      END IF

      neumannConditions%integrationMatrix=0.0_DP
      neumannConditions%pointValues=0.0_DP
      neumannConditions%localDofs=0

      numberOfNeumann=rhsBoundaryConditions%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)

      ! Initialise field interpolation parameters for the geometric field, which are required for the
      ! face/line Jacobian and scale factors
      CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,interpolationParameters,err,error,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,scalingParameters,err,error,*999)
      CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999)
      CALL FIELD_INTERPOLATED_POINTS_METRICS_INITIALISE(interpolatedPoints,interpolatedPointMetrics,err,error,*999)

      ! Loop over all components in the RHS variable
      DO componentIdx=1,rhsVariable%NUMBER_OF_COMPONENTS
        ! Get topology for finding faces/lines
        topology=>rhsVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY
        IF(.NOT.ASSOCIATED(topology)) THEN
          CALL FLAG_ERROR("Field component topology is not associated",err,error,*999)
        END IF
        decomposition=>rhsVariable%COMPONENTS(componentIdx)%DOMAIN%DECOMPOSITION
        IF(.NOT.ASSOCIATED(decomposition)) THEN
          CALL FLAG_ERROR("Field component decomposition is not associated",err,error,*999)
        END IF

        SELECT CASE(rhsVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          SELECT CASE(rhsVariable%COMPONENTS(componentIdx)%DOMAIN%NUMBER_OF_DIMENSIONS)
          CASE(3)
            ! 3D, so integrating faces
            ! Loop over all faces in the decomposition, ignoring those not on the boundary
            IF(.NOT.decomposition%CALCULATE_FACES) THEN
              CALL FLAG_ERROR("Decomposition does not have faces calculated.",err,error,*999)
            END IF
            faces=>topology%FACES
            IF(.NOT.ASSOCIATED(faces)) THEN
              CALL FLAG_ERROR("Mesh topology faces is not associated.",err,error,*999)
            END IF
            facesLoop: DO faceIdx=1,faces%NUMBER_OF_FACES
              face=>faces%faces(faceIdx)
              IF(.NOT.face%BOUNDARY_FACE) CYCLE

              basis=>face%basis
              IF(.NOT.ASSOCIATED(basis)) THEN
                CALL FLAG_ERROR("Face basis is not associated.",err,error,*999)
              END IF

              ! Check if any of the dofs in this face have a Neumann condition set, and keep
              ! track of their dof numbers as well as local element node and derivative numbers
              faceNumberOfNeumann=0
              DO nodeIdx=1,basis%NUMBER_OF_NODES
                localNode=face%NODES_IN_FACE(nodeIdx)
                DO derivIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                  version=face%DERIVATIVES_IN_FACE(2,derivIdx,nodeIdx)
                  localDof=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & NODES(localNode)%DERIVATIVES(derivIdx)%VERSIONS(version)
                  globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDof)
                  IF(rhsBoundaryConditions%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT) THEN
                    faceNumberOfNeumann=faceNumberOfNeumann+1
                    IF(faceNumberOfNeumann>16) CALL FLAG_ERROR( &
                      & "Got more than the expected maximum number of face Neumann DOFs.",err,error,*999)
                    ! Find the Neumann condition index that matches the global dof
                    neumannConditionNumber(faceNumberOfNeumann)=0
                    DO neumannConditionIdx=1,numberOfNeumann
                      IF(neumannConditions%setDofs(neumannConditionIdx)==globalDof) THEN
                        neumannConditionNumber(faceNumberOfNeumann)=neumannConditionIdx
                        CYCLE
                      END IF
                    END DO
                    IF(neumannConditionNumber(faceNumberOfNeumann)==0) THEN
                      CALL FLAG_ERROR("Could not find matching Neumann condition index for local DOF "// &
                        & TRIM(NUMBER_TO_VSTRING(localDof,"*",err,error))//".",err,error,*999)
                    END IF
                    neumannNodeAndDeriv(faceNumberOfNeumann,1)=nodeIdx
                    neumannNodeAndDeriv(faceNumberOfNeumann,2)=derivIdx
                    ! We will also update the point values vector while we know the local DOF number, and
                    ! keep track of the local DOF numbers.
                    ! Check that this isn't a ghost DOF so that we can later do an MPI reduce to
                    ! get all point values.
                    ! This will be called multiple times for DOFs in multiple faces, but that's OK
                    IF(localDof<=rhsVariable%DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE, &
                        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDof,pointValue,err,error,*999)
                      neumannConditions%pointValues( &
                        & neumannConditionNumber(faceNumberOfNeumann))=pointValue
                      neumannConditions%localDofs( &
                        & neumannConditionNumber(faceNumberOfNeumann))=localDof
                    END IF
                  ELSE IF(rhsBoundaryConditions%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED) THEN
                    ! If there is an integrated Neumann value on a DOF in this face,
                    ! then we don't integrate any other point values over this face.
                    ! If this isn't sufficient then we might need to provide some other
                    ! way to specify that a face shouldn't be integrated over
                    CYCLE facesLoop
                  END IF
                END DO
              END DO
              IF(faceNumberOfNeumann==0) CYCLE

              quadratureScheme=>basis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                CALL FLAG_ERROR("Face basis default quadrature scheme is not associated.",err,error,*999)
              END IF
              CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,face%number, &
                & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(face%ELEMENT_NUMBER, &
                  & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              END IF

              ! Calculate the face integrals, then add the calculated values into the integration matrix N
              faceBasisNodesLoop: DO nodeIdx=1,basis%NUMBER_OF_NODES
                localNode=face%NODES_IN_FACE(nodeIdx)
                DO derivIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                  version=face%DERIVATIVES_IN_FACE(2,derivIdx,nodeIdx)
                  localDof=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & NODES(localNode)%DERIVATIVES(derivIdx)%VERSIONS(version)
                  IF(localDof>rhsVariable%DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                    ! We might have a face in this domain but nodes in the face that are ghosted
                    CYCLE faceBasisNodesLoop
                  END IF
                  ! localDof is the weighting DOF
                  ! Loop over all dofs on this face with a Neumann condition set
                  DO elementNeumannDofIdx=1,faceNumberOfNeumann
                    ! Calculate the integral term
                    ! Calculates the term: S(m,a) * S(o,y) * face_integral( phi_m^a phi_o^y )
                    ! Where m and o are the local node and derivative indices for the weighting DOF
                    ! and o and y are the local node and derivative indices for a DOF with a Neumann condition set

                    neumannNode=neumannNodeAndDeriv(elementNeumannDofIdx,1)
                    neumannDeriv=neumannNodeAndDeriv(elementNeumannDofIdx,2)

                    ! Get element parameter local dofs
                    ms=basis%ELEMENT_PARAMETER_INDEX(derivIdx,nodeIdx)
                    os=basis%ELEMENT_PARAMETER_INDEX(neumannDeriv,neumannNode)

                    integratedValue=0.0_DP
                    ! Loop over face gauss points, adding gauss weighted terms to the integral
                    DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                      CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                        & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                      CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

                      !Get basis function values at guass points
                      phim=quadratureScheme%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussIdx)
                      phio=quadratureScheme%GAUSS_BASIS_FNS(os,NO_PART_DERIV,gaussIdx)

                      !Add gauss point value to total face integral
                      integratedValue=integratedValue+phim*phio* &
                        & quadratureScheme%GAUSS_WEIGHTS(gaussIdx)* &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian
                    END DO

                    ! Multiply by scale factors
                    IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                      integratedValue=integratedValue* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(ms,componentIdx)* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(os,componentIdx)
                    END IF

                    ! Add integral term to N matrix
                    matrixNeumannDofIdx=neumannConditionNumber(elementNeumannDofIdx)
                    neumannConditions%integrationMatrix(localDof,matrixNeumannDofIdx)= &
                      & neumannConditions%integrationMatrix(localDof,matrixNeumannDofIdx)+integratedValue
                  END DO !face Neumann points
                END DO !derivIdx
              END DO faceBasisNodesLoop
            END DO facesLoop
          CASE(2)
            ! 2D, so integrating lines
            ! Loop over all lines in the decomposition, ignoring those not on the boundary
            IF(.NOT.decomposition%CALCULATE_LINES) THEN
              CALL FLAG_ERROR("Decomposition does not have lines calculated.",err,error,*999)
            END IF
            lines=>topology%LINES
            IF(.NOT.ASSOCIATED(lines)) THEN
              CALL FLAG_ERROR("Mesh topology lines is not associated.",err,error,*999)
            END IF
            linesLoop: DO lineIdx=1,lines%NUMBER_OF_LINES
              line=>lines%lines(lineIdx)
              IF(.NOT.line%BOUNDARY_LINE) CYCLE

              basis=>line%basis
              IF(.NOT.ASSOCIATED(basis)) THEN
                CALL FLAG_ERROR("Line basis is not associated.",err,error,*999)
              END IF

              ! Check if any of the dofs in this line have a Neumann condition set, and keep
              ! track of their dof numbers as well as local element node and derivative numbers
              lineNumberOfNeumann=0
              DO nodeIdx=1,basis%NUMBER_OF_NODES
                localNode=line%NODES_IN_LINE(nodeIdx)
                DO derivIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                  version=line%DERIVATIVES_IN_LINE(2,derivIdx,nodeIdx)
                  localDof=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & NODES(localNode)%DERIVATIVES(derivIdx)%VERSIONS(version)
                  globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDof)
                  IF(rhsBoundaryConditions%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT) THEN
                    lineNumberOfNeumann=lineNumberOfNeumann+1
                    IF(lineNumberOfNeumann>16) CALL FLAG_ERROR( &
                      & "Got more than the expected maximum number of line Neumann DOFs.",err,error,*999)
                    ! Find the Neumann condition index that matches the global dof
                    neumannConditionNumber(lineNumberOfNeumann)=0
                    DO neumannConditionIdx=1,numberOfNeumann
                      IF(neumannConditions%setDofs(neumannConditionIdx)==globalDof) THEN
                        neumannConditionNumber(lineNumberOfNeumann)=neumannConditionIdx
                        CYCLE
                      END IF
                    END DO
                    IF(neumannConditionNumber(lineNumberOfNeumann)==0) THEN
                      CALL FLAG_ERROR("Could not find matching Neumann condition index for local DOF "// &
                        & TRIM(NUMBER_TO_VSTRING(localDof,"*",err,error))//".",err,error,*999)
                    END IF
                    neumannNodeAndDeriv(lineNumberOfNeumann,1)=nodeIdx
                    neumannNodeAndDeriv(lineNumberOfNeumann,2)=derivIdx
                    ! We will also update the point values vector while we know the local DOF number, and
                    ! keep track of the local DOF numbers.
                    ! Check that this isn't a ghost DOF so that we can later do an MPI reduce to
                    ! get all point values.
                    ! This will be called multiple times for DOFs in multiple lines, but that's OK
                    IF(localDof<=rhsVariable%DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE, &
                        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDof,pointValue,err,error,*999)
                      neumannConditions%pointValues( &
                        & neumannConditionNumber(lineNumberOfNeumann))=pointValue
                      neumannConditions%localDofs( &
                        & neumannConditionNumber(lineNumberOfNeumann))=localDof
                    END IF
                  ELSE IF(rhsBoundaryConditions%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED) THEN
                    ! If there is an integrated Neumann value on a DOF in this line,
                    ! then we don't integrate any other point values over this line.
                    CYCLE linesLoop
                  END IF
                END DO
              END DO
              IF(lineNumberOfNeumann==0) CYCLE

              quadratureScheme=>basis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                CALL FLAG_ERROR("Line basis default quadrature scheme is not associated.",err,error,*999)
              END IF
              CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,lineIdx, &
                & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(line%ELEMENT_NUMBER, &
                  & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              END IF

              ! Calculate the line integrals, then add the calculated values into the integration matrix N
              lineBasisNodesLoop: DO nodeIdx=1,basis%NUMBER_OF_NODES
                localNode=line%NODES_IN_LINE(nodeIdx)
                DO derivIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                  version=line%DERIVATIVES_IN_LINE(2,derivIdx,nodeIdx)
                  localDof=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & NODES(localNode)%DERIVATIVES(derivIdx)%VERSIONS(version)
                  ! localDof is the weighting DOF
                  IF(localDof>rhsVariable%DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                    ! We might have a face in this domain but nodes in the line that are ghosted
                    CYCLE lineBasisNodesLoop
                  END IF
                  ! Loop over all dofs on this line with a Neumann condition set
                  DO elementNeumannDofIdx=1,lineNumberOfNeumann
                    ! Calculate the integral term
                    ! Calculates the term: S(m,a) * S(o,y) * line_integral( phi_m^a phi_o^y )
                    ! Where m and o are the local node and derivative indices for the weighting DOF
                    ! and o and y are the local node and derivative indices for a DOF with a Neumann condition set

                    neumannNode=neumannNodeAndDeriv(elementNeumannDofIdx,1)
                    neumannDeriv=neumannNodeAndDeriv(elementNeumannDofIdx,2)

                    ! Get element parameter local dofs
                    ms=basis%ELEMENT_PARAMETER_INDEX(derivIdx,nodeIdx)
                    os=basis%ELEMENT_PARAMETER_INDEX(neumannDeriv,neumannNode)

                    integratedValue=0.0_DP
                    ! Loop over line gauss points, adding gauss weighted terms to the integral
                    DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                      CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                        & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                      CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_LINE_TYPE, &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

                      !Get basis function values at guass points
                      phim=quadratureScheme%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussIdx)
                      phio=quadratureScheme%GAUSS_BASIS_FNS(os,NO_PART_DERIV,gaussIdx)

                      !Add gauss point value to total line integral
                      integratedValue=integratedValue+phim*phio* &
                        & quadratureScheme%GAUSS_WEIGHTS(gaussIdx)* &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian
                    END DO

                    ! Multiply by scale factors
                    IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                      integratedValue=integratedValue* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(ms,componentIdx)* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(os,componentIdx)
                    END IF

                    ! Add integral term to N matrix
                    matrixNeumannDofIdx=neumannConditionNumber(elementNeumannDofIdx)
                    neumannConditions%integrationMatrix(localDof,matrixNeumannDofIdx)= &
                      & neumannConditions%integrationMatrix(localDof,matrixNeumannDofIdx)+integratedValue
                  END DO !line Neumann points
                END DO !derivIdx
              END DO lineBasisNodesLoop
            END DO linesLoop
          CASE DEFAULT
            CALL FLAG_ERROR("The dimension is invalid for point Neumann conditions",err,error,*999)
          END SELECT !number of dimensions
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE DEFAULT
          CALL FLAG_ERROR("The interpolation type of "// &
            & TRIM(NUMBER_TO_VSTRING(rhsVariable%COMPONENTS(componentIdx) &
            & %INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
            & TRIM(NUMBER_TO_VSTRING(componentIdx,"*",ERR,ERROR))//".", &
            & err,error,*999)
        END SELECT
      END DO !rhs field variable components

      ! Update remaining point values from other computational nodes
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,neumannConditions%pointValues, &
        & numberOfNeumann,MPI_DOUBLE,MPI_SUM, &
        & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,mpiError)
      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",mpiError,err,error,*999)

      ! Perform matrix multiplication, f = N q, to calculate force vector from integration matrix and point values
      neumannConditions%integratedValues = MATMUL(neumannConditions%integrationMatrix,neumannConditions%pointValues)

      ! Update RHS field variable with integrated values, these will be later transferred to the solver RHS vector
      DO localDof=1,rhsVariable%DOMAIN_MAPPING%NUMBER_OF_LOCAL
        globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDof)
        ! Check that the DOF doesn't have an integrated value already set
        IF(rhsBoundaryConditions%CONDITION_TYPES(globalDof)/=BOUNDARY_CONDITION_NEUMANN_INTEGRATED) THEN
          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,localDof,neumannConditions%integratedValues(localDof),err,error,*999)
        END IF
      END DO
      CALL FIELD_PARAMETER_SET_UPDATE_START(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & err,error,*999)
      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNeumann,6,6,neumannConditions%setDofs, &
          & '("  setDofs:",6(X,I8))', '(10X,6(X,I8))',err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point values",err,error,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNeumann,6,6,neumannConditions%pointValues, &
          & '("    p:",6(X,E13.6))', '(6X,6(X,E13.6))',err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann integration matrix",err,error,*999)
        CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,rhsBoundaryConditions%VARIABLE%NUMBER_OF_DOFS, &
          & 1,1,numberOfNeumann,6,6,neumannConditions%integrationMatrix,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("    N','(",I2,",:)',':",6(X,E13.6))','(12X,6(X,E13.6))',err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Integrated values",err,error,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,rhsBoundaryConditions%VARIABLE%NUMBER_OF_DOFS,6,6, &
          & neumannConditions%integratedValues,'("    f:",6(X,E13.6))', '(6X,6(X,E13.6))',err,error,*999)
      END IF
      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & err,error,*999)

    END IF !Neumann conditions associated

    CALL EXITS("BoundaryConditions_NeumannIntegrate")
    RETURN
999 CALL ERRORS("BoundaryConditions_NeumannIntegrate",err,error)
    CALL EXITS("BoundaryConditions_NeumannIntegrate")
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannIntegrate

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
      CALL BoundaryConditions_NeumannFinalise(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
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
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%neumannBoundaryConditions)
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
