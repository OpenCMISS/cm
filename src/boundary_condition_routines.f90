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
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE INTERFACE_CONDITIONS_CONSTANTS
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

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_DOFTypes BOUNDARY_CONDITIONS_ROUTINES::DOFTypes
  !> \brief DOF type for boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FREE=0 !<The dof is free. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FIXED=1 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_MIXED=2 !<The dof is set as a mixed boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_CONSTRAINED=3 !<The dof is constrained to be a linear combination of other DOFs. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
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
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY=20 !<A Neumann integrated boundary condition, and no point values will be integrated over a face or line that includes this dof. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_LINEAR_CONSTRAINT=21 !<The dof is constrained to be a linear combination of other DOFs. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED=22!<A Neumann point boundary condition that is incremented inside a load increment control loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_FITTED=23 !<The dof is fixed as a boundary condition to be updated from fitting data \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_NONREFLECTING=24 !<The dof is fixed and set to a non-reflecting type for 1D wave propagation problems. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_CELLML=25 !<The dof is fixed and set to values specified based on the coupled CellML solution at the dof. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_STREE=26 !<The dof is fixed and set to values specified based on the transmission line theory at the dof. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  !>@}

  INTEGER(INTG), PARAMETER :: MAX_BOUNDARY_CONDITION_NUMBER=26 !The maximum boundary condition type identifier, used for allocating an array with an entry for each type

  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_SparsityTypes BOUNDARY_CONDITIONS_ROUTINES::BoundaryConditions
  !> \brief Storage type for matrices used by boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_SPARSE_MATRICES=1 !<The matrices are stored as sparse matrices.
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FULL_MATRICES=2 !<The matrices are stored as full matrices.
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

  PUBLIC BOUNDARY_CONDITION_DOF_FREE,BOUNDARY_CONDITION_DOF_FIXED,BOUNDARY_CONDITION_DOF_MIXED,BOUNDARY_CONDITION_DOF_CONSTRAINED

  PUBLIC BOUNDARY_CONDITION_FREE,BOUNDARY_CONDITION_FIXED,BOUNDARY_CONDITION_FIXED_INLET,&
    & BOUNDARY_CONDITION_FIXED_OUTLET,BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL,BOUNDARY_CONDITION_FREE_WALL,&
    & BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_DIRICHLET,BOUNDARY_CONDITION_NEUMANN_POINT, &
    & BOUNDARY_CONDITION_CAUCHY,BOUNDARY_CONDITION_ROBIN,BOUNDARY_CONDITION_FIXED_INCREMENTED,BOUNDARY_CONDITION_PRESSURE,&
    & BOUNDARY_CONDITION_PRESSURE_INCREMENTED,BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED, &
    & BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE,BOUNDARY_CONDITION_IMPERMEABLE_WALL,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
    & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED,BOUNDARY_CONDITION_FIXED_STREE, &
    & BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML

  PUBLIC BOUNDARY_CONDITION_SPARSE_MATRICES,BOUNDARY_CONDITION_FULL_MATRICES

  PUBLIC BOUNDARY_CONDITIONS_CREATE_FINISH,BOUNDARY_CONDITIONS_CREATE_START,BOUNDARY_CONDITIONS_DESTROY
  
  PUBLIC BOUNDARY_CONDITIONS_ADD_CONSTANT,BOUNDARY_CONDITIONS_ADD_LOCAL_DOF,BOUNDARY_CONDITIONS_ADD_ELEMENT, &
    & BOUNDARY_CONDITIONS_ADD_NODE,BOUNDARY_CONDITIONS_VARIABLE_GET

  PUBLIC BOUNDARY_CONDITIONS_SET_CONSTANT,BOUNDARY_CONDITIONS_SET_LOCAL_DOF,BOUNDARY_CONDITIONS_SET_ELEMENT, &
    & BOUNDARY_CONDITIONS_SET_NODE,BoundaryConditions_NeumannIntegrate,BoundaryConditions_NeumannSparsityTypeSet

  PUBLIC BoundaryConditions_ConstrainNodeDofsEqual

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

    ENTERS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
          IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>0) THEN
            !Transfer all the boundary conditions to all the computational nodes.
            !\todo Look at this.
            DO variable_idx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
              BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
                FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(VARIABLE_DOMAIN_MAPPING)) THEN
                    SEND_COUNT=VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                    !\todo This operation is a little expensive as we are doing an unnecessary sum across all the ranks in order to combin
                    !\todo the data from each rank into all ranks. We will see how this goes for now.
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%DOF_TYPES, &
                      & SEND_COUNT,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES, &
                      & SEND_COUNT,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="Field variable domain mapping is not associated for variable type "// &
                      & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
                  IF(BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)>0.OR. &
                      & BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
                    CALL BoundaryConditions_NeumannInitialise(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                  END IF

                  ! Loop over all global DOFs, keeping track of the dof indices of specific BC types where required
                  pressureIdx=1
                  neumannIdx=1
                  DO dof_idx=1,FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                    IF(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES(dof_idx)== BOUNDARY_CONDITION_PRESSURE_INCREMENTED) THEN
                      BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED%PRESSURE_INCREMENTED_DOF_INDICES(pressureIdx)=dof_idx
                      pressureIdx=pressureIdx+1
                    ELSE IF(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES(dof_idx)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                        & BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES(dof_idx)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      BOUNDARY_CONDITION_VARIABLE%neumannBoundaryConditions%setDofs(neumannIdx)=dof_idx
                      neumannIdx=neumannIdx+1
                    END IF
                  END DO

                  ! Now that we know where Neumann point DOFs are, we can calculate matrix structure
                  IF(BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)>0.OR. &
                      & BOUNDARY_CONDITION_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
                    CALL BoundaryConditions_NeumannMatricesInitialise(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
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

                      !Store Dirichlet dof indices
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
                                    !Iterate through equations matrices
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
                                          CALL FlagError("Not implemented for column major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                          CALL FlagError("Not implemented for row major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                          !Get Sparsity pattern, number of non zeros, number of rows
                                          CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATION_MATRIX%MATRIX,ROW_INDICES, &
                                            & COLUMN_INDICES,ERR,ERROR,*999)
                                          CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET(EQUATION_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
                                            & ERR,ERROR,*999)
                                          !Get the matrix stored as a linked list
                                          CALL DISTRIBUTED_MATRIX_LINKLIST_GET(EQUATION_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                                          NUMBER_OF_ROWS=EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS
                                          !Initialise sparsity indices arrays
                                          CALL BoundaryConditions_SparsityIndicesInitialise(BOUNDARY_CONDITIONS_DIRICHLET% &
                                            & LINEAR_SPARSITY_INDICES(equations_set_idx,equ_matrix_idx)%PTR, &
                                            & BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS,ERR,ERROR,*999)
                                          !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                                          NULLIFY(SPARSITY_INDICES)
                                          SPARSITY_INDICES=>BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES( &
                                            & equations_set_idx,equ_matrix_idx)%PTR
                                          IF(ASSOCIATED(SPARSITY_INDICES)) THEN
                                            !Setup list for storing dirichlet non zero indices
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
                                              CALL LinkedList_to_Array(list(DIRICHLET_DOF),column_array,ERR,ERROR,*999)
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
                                              CALL LINKEDLIST_DESTROY(list(col_idx),ERR,ERROR,*999)
                                            ENDDO
                                          ELSE
                                            LOCAL_ERROR="Sparsity indices arrays are not associated for this equations matrix."
                                            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                          ENDIF
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                          CALL FlagError("Not implemented for compressed column storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                          CALL FlagError("Not implemented for row column storage.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE,"*",ERR,ERROR)) &
                                            //" is invalid."
                                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                      ELSE
                                        CALL FlagError("The equation matrix is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ENDDO
                                  ENDIF

                                  DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
                                  IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                                    !Iterate through equations matrices
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
                                          CALL FlagError("Not implemented for column major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                          CALL FlagError("Not implemented for row major storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                          !Get Sparsity pattern, number of non zeros, number of rows
                                          CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATION_MATRIX%MATRIX,ROW_INDICES, &
                                            & COLUMN_INDICES,ERR,ERROR,*999)
                                          CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET(EQUATION_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
                                            & ERR,ERROR,*999)
                                          !Sparse matrix in a list
                                          CALL DISTRIBUTED_MATRIX_LINKLIST_GET(EQUATION_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                                          NUMBER_OF_ROWS=EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS
                                          !Intialise sparsity indices arrays
                                          CALL BoundaryConditions_SparsityIndicesInitialise(BOUNDARY_CONDITIONS_DIRICHLET% &
                                            & DYNAMIC_SPARSITY_INDICES(equations_set_idx,equ_matrix_idx)%PTR, &
                                            & BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS,ERR,ERROR,*999)
                                          !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
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
                                              !Dirichlet columns
                                              DIRICHLET_DOF=BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)
                                              CALL LinkedList_to_Array(list(DIRICHLET_DOF),column_array,ERR,ERROR,*999)
                                              !The row indices
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
                                              CALL LINKEDLIST_DESTROY(list(col_idx),ERR,ERROR,*999)
                                            ENDDO
                                          ELSE
                                            LOCAL_ERROR="Sparsity indices arrays are not associated for this equations matrix."
                                            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                          ENDIF
                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                          CALL FlagError("Not implemented for compressed column storage.",ERR,ERROR,*999)
                                        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                          CALL FlagError("Not implemented for row column storage.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE,"*",ERR,ERROR)) &
                                            //" is invalid."
                                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                      ELSE
                                        CALL FlagError("The equation matrix is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ENDDO
                                  ENDIF
                                ELSE
                                  LOCAL_ERROR="Equations Matrices is not associated for these Equations."
                                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                LOCAL_ERROR="Equations is not associated for this Equations Set."
                                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              LOCAL_ERROR="Equations Set is not associated for boundary conditions variable "// &
                                & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !equations_set_idx
                          !\todo Update interface sparsity structure calculate first then update code below.
!                          !Loop over interface conditions. Note that only linear interface matrices implemented so far.
!                          DO interface_condition_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
!                            INTERFACE_CONDITION=>SOLVER_EQUATIONS%SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
!                            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
!                              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
!                              IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
!                                INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
!                                IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
!                                  !Iterate through equations matrices
!                                  DO interface_matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
!                                    INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR
!                                    IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
!                                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(INTERFACE_MATRIX%MATRIX,STORAGE_TYPE,ERR,ERROR,*999)
!                                      SELECT CASE(STORAGE_TYPE)
!                                      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
!                                        !Do nothing
!                                      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
!                                        !Do nothing
!                                      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
!                                        CALL FlagError("Not implemented for column major storage.",ERR,ERROR,*999)
!                                      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
!                                        CALL FlagError("Not implemented for row major storage.",ERR,ERROR,*999)
!                                      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
!                                        !Get Sparsity pattern, number of non zeros, number of rows
!                                        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(INTERFACE_MATRIX%MATRIX,ROW_INDICES, &
!                                          & COLUMN_INDICES,ERR,ERROR,*999)
!                                        CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET(INTERFACE_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
!                                          & ERR,ERROR,*999)
!                                        !Get the matrix stored as a linked list
!                                        CALL DISTRIBUTED_MATRIX_LINKLIST_GET(INTERFACE_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
!                                        NUMBER_OF_ROWS=EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS
!                                        !Initialise sparsity indices arrays
!                                        CALL BoundaryConditions_SparsityIndicesInitialise(BOUNDARY_CONDITIONS_DIRICHLET% &
!                                          & LINEAR_SPARSITY_INDICES(interface_condition_idx,interface_matrix_idx)%PTR, &
!                                          & BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS,ERR,ERROR,*999)
!                                        !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
!                                        NULLIFY(SPARSITY_INDICES)
!                                        SPARSITY_INDICES=>BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES( &
!                                            & interface_condition_idx,interface_matrix_idx)%PTR
!                                        IF(ASSOCIATED(SPARSITY_INDICES)) THEN
!                                          !Setup list for storing dirichlet non zero indices
!                                          NULLIFY(SPARSE_INDICES)
!                                          CALL LIST_CREATE_START(SPARSE_INDICES,ERR,ERROR,*999)
!                                          CALL LIST_DATA_TYPE_SET(SPARSE_INDICES,LIST_INTG_TYPE,ERR,ERROR,*999)
!                                          CALL LIST_INITIAL_SIZE_SET(SPARSE_INDICES, &
!                                            & NUMBER_OF_DIRICHLET_CONDITIONS*(NUMBER_OF_NON_ZEROS/NUMBER_OF_ROWS),ERR,ERROR,*999)
!                                          CALL LIST_CREATE_FINISH(SPARSE_INDICES,ERR,ERROR,*999)
!                                          COUNT=0
!                                          SPARSITY_INDICES%SPARSE_COLUMN_INDICES(1)=1
!                                          LAST=1
!                                          DO dirichlet_idx=1,BOUNDARY_CONDITION_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
!                                            DIRICHLET_DOF=BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(dirichlet_idx)
!                                            CALL LinkedList_to_Array(list(DIRICHLET_DOF),column_array)
!                                              DO row_idx=1,size(column_array)
!                                                CALL LIST_ITEM_ADD(SPARSE_INDICES,column_array(row_idx),ERR,ERROR,*999)
!                                                COUNT=COUNT+1
!                                                LAST=row_idx+1
!                                              ENDDO
!                                            SPARSITY_INDICES%SPARSE_COLUMN_INDICES(dirichlet_idx+1)=COUNT+1
!                                          ENDDO
!                                          CALL LIST_DETACH_AND_DESTROY(SPARSE_INDICES,DUMMY,SPARSITY_INDICES%SPARSE_ROW_INDICES, &
!                                            & ERR,ERROR,*999)
!                                          DO col_idx =1,NUMBER_OF_ROWS
!                                            CALL LINKEDLIST_DESTROY(list(col_idx))
!                                          ENDDO
!                                        ELSE
!                                          LOCAL_ERROR="Sparsity indices arrays are not associated for this interface matrix."
!                                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
!                                        ENDIF
!                                      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
!                                        CALL FlagError("Not implemented for compressed column storage.",ERR,ERROR,*999)
!                                      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
!                                        CALL FlagError("Not implemented for row column storage.",ERR,ERROR,*999)
!                                      CASE DEFAULT
!                                        LOCAL_ERROR="The storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE,"*",ERR,ERROR)) &
!                                          //" is invalid."
!                                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
!                                      END SELECT
!                                    ELSE
!                                      CALL FlagError("The interface matrix is not associated.",ERR,ERROR,*999)
!                                    ENDIF
!                                  ENDDO
!                                ELSE
!                                  LOCAL_ERROR="Interface matrices is not associated for these interface equations."
!                                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
!                                ENDIF
!                              ELSE
!                                LOCAL_ERROR="Interface equations is not associated for this interface condition."
!                                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
!                              ENDIF
!                            ELSE
!                              LOCAL_ERROR="Interface condition is not associated for boundary conditions variable "// &
!                                & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
!                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
!                            ENDIF
!                          ENDDO !interface_condition_idx
                        ELSE
                          LOCAL_ERROR="Solver equations solver mapping is not associated."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Solver equations is not associated."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Dirichlet Boundary Conditions type is not associated for boundary condition variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                  ! Finish field update
                  DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
                    IF(BOUNDARY_CONDITION_VARIABLE%parameterSetRequired(parameterSetIdx)) THEN
                      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD_VARIABLE%FIELD,FIELD_VARIABLE%VARIABLE_TYPE, &
                        & parameterSetIdx,ERR,ERROR,*999)
                    END IF
                  END DO

                  !Finish creating the boundary conditions DOF constraints
                  CALL BoundaryConditions_DofConstraintsCreateFinish(boundary_condition_variable,err,error,*999)
                ELSE
                  LOCAL_ERROR="Field variable is not associated for variable index "// &
                    & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Boundary conditions variable is not associated for variable index "// &
                    & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//".",ERR,ERROR,*999)
              ENDIF
            ENDDO ! variable_idx

          ENDIF
          !Set the finished flag
          BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.TRUE.
        ELSE
          CALL FlagError("Boundary conditions variables array is not allocated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
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
          CALL FlagError("Boundary condition variable is not associated",ERR,ERROR,*999)
        ENDIF
      ENDDO !variable_idx
    ENDIF
    
    EXITS("BOUNDARY_CONDITIONS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS)) THEN
        CALL FlagError("Boundary conditions are already associated for the solver equations.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
          CALL FlagError("Boundary conditions is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
            !Initialise the boundary conditions
            CALL BOUNDARY_CONDITIONS_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Solver equations solver mapping is not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
          !Return the pointer
          BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_CREATE_START")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      CALL BOUNDARY_CONDITIONS_FINALISE(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_DESTROY")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
        DO variable_idx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
          IF(ASSOCIATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR, &
                & ERR,ERROR,*999)
          ELSE
            CALL FlagError("Boundary conditions variable number "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                  & " is not associated",ERR,ERROR,*999)
          ENDIF
        ENDDO !variable_idx
        NULLIFY(BOUNDARY_CONDITIONS%SOLVER_EQUATIONS%SOLVER%SOLVER_EQUATIONS)
        !BOUNDARY_CONDITIONS%SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED = .FALSE.
        !BOUNDARY_CONDITIONS%SOLVER_EQUATIONS%SOLVER_MAPPING%SOLVER_MAPPING_FINISHED = .FALSE.
        DEALLOCATE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)
      ENDIF
      DEALLOCATE(BOUNDARY_CONDITIONS)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_FINALISE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_FINALISE",ERR,ERROR)
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
    INTEGER(INTG) :: DUMMY_ERR,variable_idx,variable_type,equations_set_idx,interface_condition_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: INTERFACE_RHS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS)) THEN
        CALL FlagError("Boundary conditions is already associated for these solver equations.",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
          ALLOCATE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate boundary conditions.",ERR,ERROR,*999)
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.FALSE.
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES=0
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%SOLVER_EQUATIONS=>SOLVER_EQUATIONS
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
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
                              variable_type=LINEAR_MAPPING%LINEAR_MATRIX_VARIABLE_TYPES(variable_idx)
                              IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(VARIABLE_TYPE)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                                CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                  & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(VARIABLE_TYPE)%VARIABLE,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !variable_idx
                          ELSE
                            CALL FlagError("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
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
                            CALL FlagError("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FlagError("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*", &
                                & ERR,ERROR))//" is invalid."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                        SELECT CASE(EQUATIONS%LINEARITY)
                        CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & DYNAMIC_MAPPING%DYNAMIC_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FlagError("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE(EQUATIONS_NONLINEAR)
                          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & DYNAMIC_MAPPING%DYNAMIC_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FlagError("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*", &
                                & ERR,ERROR))//" is invalid."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE DEFAULT
                        LOCAL_ERROR="The equations time dependence type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations mapping has not been finished.",ERR,ERROR,*998)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations equations mapping is not associated.",ERR,ERROR,*998)
                  ENDIF
                ELSE
                  CALL FlagError("Equations has not been finished.",ERR,ERROR,*998)
                ENDIF
              ELSE
                CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",ERR,ERROR,*998)
            ENDIF
          ENDDO !equations_set_idx
          DO interface_condition_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
            INTERFACE_CONDITION=>SOLVER_EQUATIONS%SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
                  INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                  IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                    IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
                      INTERFACE_CONDITION%BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      !Only linear interface equations implemented at the moment
                      SELECT CASE(INTERFACE_EQUATIONS%TIME_DEPENDENCE)
                      CASE(INTERFACE_CONDITION_STATIC,INTERFACE_CONDITION_QUASISTATIC)
                        SELECT CASE(INTERFACE_EQUATIONS%LINEARITY)
                        CASE(INTERFACE_CONDITION_LINEAR)
                          INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                          IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                            variable_type=INTERFACE_MAPPING%LAGRANGE_VARIABLE_TYPE
                            IF(INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES>0) THEN
                              CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & INTERFACE_MAPPING%LAGRANGE_VARIABLE,ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Interface mapping mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                          INTERFACE_RHS_MAPPING=>INTERFACE_MAPPING%RHS_MAPPING
                          IF(ASSOCIATED(INTERFACE_RHS_MAPPING)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                & INTERFACE_RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*", &
                                & ERR,ERROR))//" is invalid."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE DEFAULT
                        LOCAL_ERROR="The equations time dependence type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Interface mapping has not been finished.",ERR,ERROR,*998)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*998)
                  ENDIF
                ELSE
                  CALL FlagError("Interface equations has not been finished.",ERR,ERROR,*998)
                ENDIF
              ELSE
                CALL FlagError("Interface equations is not associated.",ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FlagError("Interface condition not associated.",ERR,ERROR,*998)
            ENDIF
          ENDDO !interface_condition_idx
        ELSE
          CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_FINALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR)
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
 
    ENTERS("BOUNDARY_CONDITIONS_ADD_CONSTANT",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(DEPENDENT_FIELD_VARIABLE)

    !Note: This routine is for constant interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_ADD_CONSTANT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_ADD_CONSTANT",ERR,ERROR)
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
 
    ENTERS("BOUNDARY_CONDITIONS_SET_CONSTANT",ERR,ERROR,*999)

    !Note: This routine is for constant interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_SET_CONSTANT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_SET_CONSTANT",ERR,ERROR)
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
    
    ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,(/DOF_INDEX/),(/CONDITION/),(/VALUE/), &
        & ERR,ERROR,*999)

    EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS",ERR,ERROR,*999)
    NULLIFY(dependent_variable)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
                          ! For integrated Neumann condition, integration is already done, so set the RHS
                          ! dof value directly
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
                          &  BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The specified boundary condition type for dof index "// &
                            & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" of "// &
                            & TRIM(NUMBER_TO_VSTRING(CONDITIONS(i),"*",ERR,ERROR))//" is invalid."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        LOCAL_ERROR="The local dof of  "//&
                          & TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))//" at dof index "// &
                          & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))// &
                          & " is invalid. The dof should be between 1 and "// &
                          & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !i
                  ELSE
                    LOCAL_ERROR="The size of the dof indices array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                      & ") does not match the size of the values array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The size of the dof indices array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                    & ") does not match the size of the fixed conditions array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(CONDITIONS,1),"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Boundary conditions variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The dependent field variable domain mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated..",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,(/DOF_INDEX/),(/CONDITION/),(/VALUE/), &
      & ERR,ERROR,*999)

    EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1",ERR,ERROR)
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
 
    ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
                            & BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The specified boundary condition type for dof index "// &
                            & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" of "// &
                            & TRIM(NUMBER_TO_VSTRING(CONDITIONS(i),"*",ERR,ERROR))//" is invalid."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        LOCAL_ERROR="The local dof of  "//&
                          & TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))//" at dof index "// &
                          & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))// &
                          & " is invalid. The dof should be between 1 and "// &
                          & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !i
                  ELSE
                    LOCAL_ERROR="The size of the dof indices array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                      & ") does not match the size of the values array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The size of the dof indices array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                    & ") does not match the size of the fixed conditions array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(CONDITIONS,1),"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Boundary conditions variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The dependent field variable domain mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS

  !
  !================================================================================================================================
  !

  !> Checks the boundary condition type and sets the boundary condition type and dof type for the boundary conditions.
  !> Makes sure any field parameter sets required are created, and sets the parameter set required array value.
  SUBROUTINE BoundaryConditions_SetConditionType(boundaryConditionsVariable,globalDof,condition,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: globalDof !<The globalDof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dofType, previousCondition, previousDof

    ENTERS("BoundaryConditions_SetConditionType",err,error,*999)

    ! We won't do much checking here as this is only used internally and everything has been checked for
    ! association already
    ! Don't need to make sure field values set type is available as this will always be there, but need
    ! to make sure any other parameter sets required are.
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_FREE)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_FIXED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
      dofType=BOUNDARY_CONDITION_DOF_CONSTRAINED
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
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%VARIABLE_TYPE, &
        & FIELD_INTEGRATED_NEUMANN_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
      boundaryConditionsVariable%parameterSetRequired(FIELD_INTEGRATED_NEUMANN_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE DEFAULT
      CALL FlagError("The specified boundary condition type for dof number "// &
        & TRIM(NUMBER_TO_VSTRING(globalDof,"*",err,error))//" of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT

    !We have a valid boundary condition type
    !Update condition type counts
    previousCondition=boundaryConditionsVariable%CONDITION_TYPES(globalDof)
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
    previousDof=boundaryConditionsVariable%DOF_TYPES(globalDof)
    IF(dofType==BOUNDARY_CONDITION_DOF_FIXED.AND.previousDof/=BOUNDARY_CONDITION_DOF_FIXED) THEN
      boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS= &
        & boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS+1
    ELSE IF(dofType/=BOUNDARY_CONDITION_DOF_FIXED.AND.previousDof==BOUNDARY_CONDITION_DOF_FIXED) THEN
      boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS= &
        & boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS-1
    END IF

    !Set the boundary condition type and DOF type
    boundaryConditionsVariable%CONDITION_TYPES(globalDof)=condition
    boundaryConditionsVariable%DOF_TYPES(globalDof)=dofType

    EXITS("BoundaryConditions_SetConditionType")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetConditionType",err,error)
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

    ENTERS("BOUNDARY_CONDITIONS_ADD_ELEMENT",ERR,ERROR,*999)

    !Note: this routine is for element based interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_ADD_ELEMENT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_ADD_ELEMENT",ERR,ERROR)
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

    ENTERS("BoundaryConditions_CheckInterpolationType",err,error,*999)

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
        & BOUNDARY_CONDITION_NEUMANN_INTEGRATED, &
        & BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
        & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE DEFAULT
      CALL FlagError("The specified boundary condition type of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT
    IF(.NOT.validCondition) THEN
      CALL FlagError("The specified boundary condition type of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is not valid for the field component "// &
        & "interpolation type of "//TRIM(NUMBER_TO_VSTRING(interpolationType,"*",err,error))//".", &
        & err,error,*999)
    END IF

    EXITS("BoundaryConditions_CheckInterpolationType")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CheckInterpolationType",err,error)
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
    INTEGER(INTG) :: boundaryConditionType,equationsSetIdx,specificationSize
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    LOGICAL :: validEquationsSetFound

    ENTERS("BoundaryConditions_CheckEquations",err,error,*999)

    !Get and check pointers we need
    solverEquations=>boundaryConditionsVariable%BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
    IF(.NOT.ASSOCIATED(solverEquations)) THEN
      CALL FlagError("Boundary conditions solver equations are not associated.",err,error,*999)
    END IF
    solverMapping=>solverEquations%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(solverMapping)) THEN
      CALL FlagError("Solver equations solver mapping is not associated.",err,error,*999)
    END IF

    DO boundaryConditionType=1,MAX_BOUNDARY_CONDITION_NUMBER
      !Check if any DOFs have been set to this BC type
      IF(boundaryConditionsVariable%DOF_COUNTS(boundaryConditionType)>0) THEN
        validEquationsSetFound=.FALSE.
        DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
          equationsSet=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%PTR
          IF(.NOT.ASSOCIATED(equationsSet)) THEN
            CALL FlagError("Solver equations equations set is not associated.",err,error,*999)
          END IF
          IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
            CALL FlagError("Equations set specification is not allocated.",err,error,*999)
          END IF
          specificationSize=SIZE(equationsSet%specification,1)

          SELECT CASE(boundaryConditionType)
          CASE(BOUNDARY_CONDITION_FREE)
            ! Valid for any equations set
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_INLET, &
              & BOUNDARY_CONDITION_FIXED_OUTLET)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                  & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
                validEquationsSetFound=.TRUE.
              END IF
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL, &
              & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED,BOUNDARY_CONDITION_FREE_WALL)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                  & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_DARCY_EQUATION_TYPE)) THEN
                validEquationsSetFound=.TRUE.
              ELSE IF(specificationSize==3) THEN
                IF(equationsSet%specification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                    & equationsSet%specification(2)==EQUATIONS_SET_LAPLACE_EQUATION_TYPE.AND. &
                    & equationsSet%specification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
                  validEquationsSetFound=.TRUE.
                END IF
              END IF
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_PRESSURE, &
              & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                & equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
                validEquationsSetFound=.TRUE.
              ELSE IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS .AND. &
                & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE) THEN
                validEquationsSetFound=.TRUE.
              END IF
            ENDIF
          CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
            !Not actually used anywhere? So keep it as invalid, although maybe it should be removed?
            validEquationsSetFound=.FALSE.
          CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
            IF(specificationSize>=3) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                  & equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
                  & (equationsSet%specification(3)==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE.OR. &
                  & equationsSet%specification(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
                  & equationsSet%specification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)) THEN
                validEquationsSetFound=.TRUE.
              END IF
            END IF
          CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_FITTED)
            IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
              & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
              & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
            IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                & (equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_FixedNonreflecting,BOUNDARY_CONDITION_FixedCellml,BOUNDARY_CONDITION_FixedStree)
            IF(equationsSet%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                & (equationsSet%TYPE==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                & equationsSet%TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE DEFAULT
            CALL FlagError("The specified boundary condition type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
              & " is invalid.",err,error,*999)
          END SELECT
        END DO

        IF(.NOT.validEquationsSetFound) THEN
            CALL FlagError("The specified boundary condition type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
              & " is invalid for the equations sets in the solver equations.",err,error,*999)
        END IF
      END IF
    END DO

    EXITS("BoundaryConditions_CheckEquations")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CheckEquations",err,error)
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
 
    ENTERS("BOUNDARY_CONDITIONS_SET_ELEMENT",ERR,ERROR,*999)

    !Note: this routine is for element based interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_SET_ELEMENT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_SET_ELEMENT",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_ADD_NODE",ERR,ERROR,*999)

    NULLIFY(FIELD_VARIABLE)
    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_ADD_NODE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_ADD_NODE",ERR,ERROR)
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
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditions_NeumannInitialise",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      numberOfValues=boundaryConditionsVariable%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)+ &
        & boundaryConditionsVariable%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      ALLOCATE(boundaryConditionsVariable%neumannBoundaryConditions,stat=err)
      IF(err/=0) CALL FlagError("Could not allocate Neumann Boundary Conditions",err,error,*998)
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        NULLIFY(boundaryConditionsNeumann%integrationMatrix)
        NULLIFY(boundaryConditionsNeumann%pointValues)
        NULLIFY(boundaryConditionsNeumann%pointDofMapping)

        numberOfLocalDofs=boundaryConditionsVariable%VARIABLE%NUMBER_OF_DOFS
        ALLOCATE(boundaryConditionsNeumann%setDofs(numberOfValues),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate Neumann set DOFs.",err,error,*999)
        boundaryConditionsNeumann%setDofs=0
      ELSE
        CALL FlagError("The boundary condition Neumann is not associated",err,error,*998)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*998)
    END IF

    EXITS("BoundaryConditions_NeumannInitialise")
    RETURN
999 CALL BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditions_NeumannInitialise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannInitialise

  !
  !================================================================================================================================
  !

  !>Initialise the Neumann boundary conditions matrices and vectors.
  !>This must be done after we know which DOFs have Neumann point conditions so
  !>that we can work out the matrix sparsity pattern.
  SUBROUTINE BoundaryConditions_NeumannMatricesInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann condition matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowMapping, pointDofMapping
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology
    TYPE(DOMAIN_LINE_TYPE), POINTER :: line
    TYPE(DOMAIN_FACE_TYPE), POINTER :: face
    TYPE(LIST_TYPE), POINTER :: columnIndicesList, rowColumnIndicesList
    INTEGER(INTG) :: myComputationalNodeNumber
    INTEGER(INTG) :: numberOfPointDofs, numberNonZeros, numberRowEntries, neumannConditionNumber, localNeumannConditionIdx
    INTEGER(INTG) :: neumannIdx, globalDof, localDof, localDofNyy, domainIdx, numberOfDomains, domainNumber, componentNumber
    INTEGER(INTG) :: nodeIdx, derivIdx, nodeNumber, versionNumber, derivativeNumber, columnNodeNumber, lineIdx, faceIdx, columnDof
    INTEGER(INTG), ALLOCATABLE :: rowIndices(:), columnIndices(:), localDofNumbers(:)
    REAL(DP) :: pointValue
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditions_NeumannMatricesInitialise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      rhsVariable=>boundaryConditionsVariable%variable
      IF(.NOT.ASSOCIATED(rhsVariable)) &
        & CALL FlagError("RHS boundary conditions variable field variable is not associated.",err,error,*999)
      numberOfPointDofs=boundaryConditionsVariable%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT) + &
        & boundaryConditionsVariable%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        ! For rows we can re-use the RHS variable row mapping
        rowMapping=>rhsVariable%DOMAIN_MAPPING
        IF(.NOT.ASSOCIATED(rowMapping)) &
          & CALL FlagError("RHS field variable mapping is not associated.",err,error,*998)

        ! Create a domain mapping for the Neumann point DOFs, required for the distributed matrix columns
        ALLOCATE(pointDofMapping,stat=err)
        IF(err/=0) CALL FlagError("Could not allocate Neumann DOF domain mapping.",err,error,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(pointDofMapping,rowMapping%NUMBER_OF_DOMAINS,err,error,*999)
        boundaryConditionsNeumann%pointDofMapping=>pointDofMapping
        ! Calculate global to local mapping for Neumann DOFs
        pointDofMapping%NUMBER_OF_GLOBAL=numberOfPointDofs
        ALLOCATE(pointDofMapping%GLOBAL_TO_LOCAL_MAP(numberOfPointDofs),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate Neumann point DOF global to local mapping.",err,error,*999)
        ALLOCATE(localDofNumbers(0:rowMapping%NUMBER_OF_DOMAINS-1),stat=err)
        IF(ERR/=0) CALL FlagError("Could not allocate local Neumann DOF numbers.",err,error,*999)
        localDofNumbers=0

        IF(DIAGNOSTICS2) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Local numbering",err,error,*999)
        END IF
        DO neumannIdx=1,numberOfPointDofs
          globalDof=boundaryConditionsNeumann%setDofs(neumannIdx)
          ! Get domain information from the RHS variable domain mapping, but set new local numbers.
          numberOfDomains=rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(globalDof)%NUMBER_OF_DOMAINS
          pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%NUMBER_OF_DOMAINS=numberOfDomains
          ALLOCATE(pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_NUMBER(numberOfDomains),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map local number.",err,error,*999)
          ALLOCATE(pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%DOMAIN_NUMBER(numberOfDomains),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map domain number.",err,error,*999)
          ALLOCATE(pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_TYPE(numberOfDomains),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map local type.",err,error,*999)
          IF(DIAGNOSTICS2) THEN
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point DOF index = ",neumannIdx,err,error,*999)
          END IF
          DO domainIdx=1,numberOfDomains
            domainNumber=rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(globalDof)%DOMAIN_NUMBER(domainIdx)
            pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%DOMAIN_NUMBER(domainIdx)=domainNumber
            pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_TYPE(domainIdx)= &
              & rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_TYPE(domainIdx)
            IF(pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_TYPE(domainIdx)==DOMAIN_LOCAL_INTERNAL.OR. &
                & pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_TYPE(domainIdx)==DOMAIN_LOCAL_BOUNDARY) THEN
              localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
              pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_NUMBER(domainIdx)=localDofNumbers(domainNumber)
              IF(DIAGNOSTICS2) THEN
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global rhs var DOF = ",globalDof,err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainNumber,err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local type = ", &
                  & pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_TYPE(domainIdx),err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local number = ",localDofNumbers(domainNumber),err,error,*999)
              END IF
            ENDIF
          END DO
        END DO
        !Local DOFs must be numbered before ghost DOFs, so loop though again, this time numbering GHOST DOFs
        IF(DIAGNOSTICS2) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Ghost numbering",err,error,*999)
        END IF
        DO neumannIdx=1,numberOfPointDofs
          globalDof=boundaryConditionsNeumann%setDofs(neumannIdx)
          numberOfDomains=rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(globalDof)%NUMBER_OF_DOMAINS
          IF(DIAGNOSTICS2) THEN
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point DOF index = ",neumannIdx,err,error,*999)
          END IF
          DO domainIdx=1,numberOfDomains
            IF(pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_TYPE(domainIdx)==DOMAIN_LOCAL_GHOST) THEN
              domainNumber=rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(globalDof)%DOMAIN_NUMBER(domainIdx)
              localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
              pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_NUMBER(domainIdx)=localDofNumbers(domainNumber)
              IF(DIAGNOSTICS2) THEN
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global rhs var DOF = ",globalDof,err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainNumber,err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local number = ",localDofNumbers(domainNumber),err,error,*999)
              END IF
            ENDIF
          END DO
        END DO

        CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(pointDofMapping,err,error,*999)

        CALL DISTRIBUTED_MATRIX_CREATE_START(rowMapping,pointDofMapping,boundaryConditionsNeumann%integrationMatrix,err,error,*999)
        SELECT CASE(boundaryConditionsVariable%BOUNDARY_CONDITIONS%neumannMatrixSparsity)
        CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
          ! Work out integration matrix sparsity structure
          ! For a single process, compressed column would be more memory efficient, but with
          ! multiple processes the number of Neumann point DOFs could be more than the number
          ! of local row DOFs, and multiplying a compressed row matrix by a vector is faster,
          ! so we will use compressed row storage
          ALLOCATE(rowIndices(rowMapping%TOTAL_NUMBER_OF_LOCAL+1),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann integration matrix column indices.",err,error,*999)
          ! We don't know the number of non zeros before hand, so use a list to keep track of column indices
          NULLIFY(columnIndicesList)
          CALL LIST_CREATE_START(columnIndicesList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(columnIndicesList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_CREATE_FINISH(columnIndicesList,err,error,*999)
          ! Stores the column indices for the current row
          NULLIFY(rowColumnIndicesList)
          CALL LIST_CREATE_START(rowColumnIndicesList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(rowColumnIndicesList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_MUTABLE_SET(rowColumnIndicesList,.TRUE.,err,error,*999)
          CALL LIST_CREATE_FINISH(rowColumnIndicesList,err,error,*999)
          rowIndices(1)=1

          DO localDof=1,rhsVariable%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
            localDofNyy=rhsVariable%DOF_TO_PARAM_MAP%DOF_TYPE(2,localDof)
            componentNumber=rhsVariable%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(4,localDofNyy)
            ! Get topology for finding faces/lines
            topology=>rhsVariable%COMPONENTS(componentNumber)%DOMAIN%TOPOLOGY
            IF(.NOT.ASSOCIATED(topology)) THEN
              CALL FlagError("Field component topology is not associated.",err,error,*999)
            END IF

            SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%INTERPOLATION_TYPE)
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              nodeNumber=rhsVariable%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,localDofNyy)
              IF(.NOT.ASSOCIATED(topology%NODES%NODES)) THEN
                CALL FlagError("Topology nodes are not associated.",err,error,*999)
              END IF
              IF(topology%NODES%NODES(nodeNumber)%BOUNDARY_NODE) THEN
                SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%DOMAIN%NUMBER_OF_DIMENSIONS)
                CASE(1)
                  ! Only one column used, as this is the same as setting an integrated
                  ! value so no other DOFs are affected
                  globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDof)
                  IF(boundaryConditionsVariable%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & boundaryConditionsVariable%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                    ! Find the Neumann condition number
                    neumannConditionNumber=0
                    DO neumannIdx=1,numberOfPointDofs
                      IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDof) THEN
                        neumannConditionNumber=neumannIdx
                      END IF
                    END DO
                    IF(neumannConditionNumber==0) THEN
                      CALL FlagError("Could not find matching Neuamann condition number for global DOF "// &
                        & TRIM(NUMBER_TO_VSTRING(globalDof,"*",err,error))//" with Neumann condition set.",err,error,*999)
                    ELSE
                      CALL LIST_ITEM_ADD(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                    END IF
                  END IF
                CASE(2)
                  ! Loop over all lines for this node and find any DOFs that have a Neumann point condition set
                  DO lineIdx=1,topology%NODES%NODES(nodeNumber)%NUMBER_OF_NODE_LINES
                    IF(.NOT.ALLOCATED(topology%LINES%LINES)) THEN
                      CALL FlagError("Topology lines have not been calculated.",err,error,*999)
                    END IF
                    line=>topology%LINES%LINES(topology%NODES%NODES(nodeNumber)%NODE_LINES(lineIdx))
                    IF(.NOT.line%BOUNDARY_LINE) CYCLE
                    DO nodeIdx=1,line%BASIS%NUMBER_OF_NODES
                      columnNodeNumber=line%NODES_IN_LINE(nodeIdx)
                      DO derivIdx=1,line%BASIS%NUMBER_OF_DERIVATIVES(nodeIdx)
                        derivativeNumber=line%DERIVATIVES_IN_LINE(1,derivIdx,nodeIdx)
                        versionNumber=line%DERIVATIVES_IN_LINE(2,derivIdx,nodeIdx)
                        columnDof=rhsVariable%COMPONENTS(componentNumber)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                          & NODES(columnNodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                        globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(columnDof)
                        IF(boundaryConditionsVariable%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                            & boundaryConditionsVariable%CONDITION_TYPES(globalDof)== &
                            & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                          neumannConditionNumber=0
                          DO neumannIdx=1,numberOfPointDofs
                            IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDof) THEN
                              neumannConditionNumber=neumannIdx
                            END IF
                          END DO
                          IF(neumannConditionNumber==0) THEN
                            CALL FlagError("Could not find matching Neuamann condition number for global DOF "// &
                              & TRIM(NUMBER_TO_VSTRING(globalDof,"*",err,error))//" with Neumann condition set.",err,error,*999)
                          ELSE
                            CALL LIST_ITEM_ADD(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                CASE(3)
                  ! Loop over all faces for this node and find any DOFs that have a Neumann point condition set 
                  DO faceIdx=1,topology%NODES%NODES(nodeNumber)%NUMBER_OF_NODE_FACES
                    IF(.NOT.ALLOCATED(topology%faces%faces)) THEN
                      CALL FlagError("Topology faces have not been calculated.",err,error,*999)
                    END IF
                    face=>topology%FACES%FACES(topology%NODES%NODES(nodeNumber)%NODE_FACES(faceIdx))
                    IF(.NOT.face%BOUNDARY_FACE) CYCLE
                    DO nodeIdx=1,face%BASIS%NUMBER_OF_NODES
                      columnNodeNumber=face%NODES_IN_FACE(nodeIdx)
                      DO derivIdx=1,face%BASIS%NUMBER_OF_DERIVATIVES(nodeIdx)
                        derivativeNumber=face%DERIVATIVES_IN_FACE(1,derivIdx,nodeIdx)
                        versionNumber=face%DERIVATIVES_IN_FACE(2,derivIdx,nodeIdx)
                        columnDof=rhsVariable%COMPONENTS(componentNumber)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                          & NODES(columnNodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                        globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(columnDof)
                        IF(boundaryConditionsVariable%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                            & boundaryConditionsVariable%CONDITION_TYPES(globalDof)== &
                            & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                          neumannConditionNumber=0
                          DO neumannIdx=1,numberOfPointDofs
                            IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDof) THEN
                              neumannConditionNumber=neumannIdx
                            END IF
                          END DO
                          IF(neumannConditionNumber==0) THEN
                            CALL FlagError("Could not find matching Neuamann condition number for global DOF "// &
                              & TRIM(NUMBER_TO_VSTRING(globalDof,"*",err,error))//" with Neumann condition set.",err,error,*999)
                          ELSE
                            CALL LIST_ITEM_ADD(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                CASE DEFAULT
                  CALL FlagError("The dimension is invalid for point Neumann conditions",err,error,*999)
                END SELECT !number of dimensions
              END IF
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              CALL FlagError("The interpolation type of "// &
                & TRIM(NUMBER_TO_VSTRING(rhsVariable%COMPONENTS(componentNumber) &
                & %INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                & TRIM(NUMBER_TO_VSTRING(componentNumber,"*",ERR,ERROR))//".", &
                & err,error,*999)
            END SELECT

            !Sort and remove duplicates
            CALL LIST_REMOVE_DUPLICATES(rowColumnIndicesList,err,error,*999)
            !Now add all column DOFs in this row that use Neumann conditions to the overall column indices
            CALL List_AppendList(columnIndicesList,rowColumnIndicesList,err,error,*999)
            CALL LIST_NUMBER_OF_ITEMS_GET(rowColumnIndicesList,numberRowEntries,err,error,*999)
            rowIndices(localDof+1)=rowIndices(localDof)+numberRowEntries
            CALL List_ClearItems(rowColumnIndicesList,err,error,*999)
          END DO !local DOFs

          CALL LIST_DESTROY(rowColumnIndicesList,err,error,*999)
          CALL LIST_DETACH_AND_DESTROY(columnIndicesList,numberNonZeros,columnIndices,err,error,*999)
          IF(DIAGNOSTICS1) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Neumann integration matrix sparsity",err,error,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number non-zeros = ", numberNonZeros,err,error,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number columns = ",numberOfPointDofs,err,error,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number rows = ", &
              & rhsVariable%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,err,error,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfPointDofs+1,6,6, &
              & rowIndices,'("  Row indices: ",6(X,I6))', '(6X,6(X,I6))',err,error,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberNonZeros,6,6, &
              & columnIndices,'("  Column indices: ",6(X,I6))', '(6X,6(X,I6))',err,error,*999)
          END IF

          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(boundaryConditionsNeumann%integrationMatrix, &
            & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
          CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(boundaryConditionsNeumann%integrationMatrix,numberNonZeros,err,error,*999)
          CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(boundaryConditionsNeumann%integrationMatrix, &
            & rowIndices,columnIndices(1:numberNonZeros),err,error,*999)

          DEALLOCATE(localDofNumbers)
          DEALLOCATE(rowIndices)
          DEALLOCATE(columnIndices)
        CASE(BOUNDARY_CONDITION_FULL_MATRICES)
          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(boundaryConditionsNeumann%integrationMatrix, &
            & DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
        CASE DEFAULT
          CALL FlagError("The Neumann matrix sparsity type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionsVariable%BOUNDARY_CONDITIONS%neumannMatrixSparsity,"*",err,error))// &
              & " is invalid.",err,error,*999)
        END SELECT

        CALL DISTRIBUTED_MATRIX_CREATE_FINISH(boundaryConditionsNeumann%integrationMatrix,err,error,*999)

        !Set up vector of Neumann point values
        CALL DISTRIBUTED_VECTOR_CREATE_START(pointDofMapping,boundaryConditionsNeumann%pointValues,err,error,*999)
        CALL DISTRIBUTED_VECTOR_CREATE_FINISH(boundaryConditionsNeumann%pointValues,err,error,*999)
        myComputationalNodeNumber=COMPUTATIONAL_NODE_NUMBER_GET(err,error)
        !Set point values vector from boundary conditions field parameter set
        DO neumannIdx=1,numberOfPointDofs
          globalDof=boundaryConditionsNeumann%setDofs(neumannIdx)
          IF(rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(globalDof)%DOMAIN_NUMBER(1)==myComputationalNodeNumber) THEN
            localDof=rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_NUMBER(1)
            ! Set point DOF vector value
            localNeumannConditionIdx=boundaryConditionsNeumann%pointDofMapping%GLOBAL_TO_LOCAL_MAP(neumannIdx)%LOCAL_NUMBER(1)
            CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE, &
              & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDof,pointValue,err,error,*999)
            CALL DISTRIBUTED_VECTOR_VALUES_SET(boundaryConditionsNeumann%pointValues, &
              & localNeumannConditionIdx,pointValue,err,error,*999)
          END IF
        END DO
        CALL DISTRIBUTED_VECTOR_UPDATE_START(boundaryConditionsNeumann%pointValues,err,error,*999)
        CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(boundaryConditionsNeumann%pointValues,err,error,*999)

      ELSE
        CALL FlagError("The boundary condition Neumann is not associated",err,error,*998)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*998)
    END IF

    EXITS("BoundaryConditions_NeumannMatricesInitialise")
    RETURN
999 IF(ALLOCATED(rowIndices)) THEN
      DEALLOCATE(rowIndices)
    END IF
    IF(ALLOCATED(columnIndices)) THEN
      DEALLOCATE(columnIndices)
    END IF
    IF(ALLOCATED(localDofNumbers)) THEN
      DEALLOCATE(localDofNumbers)
    END IF
    CALL BoundaryConditions_NeumannMatricesFinalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
998 ERRORS("BoundaryConditions_NeumannMatricesInitialise",err,error)
    EXITS("BoundaryConditions_NeumannMatricesInitialise")
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_NeumannMatricesInitialise

  !
  !================================================================================================================================
  !

  !Finalise the Neumann condition information for a boundary conditions variable
  SUBROUTINE BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise the Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann

    ENTERS("BoundaryConditions_NeumannFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        IF(ALLOCATED(boundaryConditionsNeumann%setDofs)) &
          & DEALLOCATE(boundaryConditionsNeumann%setDofs)
        CALL BoundaryConditions_NeumannMatricesFinalise(boundaryConditionsVariable,err,error,*999)
        DEALLOCATE(boundaryConditionsNeumann)
        NULLIFY(boundaryConditionsVariable%neumannBoundaryConditions)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannFinalise

  !
  !================================================================================================================================
  !

  !Finalise the Neumann condition matrices for a boundary conditions variable
  SUBROUTINE BoundaryConditions_NeumannMatricesFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise Neumann condition matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann

    ENTERS("BoundaryConditions_NeumannMatricesFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        IF(ASSOCIATED(boundaryConditionsNeumann%integrationMatrix)) &
          & CALL DISTRIBUTED_MATRIX_DESTROY(boundaryConditionsNeumann%integrationMatrix,err,error,*999)
        IF(ASSOCIATED(boundaryConditionsNeumann%pointValues)) &
          & CALL DISTRIBUTED_VECTOR_DESTROY(boundaryConditionsNeumann%pointValues,err,error,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(boundaryConditionsNeumann%pointDofMapping,err,error,*999)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannMatricesFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannMatricesFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannMatricesFinalise

  !
  !================================================================================================================================
  !

  !>Calculates integrated Neumann condition values from point values for a boundary conditions variable and
  !>updates the FIELD_INTEGRATED_NEUMANN_SET_TYPE parameter set for the field variable.
  SUBROUTINE BoundaryConditions_NeumannIntegrate(rhsBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER, INTENT(IN) :: rhsBoundaryConditions !<The boundary conditions for the right hand side field variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
    INTEGER(INTG) :: componentNumber,globalDof,localDof,neumannDofIdx,myComputationalNodeNumber
    INTEGER(INTG) :: numberOfNeumann,neumannLocalDof,neumannDofNyy
    INTEGER(INTG) :: neumannGlobalDof,neumannNodeNumber,neumannLocalNodeNumber,neumannLocalDerivNumber
    INTEGER(INTG) :: faceIdx,lineIdx,nodeIdx,derivIdx,gaussIdx
    INTEGER(INTG) :: faceNumber,lineNumber
    INTEGER(INTG) :: ms,os,nodeNumber,derivativeNumber,versionNumber
    LOGICAL :: dependentGeometry
    REAL(DP) :: integratedValue,phim,phio
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannConditions
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointMetrics(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:), scalingParameters(:)
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: integratedValues
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology
    TYPE(DOMAIN_FACES_TYPE), POINTER :: faces
    TYPE(DOMAIN_LINES_TYPE), POINTER :: lines
    TYPE(DOMAIN_FACE_TYPE), POINTER :: face
    TYPE(DOMAIN_LINE_TYPE), POINTER :: line
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme

    ENTERS("BoundaryConditions_NeumannIntegrate",err,error,*999)

    NULLIFY(scalingParameters)
    NULLIFY(interpolationParameters)
    NULLIFY(interpolatedPoints)
    NULLIFY(interpolatedPointMetrics)
    NULLIFY(integratedValues)

    neumannConditions=>rhsBoundaryConditions%neumannBoundaryConditions
    !Check that Neumann conditions are associated, otherwise do nothing
    IF(ASSOCIATED(neumannConditions)) THEN
      rhsVariable=>rhsBoundaryConditions%VARIABLE
      IF(.NOT.ASSOCIATED(rhsVariable)) THEN
        CALL FlagError("Field variable for RHS boundary conditions is not associated.",err,error,*999)
      END IF

      CALL Field_GeometricGeneralFieldGet(rhsVariable%field,geometricField,dependentGeometry,err,error,*999)

      CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(neumannConditions%integrationMatrix,0.0_DP,err,error,*999)

      numberOfNeumann=rhsBoundaryConditions%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT) + &
        & rhsBoundaryConditions%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      myComputationalNodeNumber=COMPUTATIONAL_NODE_NUMBER_GET(err,error)

      ! Initialise field interpolation parameters for the geometric field, which are required for the
      ! face/line Jacobian and scale factors
      CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,interpolationParameters,err,error,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(rhsVariable%field,scalingParameters,err,error,*999)
      CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999)
      CALL Field_InterpolatedPointsMetricsInitialise(interpolatedPoints,interpolatedPointMetrics,err,error,*999)

      ! Loop over all Neumann point DOFs, finding the boundary lines or faces they are on
      ! and integrating over them
      DO neumannDofIdx=1,numberOfNeumann
        neumannGlobalDof=neumannConditions%setDofs(neumannDofIdx)
        IF(rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(neumannGlobalDof)%DOMAIN_NUMBER(1)==myComputationalNodeNumber) THEN
          neumannLocalDof=rhsVariable%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(neumannGlobalDof)%LOCAL_NUMBER(1)
          ! Get Neumann DOF component and topology for that component
          neumannDofNyy=rhsVariable%DOF_TO_PARAM_MAP%DOF_TYPE(2,neumannLocalDof)
          componentNumber=rhsVariable%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(4,neumannDofNyy)
          topology=>rhsVariable%COMPONENTS(componentNumber)%DOMAIN%TOPOLOGY
          IF(.NOT.ASSOCIATED(topology)) THEN
            CALL FlagError("Field component topology is not associated.",err,error,*999)
          END IF
          decomposition=>rhsVariable%COMPONENTS(componentNumber)%DOMAIN%DECOMPOSITION
          IF(.NOT.ASSOCIATED(decomposition)) THEN
            CALL FlagError("Field component decomposition is not associated.",err,error,*999)
          END IF
          SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%INTERPOLATION_TYPE)
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            neumannNodeNumber=rhsVariable%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,neumannDofNyy)
            SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%DOMAIN%NUMBER_OF_DIMENSIONS)
            CASE(1)
              CALL DISTRIBUTED_MATRIX_VALUES_SET(neumannConditions%integrationMatrix,neumannLocalDof,neumannDofIdx, &
                & 1.0_DP,err,error,*999)
            CASE(2)
              IF(.NOT.decomposition%CALCULATE_LINES) THEN
                CALL FlagError("Decomposition does not have lines calculated.",err,error,*999)
              END IF
              lines=>topology%LINES
              IF(.NOT.ASSOCIATED(lines)) THEN
                CALL FlagError("Mesh topology lines is not associated.",err,error,*999)
              END IF
              linesLoop: DO lineIdx=1,topology%NODES%NODES(neumannNodeNumber)%NUMBER_OF_NODE_LINES
                lineNumber=topology%NODES%NODES(neumannNodeNumber)%NODE_LINES(lineIdx)
                line=>topology%lines%lines(lineNumber)
                IF(.NOT.line%BOUNDARY_LINE) &
                  CYCLE linesLoop
                basis=>line%basis
                IF(.NOT.ASSOCIATED(basis)) THEN
                  CALL FlagError("Line basis is not associated.",err,error,*999)
                END IF
                neumannLocalNodeNumber=0
                neumannLocalDerivNumber=0
                ! Check all nodes in line to find the local numbers for the Neumann DOF, and
                ! make sure we don't have an integrated_only condition set on the line
                DO nodeIdx=1,line%BASIS%NUMBER_OF_NODES
                  nodeNumber=line%NODES_IN_LINE(nodeIdx)
                  DO derivIdx=1,line%BASIS%NUMBER_OF_DERIVATIVES(nodeIdx)
                    derivativeNumber=line%DERIVATIVES_IN_LINE(1,derivIdx,nodeIdx)
                    versionNumber=line%DERIVATIVES_IN_LINE(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                    globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDof)
                    IF(globalDof==neumannGlobalDof) THEN
                      neumannLocalNodeNumber=nodeIdx
                      neumannLocalDerivNumber=derivIdx
                    ELSE IF(rhsBoundaryConditions%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE linesLoop
                    END IF
                  END DO
                END DO
                IF(neumannLocalNodeNumber==0) THEN
                  CALL FlagError("Could not find local Neumann node and derivative numbers in line.",err,error,*999)
                END IF

                ! Now perform actual integration
                quadratureScheme=>basis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                  CALL FlagError("Line basis default quadrature scheme is not associated.",err,error,*999)
                END IF
                CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,lineNumber, &
                  & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                  CALL Field_InterpolationParametersScaleFactorsLineGet(lineNumber, &
                    & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                END IF

                DO nodeIdx=1,line%BASIS%NUMBER_OF_NODES
                  nodeNumber=line%NODES_IN_LINE(nodeIdx)
                  DO derivIdx=1,line%BASIS%NUMBER_OF_DERIVATIVES(nodeIdx)
                    derivativeNumber=line%DERIVATIVES_IN_LINE(1,derivIdx,nodeIdx)
                    versionNumber=line%DERIVATIVES_IN_LINE(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)

                    ms=basis%ELEMENT_PARAMETER_INDEX(derivIdx,nodeIdx)
                    os=basis%ELEMENT_PARAMETER_INDEX(neumannLocalDerivNumber,neumannLocalNodeNumber)

                    integratedValue=0.0_DP
                    ! Loop over line gauss points, adding gauss weighted terms to the integral
                    DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                      CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                        & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
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

                    ! Multiply by scale factors for dependent variable
                    IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                      integratedValue=integratedValue* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(ms,componentNumber)* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(os,componentNumber)
                    END IF

                    ! Add integral term to N matrix
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(neumannConditions%integrationMatrix,localDof,neumannDofIdx, &
                      & integratedValue,err,error,*999)
                  END DO
                END DO
              END DO linesLoop
            CASE(3)
              IF(.NOT.decomposition%CALCULATE_FACES) THEN
                CALL FlagError("Decomposition does not have faces calculated.",err,error,*999)
              END IF
              faces=>topology%FACES
              IF(.NOT.ASSOCIATED(faces)) THEN
                CALL FlagError("Mesh topology faces is not associated.",err,error,*999)
              END IF
              facesLoop: DO faceIdx=1,topology%NODES%NODES(neumannNodeNumber)%NUMBER_OF_NODE_FACES 
                faceNumber=topology%NODES%NODES(neumannNodeNumber)%NODE_FACES(faceIdx)
                face=>topology%FACES%FACES(faceNumber)
                IF(.NOT.face%BOUNDARY_FACE) &
                  CYCLE facesLoop
                basis=>face%BASIS
                IF(.NOT.ASSOCIATED(basis)) THEN
                  CALL FlagError("Line face is not associated.",err,error,*999)
                END IF
                neumannLocalNodeNumber=0
                neumannLocalDerivNumber=0
                ! Check all nodes in the face to find the local numbers for the Neumann DOF, and
                ! make sure we don't have an integrated_only condition set on the face
                DO nodeIdx=1,basis%NUMBER_OF_NODES
                  nodeNumber=face%NODES_IN_FACE(nodeIdx)
                  DO derivIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                    derivativeNumber=face%DERIVATIVES_IN_FACE(1,derivIdx,nodeIdx)
                    versionNumber=face%DERIVATIVES_IN_FACE(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                    globalDof=rhsVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDof)
                    IF(globalDof==neumannGlobalDof) THEN
                      neumannLocalNodeNumber=nodeIdx
                      neumannLocalDerivNumber=derivIdx
                    ELSE IF(rhsBoundaryConditions%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE facesLoop
                    END IF
                  END DO
                END DO
                IF(neumannLocalNodeNumber==0) THEN
                  CALL FlagError("Could not find local Neumann node and derivative numbers in line.",err,error,*999)
                END IF

                ! Now perform actual integration
                quadratureScheme=>basis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                  CALL FlagError("Face basis default quadrature scheme is not associated.",err,error,*999)
                END IF
                CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,faceNumber, &
                  & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                  CALL Field_InterpolationParametersScaleFactorsFaceGet(faceNumber, &
                    & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                END IF

                DO nodeIdx=1,basis%NUMBER_OF_NODES
                  nodeNumber=face%NODES_IN_FACE(nodeIdx)
                  DO derivIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                    derivativeNumber=face%DERIVATIVES_IN_FACE(1,derivIdx,nodeIdx)
                    versionNumber=face%DERIVATIVES_IN_FACE(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)

                    ms=basis%ELEMENT_PARAMETER_INDEX(derivIdx,nodeIdx)
                    os=basis%ELEMENT_PARAMETER_INDEX(neumannLocalDerivNumber,neumannLocalNodeNumber)

                    integratedValue=0.0_DP
                    ! Loop over line gauss points, adding gauss weighted terms to the integral
                    DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                      CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                        & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
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
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(ms,componentNumber)* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(os,componentNumber)
                    END IF

                    ! Add integral term to N matrix
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(neumannConditions%integrationMatrix,localDof,neumannDofIdx, &
                      & integratedValue,err,error,*999)
                  END DO
                END DO
              END DO facesLoop
            CASE DEFAULT
              CALL FlagError("The dimension is invalid for point Neumann conditions",err,error,*999)
            END SELECT
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            CALL FlagError("The interpolation type of "// &
              & TRIM(NUMBER_TO_VSTRING(rhsVariable%COMPONENTS(componentNumber) &
              & %INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
              & TRIM(NUMBER_TO_VSTRING(componentNumber,"*",ERR,ERROR))//".", &
              & err,error,*999)
          END SELECT
        END IF
      END DO

      CALL DISTRIBUTED_MATRIX_UPDATE_START(neumannConditions%integrationMatrix,err,error,*999)
      CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(neumannConditions%integrationMatrix,err,error,*999)

      CALL FIELD_PARAMETER_SET_VECTOR_GET(rhsVariable%field,rhsVariable%variable_type,FIELD_INTEGRATED_NEUMANN_SET_TYPE, &
        & integratedValues,err,error,*999)
      CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(integratedValues,0.0_DP,err,error,*999)
      ! Perform matrix multiplication, f = N q, to calculate force vector from integration matrix and point values
      CALL DISTRIBUTED_MATRIX_BY_VECTOR_ADD(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
        & neumannConditions%integrationMatrix,neumannConditions%pointValues,integratedValues, &
        & err,error,*999)

      CALL FIELD_PARAMETER_SET_UPDATE_START(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE,FIELD_INTEGRATED_NEUMANN_SET_TYPE, &
        & err,error,*999)
      IF(DIAGNOSTICS1) THEN
        IF(dependentGeometry) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Using dependent field geometry",err,error,*999)
        ELSE
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Using undeformed geometry",err,error,*999)
        END IF
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNeumann,6,6,neumannConditions%setDofs, &
          & '("  setDofs:",6(X,I8))', '(10X,6(X,I8))',err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point values",err,error,*999)
        CALL DISTRIBUTED_VECTOR_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%pointValues,err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann integration matrix",err,error,*999)
        CALL DISTRIBUTED_MATRIX_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%integrationMatrix,err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Integrated values",err,error,*999)
        CALL DISTRIBUTED_VECTOR_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,integratedValues,err,error,*999)
      END IF
      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE,FIELD_INTEGRATED_NEUMANN_SET_TYPE, &
        & err,error,*999)

    END IF !Neumann conditions associated

    EXITS("BoundaryConditions_NeumannIntegrate")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannIntegrate",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannIntegrate

  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the Neumann integration matrices
  SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet(boundaryConditions,sparsityType,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The matrix sparsity type to be set \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions

    ENTERS("BoundaryConditions_NeumannSparsityTypeSet",ERR,ERROR,*999)

    IF(ASSOCIATED(boundaryConditions)) THEN
      SELECT CASE(sparsityType)
      CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
        boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
      CASE(BOUNDARY_CONDITION_FULL_MATRICES)
        boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_FULL_MATRICES
      CASE DEFAULT
        CALL FlagError("The specified Neumann integration matrix sparsity type of "// &
          & TRIM(NUMBER_TO_VSTRING(sparsityType,"*",err,error))//" is invalid.",err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Boundary conditions are not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannSparsityTypeSet")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannSparsityTypeSet",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet

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

    ENTERS("BOUNDARY_CONDITIONS_SET_NODE",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(FIELD_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
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
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_SET_NODE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_SET_NODE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_NODE
  
  !
  !================================================================================================================================
  !

  !>Constrain multiple equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainDofsEqual(boundaryConditions,fieldVariable,globalDofs,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions for the solver equations in which to constrain the DOF.
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable containing the DOFs.
    INTEGER(INTG), INTENT(IN) :: globalDofs(:) !<The global DOFs to be constrained to be equal.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfDofs,dofIdx,dofIdx2

    ENTERS("BoundaryConditions_ConstrainDofsEqual",err,error,*999)

    numberOfDofs=SIZE(globalDofs,1)
    IF(numberOfDofs<2) THEN
      CALL FlagError("Cannot constrain zero or 1 DOF to be equal.",err,error,*999)
    END IF

    !Check for duplicate DOFs
    DO dofIdx=1,numberOfDofs
      DO dofIdx2=dofIdx+1,numberOfDofs
        IF(globalDofs(dofIdx)==globalDofs(dofIdx2)) THEN
          CALL FlagError("DOF number "//TRIM(NumberToVstring(globalDofs(dofIdx),"*",err,error))// &
            & " is duplicated in the DOFs constrained to be equal.",err,error,*999)
        END IF
      END DO
    END DO

    !Add new DOF constraints
    !We set all DOFs except the first to be equal to 1.0 * the first DOF
    !The first DOF is left unconstrained
    DO dofIdx=2,numberOfDofs
      CALL BoundaryConditions_DofConstraintSet( &
        & boundaryConditions,fieldVariable,globalDofs(dofIdx),[globalDofs(1)],[1.0_DP],err,error,*999)
    END DO

    EXITS("BoundaryConditions_ConstrainDofsEqual")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstrainDofsEqual",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_ConstrainDofsEqual

  !
  !================================================================================================================================
  !

  !>Constrain multiple nodal equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainNodeDofsEqual( &
      & boundaryConditions,field,fieldVariableType,versionNumber,derivativeNumber,component,nodes,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER, INTENT(IN) :: boundaryConditions !<The solver equations boundary conditions to constrain the DOFs for.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The equations dependent field containing the field DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: fieldVariableType !<The field variable type of the DOFs to be constrained. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version number.
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number.
    INTEGER(INTG), INTENT(IN) :: component !<The field component number of the DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: nodes(:) !<The user numbers of the nodes to be constrained to be equal.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    INTEGER(INTG) :: numberOfNodes, nodeIdx, dof
    INTEGER(INTG), ALLOCATABLE :: globalDofs(:)

    ENTERS("BoundaryConditions_ConstrainNodeDofsEqual",err,error,*998)

    NULLIFY(fieldVariable)

    IF(.NOT.ASSOCIATED(boundaryConditions)) THEN
      CALL FlagError("Boundary conditions are not associated.",err,error,*998)
    END IF

    numberOfNodes=SIZE(nodes,1)
    ALLOCATE(globalDofs(numberOfNodes),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equal global DOFs array.",err,error,*998)
    !Get field DOFs for nodes
    DO nodeIdx=1,numberOfNodes
      CALL FIELD_COMPONENT_DOF_GET_USER_NODE(field,fieldVariableType,versionNumber,derivativeNumber,nodes(nodeIdx), &
        & component,dof,globalDofs(nodeIdx),err,error,*999)
    END DO
    !Get the field variable and boundary conditions variable for the field
    CALL FIELD_VARIABLE_GET(field,fieldVariableType,fieldVariable,err,error,*999)

    !Now set DOF constraint
    CALL BoundaryConditions_ConstrainDofsEqual(boundaryConditions,fieldVariable,globalDofs,err,error,*999)

    DEALLOCATE(globalDofs)

    EXITS("BoundaryConditions_ConstrainNodeDofsEqual")
    RETURN
999 IF(ALLOCATED(globalDofs)) DEALLOCATE(globalDofs)
998 ERRORSEXITS("BoundaryConditions_ConstrainNodeDofsEqual",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_ConstrainNodeDofsEqual

  !
  !================================================================================================================================
  !

  !>Constrain a DOF to be a linear combination of other DOFs.
  SUBROUTINE BoundaryConditions_DofConstraintSet(boundaryConditions,fieldVariable,globalDof,dofs,coefficients,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions for the solver equations in which to constrain the DOF.
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable containing the DOFs.
    INTEGER(INTG), INTENT(IN) :: globalDof !<The global DOF to set the constraint on.
    INTEGER(INTG), INTENT(IN) :: dofs(:) !<The global DOFs that this DOF is constrained to depend on.
    REAL(DP), INTENT(IN) :: coefficients(:) !<The coefficient values in the DOF constraint.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfDofs,dofIdx,dofIdx2
    TYPE(BoundaryConditionsDofConstraintPtrType), ALLOCATABLE :: newConstraints(:)
    TYPE(Boundary_Conditions_Variable_Type), POINTER :: boundaryConditionsVariable
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint

    NULLIFY(dofConstraint)
    NULLIFY(dofConstraints)

    ENTERS("BoundaryConditions_DofConstraintSet",err,error,*998)

    !Check pointers for association
    IF(.NOT.ASSOCIATED(boundaryConditions)) THEN
      CALL FlagError("Boundary conditions are not associated.",err,error,*998)
    END IF
    IF(boundaryConditions%boundary_conditions_finished) THEN
      CALL FlagError("The boundary conditions have already been finished.",err,error,*998)
    END IF
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      CALL FlagError("Field variable is not associated.",err,error,*998)
    END IF
    CALL boundary_conditions_variable_get(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) THEN
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*998)
    END IF
    dofConstraints=>boundaryConditionsVariable%dofConstraints
    IF(.NOT.ASSOCIATED(dofConstraints)) THEN
      CALL FlagError("Boundary conditions DOF constraints are not associated.",err,error,*998)
    END IF

    numberOfDofs=SIZE(dofs,1)
    IF(numberOfDofs==0) THEN
      CALL FlagError("Empty DOFs list.",err,error,*998)
    ELSE IF(numberOfDofs/=SIZE(coefficients,1)) THEN
      CALL FlagError("Length of coefficients does not match length of DOFs array.",err,error,*998)
    ELSE IF(numberOfDofs>1) THEN
      CALL FlagError("Support for constraining an equations DOF to be depended on multiple "// &
        & "other DOFs is not yet implemented.",err,error,*998)
    END IF

    !Check for duplicate DOFs
    DO dofIdx=1,numberOfDofs
      DO dofIdx2=dofIdx+1,numberOfDofs
        IF(dofs(dofIdx)==dofs(dofIdx2)) THEN
          CALL FlagError("DOF number "//TRIM(NumberToVstring(dofs(dofIdx),"*",err,error))// &
            & " is duplicated in the DOF constraint.",err,error,*998)
        END IF
      END DO
    END DO

    !Check DOFs are free
    DO dofIdx=1,numberOfDofs
      IF(boundaryConditionsVariable%dof_types(dofs(dofIdx))/=BOUNDARY_CONDITION_DOF_FREE) THEN
        CALL FlagError("DOF number "//TRIM(NumberToVstring(dofs(dofIdx),"*",err,error))// &
          & " is not free in the boundary conditions.",err,error,*998)
      END IF
    END DO

    !Allocate new DOF constraints and copy over old constraints
    ALLOCATE(newConstraints(dofConstraints%numberOfConstraints+1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate new DOF constraints array.",err,error,*998)
    IF(dofConstraints%numberOfConstraints>0) THEN
      newConstraints(1:dofConstraints%numberOfConstraints)= &
        & dofConstraints%constraints(1:dofConstraints%numberOfConstraints)
    END IF

    !Set the new DOF constraint
    ALLOCATE(dofConstraint,stat=err)
    IF(err/=0) CALL FlagError("Could not allocate new DOF constraint.",err,error,*999)
    ALLOCATE(dofConstraint%dofs(numberOfDofs),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate constraint DOFs array.",err,error,*999)
    ALLOCATE(dofConstraint%coefficients(numberOfDofs),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate constraint coefficients array.",err,error,*999)
    dofConstraint%globalDof=globalDof
    dofConstraint%numberOfDofs=numberOfDofs
    dofConstraint%dofs(1:numberOfDofs)=dofs(1:numberOfDofs)
    dofConstraint%coefficients(1:numberOfDofs)=coefficients(1:numberOfDofs)

    !Add new DOF constraint to new array
    newConstraints(dofConstraints%numberOfConstraints+1)%ptr=>dofConstraint
    !Replace old DOF constraints with new ones
    CALL MOVE_ALLOC(newConstraints,dofConstraints%constraints)
    dofConstraints%numberOfConstraints=dofConstraints%numberOfConstraints+1

    !Set the DOF type and BC type of the constrained DOF
    CALL BoundaryConditions_SetConditionType(boundaryConditionsVariable,globalDof,BOUNDARY_CONDITION_LINEAR_CONSTRAINT, &
      & err,error,*999)

    EXITS("BoundaryConditions_DofConstraintSet")
    RETURN
999 IF(ASSOCIATED(dofConstraint)) THEN
      IF(ALLOCATED(dofConstraint%dofs)) DEALLOCATE(dofConstraint%dofs)
      IF(ALLOCATED(dofConstraint%coefficients)) DEALLOCATE(dofConstraint%coefficients)
      DEALLOCATE(dofConstraint)
    END IF
    IF(ALLOCATED(newConstraints)) DEALLOCATE(newConstraints)
998 ERRORSEXITS("BoundaryConditions_DofConstraintSet",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_DofConstraintSet

  !
  !================================================================================================================================
  !

  !>Finish the creation of the dof constraints for a boundary conditions variable
  SUBROUTINE BoundaryConditions_DofConstraintsCreateFinish(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to boundary conditions variable to finish the dof constraints for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: constraintIdx,dofIdx,thisDofDomain,otherDofDomain
    INTEGER(INTG) :: globalDof,globalDof2,localDof,localDof2
    INTEGER(INTG) :: numberOfCoupledDofs
    INTEGER(INTG), ALLOCATABLE :: newCoupledGlobalDofs(:),newCoupledLocalDofs(:)
    REAL(DP), ALLOCATABLE :: newCoefficients(:)
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: dofCoupling
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: variableDomainMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable

    ENTERS("BoundaryConditions_DofConstraintsCreateFinish",err,error,*998)

    NULLIFY(dofCoupling)

    !We have a list of DOF constraints, which give the values for a field variable
    !DOF as a linear combination of other DOFs.
    !In order to be able to construct the solver matrices in the solver mapping routines,
    !we create a set of couplings, where a coupling is a set of field variable DOFs
    !mapped to a single solver row or column.

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      fieldVariable=>boundaryConditionsVariable%variable
      IF(ASSOCIATED(fieldVariable)) THEN
        IF(ASSOCIATED(boundaryConditionsVariable%dofConstraints)) THEN
          dofConstraints=>boundaryConditionsVariable%dofConstraints
        ELSE
          CALL FlagError("Boundary conditions DOF constraints are not associated.",err,error,*998)
        END IF

        variableDomainMapping=>fieldVariable%domain_mapping
        IF(.NOT.ASSOCIATED(variableDomainMapping)) THEN
          CALL FlagError("Field variable domain mapping is not associated for variable type "// &
            & TRIM(NumberToVstring(fieldVariable%variable_type,"*",err,error))//".",err,error,*998)
        END IF

        !Allocate an array of pointers to DOF couplings
        IF(dofConstraints%numberOfConstraints>0) THEN
          ALLOCATE(dofConstraints%dofCouplings(fieldVariable%number_of_global_dofs),stat=err)
          IF(err/=0) CALL FlagError( &
            & "Could not allocate dof constraints dof couplings array.",err,error,*998)
          dofConstraints%numberOfDofs=fieldVariable%number_of_global_dofs
          DO dofIdx=1,fieldVariable%number_of_global_dofs
            NULLIFY(dofConstraints%dofCouplings(dofIdx)%ptr)
          END DO
        END IF

        !Loop over all constraints
        DO constraintIdx=1,dofConstraints%numberOfConstraints
          dofConstraint=>dofConstraints%constraints(constraintIdx)%ptr
          IF(.NOT.ASSOCIATED(dofConstraint)) THEN
            CALL FlagError("DOF constraint number "// &
              & TRIM(NumberToVstring(constraintIdx,"*",err,error))// &
              & " is not associated.",err,error,*999)
          END IF

          globalDof=dofConstraint%globalDof
          localDof=variableDomainMapping%global_to_local_map(globalDof)%local_number(1)
          thisDofDomain=variableDomainMapping%global_to_local_map(globalDof)%domain_number(1)

          !Check that the constrained DOFs are still set to be constrained, as
          !subsequently setting a boundary condition would change the DOF type but
          !not update the DOF constraints structure.
          IF(boundaryConditionsVariable%dof_types(globalDof)/=BOUNDARY_CONDITION_DOF_CONSTRAINED) THEN
            CALL FlagError("Global DOF number "//TRIM(NumberToVstring(globalDof,"*",err,error))// &
              & " is part of a linear constraint but the DOF type has been changed"// &
              & " by applying a boundary condition.",err,error,*999)
          END IF

          DO dofIdx=1,dofConstraint%numberOfDofs
            globalDof2=dofConstraint%dofs(dofIdx)
            localDof2=variableDomainMapping%global_to_local_map(globalDof2)%local_number(1)
            !Check a Dirichlet conditions hasn't also been set on this DOF
            IF(boundaryConditionsVariable%dof_types(globalDof2)/=BOUNDARY_CONDITION_DOF_FREE) THEN
              CALL FlagError("A Dirichlet boundary condition has been set on DOF number "// &
                & TRIM(NumberToVstring(globalDof2,"*",err,error))// &
                & " which is part of a linear constraint.",err,error,*999)
            END IF

            !Check we don't have DOF constraints that are split over domains
            !\todo Implement support for DOF constraints that are split over domains
            IF(variableDomainMapping%number_of_domains>1) THEN
              otherDofDomain=variableDomainMapping%global_to_local_map(globalDof2)%domain_number(1)
              IF(thisDofDomain/=otherDofDomain) THEN
                CALL FlagError("An equal DOF constraint is split over multiple domains, "// &
                  & "support for this has not yet been implemented.",err,error,*999)
              END IF
            END IF

            !Add to the DOFs that are coupled with globalDof2
            !This might be quite inefficient if there are many dofs mapped to a single row/column
            !due to the reallocation at each step
            IF(ASSOCIATED(dofConstraints%dofCouplings(globalDof2)%ptr)) THEN
              numberOfCoupledDofs=dofConstraints%dofCouplings(globalDof2)%ptr%numberOfDofs
              ALLOCATE(newCoupledGlobalDofs(numberOfCoupledDofs+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling global DOFs.",err,error,*999)
              ALLOCATE(newCoupledLocalDofs(numberOfCoupledDofs+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling local DOFs.",err,error,*999)
              ALLOCATE(newCoefficients(numberOfCoupledDofs+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling values.",err,error,*999)
              newCoupledGlobalDofs(1:numberOfCoupledDofs)=dofCoupling%globalDofs(1:numberOfCoupledDofs)
              newCoupledLocalDofs(1:numberOfCoupledDofs)=dofCoupling%localDofs(1:numberOfCoupledDofs)
              newCoefficients(1:numberOfCoupledDofs)=dofCoupling%coefficients(1:numberOfCoupledDofs)
            ELSE
              !Set up a a new dofCoupling and set globalDof2 as the first DOF
              ALLOCATE(dofConstraints%dofCouplings(globalDof2)%ptr,stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling type.",err,error,*999)
              ALLOCATE(newCoupledGlobalDofs(2),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling global DOFs.",err,error,*999)
              ALLOCATE(newCoupledLocalDofs(2),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling local DOFs.",err,error,*999)
              ALLOCATE(newCoefficients(2),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling values.",err,error,*999)
              newCoupledGlobalDofs(1)=globalDof2
              newCoupledLocalDofs(1)=localDof2
              newCoefficients(1)=1.0_DP
              numberOfCoupledDofs=1
            END IF
            dofCoupling=>dofConstraints%dofCouplings(globalDof2)%ptr
            newCoupledGlobalDofs(numberOfCoupledDofs+1)=globalDof
            newCoupledLocalDofs(numberOfCoupledDofs+1)=localDof
            newCoefficients(numberOfCoupledDofs+1)=dofConstraint%coefficients(dofIdx)
            CALL MOVE_ALLOC(newCoupledGlobalDofs,dofCoupling%globalDofs)
            CALL MOVE_ALLOC(newCoupledLocalDofs,dofCoupling%localDofs)
            CALL MOVE_ALLOC(newCoefficients,dofCoupling%coefficients)
            dofCoupling%numberOfDofs=numberOfCoupledDofs+1
          END DO !dofIdx
        END DO !constraintIdx
      ELSE
        CALL FlagError("Field variable is not associated for this boundary conditions variable",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_DofConstraintsCreateFinish")
    RETURN
999 IF(ALLOCATED(newCoupledGlobalDofs)) DEALLOCATE(newCoupledGlobalDofs)
    IF(ALLOCATED(newCoupledLocalDofs)) DEALLOCATE(newCoupledLocalDofs)
    IF(ALLOCATED(newCoefficients)) DEALLOCATE(newCoefficients)
    CALL BoundaryConditions_DofConstraintsFinalise(dofConstraints,err,error,*998)
998 ERRORS("BoundaryConditions_DofConstraintsCreateFinish",err,error)
    EXITS("BoundaryConditions_DofConstraintsCreateFinish")
    RETURN 1

  END SUBROUTINE BoundaryConditions_DofConstraintsCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise the DOF constraints structure
  SUBROUTINE BoundaryConditions_DofConstraintsFinalise(dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the dof constraints to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: constraintIdx,dofIdx

    ENTERS("BoundaryConditions_DofConstraintsFinalise",err,error,*999)

    IF(ASSOCIATED(dofConstraints)) THEN
      IF(ALLOCATED(dofConstraints%constraints)) THEN
        DO constraintIdx=1,dofConstraints%numberOfConstraints
          IF(ASSOCIATED(dofConstraints%constraints(constraintIdx)%ptr)) THEN
            IF(ALLOCATED(dofConstraints%constraints(constraintIdx)%ptr%dofs)) THEN
              DEALLOCATE(dofConstraints%constraints(constraintIdx)%ptr%dofs)
            END IF
            IF(ALLOCATED(dofConstraints%constraints(constraintIdx)%ptr%coefficients)) THEN
              DEALLOCATE(dofConstraints%constraints(constraintIdx)%ptr%coefficients)
            END IF
            DEALLOCATE(dofConstraints%constraints(constraintIdx)%ptr)
          END IF
        END DO
        DEALLOCATE(dofConstraints%constraints)
      END IF
      IF(ALLOCATED(dofConstraints%dofCouplings)) THEN
        DO dofIdx=1,dofConstraints%numberOfDofs
          IF(ASSOCIATED(dofConstraints%dofCouplings(dofIdx)%ptr)) THEN
            IF(ALLOCATED(dofConstraints%dofCouplings(dofIdx)%ptr%globalDofs)) THEN
              DEALLOCATE(dofConstraints%dofCouplings(dofIdx)%ptr%globalDofs)
            END IF
            IF(ALLOCATED(dofConstraints%dofCouplings(dofIdx)%ptr%localDofs)) THEN
              DEALLOCATE(dofConstraints%dofCouplings(dofIdx)%ptr%localDofs)
            END IF
            IF(ALLOCATED(dofConstraints%dofCouplings(dofIdx)%ptr%coefficients)) THEN
              DEALLOCATE(dofConstraints%dofCouplings(dofIdx)%ptr%coefficients)
            END IF
          END IF
        END DO
        DEALLOCATE(dofConstraints%dofCouplings)
      END IF
    ELSE
      CALL FlagError("dofConstraints pointer is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_DofConstraintsFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_DofConstraintsFinalise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_DofConstraintsFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the DOF constraints structure
  SUBROUTINE BoundaryConditions_DofConstraintsInitialise(dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the dof constraints to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("BoundaryConditions_DofConstraintsInitialise",err,error,*999)

    IF(ASSOCIATED(dofConstraints)) THEN
      dofConstraints%numberOfConstraints=0
      dofConstraints%numberOfDofs=0
    ELSE
      CALL FlagError("dofConstraints pointer is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_DofConstraintsInitialise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_DofConstraintsInitialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_DofConstraintsInitialise

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

    ENTERS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ALLOCATED(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES))  &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES)
      IF(ALLOCATED(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES))  &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS)) THEN
        BOUNDARY_CONDITIONS_DIRICHLET=>BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS
        CALL BoundaryConditions_SparsityIndicesArrayFinalise(BOUNDARY_CONDITIONS_DIRICHLET% &
            & LINEAR_SPARSITY_INDICES,ERR,ERROR,*999)
        CALL BoundaryConditions_SparsityIndicesArrayFinalise(BOUNDARY_CONDITIONS_DIRICHLET% &
            & DYNAMIC_SPARSITY_INDICES,ERR,ERROR,*999)
        IF(ALLOCATED(BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES)) THEN
          DEALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES)
        ENDIF
        DEALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET)
      ENDIF
      CALL BoundaryConditions_NeumannFinalise(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)) &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)
      IF(ASSOCIATED(boundary_conditions_variable%dofConstraints)) THEN
        CALL BoundaryConditions_DofConstraintsFinalise(boundary_conditions_variable%dofConstraints,err,error,*999)
        DEALLOCATE(boundary_conditions_variable%dofConstraints)
      END IF
      DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalise an array of sparcity indices and deallocate all memory.
  SUBROUTINE BoundaryConditions_SparsityIndicesArrayFinalise(SPARSITY_INDICES_ARRAY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_PTR_TYPE), ALLOCATABLE :: SPARSITY_INDICES_ARRAY(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equ_set_idx, equ_matrix_idx
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES
    
    ENTERS("BoundaryConditions_SparsityIndicesArrayFinalise",ERR,ERROR,*999)
    
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
    
    EXITS("BoundaryConditions_SparsityIndicesArrayFinalise")
    RETURN
999 ERRORS("BoundaryConditions_SparsityIndicesArrayFinalise",ERR,ERROR)
    EXITS("BoundaryConditions_SparsityIndicesArrayFinalise")
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SparsityIndicesArrayFinalise

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

    ENTERS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR,*998)

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
            IF(ERR/=0) CALL FlagError("Could not allocate new boundary conditions variables array.",ERR,ERROR,*998)
            IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
              DO variable_idx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
                NEW_BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR=> &
                    & BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
              ENDDO
            ENDIF

            ALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES(BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate boundary condition variable.",ERR,ERROR,*998)
            BOUNDARY_CONDITIONS_VARIABLE=>NEW_BOUNDARY_CONDITIONS_VARIABLES( &
                & BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1)%PTR
            BOUNDARY_CONDITIONS_VARIABLE%BOUNDARY_CONDITIONS=>BOUNDARY_CONDITIONS
            BOUNDARY_CONDITIONS_VARIABLE%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            BOUNDARY_CONDITIONS_VARIABLE%VARIABLE=>FIELD_VARIABLE
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate global boundary condition types.",ERR,ERROR,*999)
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES(VARIABLE_DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate global boundary condition dof types.",ERR,ERROR,*999)
            BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES=BOUNDARY_CONDITION_FREE
            BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES=BOUNDARY_CONDITION_DOF_FREE
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(MAX_BOUNDARY_CONDITION_NUMBER),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate boundary condition DOF counts array.",ERR,ERROR,*999)
            BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS=0
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS)
            BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS=0
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%neumannBoundaryConditions)
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%parameterSetRequired(FIELD_NUMBER_OF_SET_TYPES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate boundary condition parameter set required array.",ERR,ERROR,*999)
            BOUNDARY_CONDITIONS_VARIABLE%parameterSetRequired=.FALSE.
            BOUNDARY_CONDITIONS_VARIABLE%parameterSetRequired(FIELD_VALUES_SET_TYPE)=.TRUE.

            CALL MOVE_ALLOC(NEW_BOUNDARY_CONDITIONS_VARIABLES,BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)
            BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES= &
                & BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1

            ALLOCATE(boundary_conditions_variable%DofConstraints,stat=err)
            IF(err/=0) CALL FlagError("Could not allocate boundary conditions dof constraints.",err,error,*999)
            CALL BoundaryConditions_DofConstraintsInitialise(boundary_conditions_variable%DofConstraints,err,error,*999)
          END IF
        ELSE
          CALL FlagError("Field variable domain mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FlagError("Field variable is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*998)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS_VARIABLE,DUMMY_ERR,DUMMY_ERROR,*998)
    DEALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES)
998 ERRORSEXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_VARIABLE_GET",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
          VARIABLE_FOUND=.FALSE.
          variable_idx=1
          DO WHILE(variable_idx<=BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES.AND..NOT.VARIABLE_FOUND)
            VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR%VARIABLE
            IF(ASSOCIATED(VARIABLE)) THEN
              IF(VARIABLE%VARIABLE_TYPE==FIELD_VARIABLE%VARIABLE_TYPE.AND. &
                & VARIABLE%FIELD%USER_NUMBER==FIELD_VARIABLE%FIELD%USER_NUMBER) THEN
                IF(ASSOCIATED(VARIABLE%FIELD%REGION)) THEN
                  IF(VARIABLE%FIELD%REGION%USER_NUMBER==FIELD_VARIABLE%FIELD%REGION%USER_NUMBER) THEN
                    VARIABLE_FOUND=.TRUE.
                    BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
                  ENDIF
                ELSEIF(ASSOCIATED(VARIABLE%FIELD%INTERFACE)) THEN
                  IF(VARIABLE%FIELD%INTERFACE%USER_NUMBER==FIELD_VARIABLE%FIELD%INTERFACE%USER_NUMBER) THEN
                    VARIABLE_FOUND=.TRUE.
                    BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR
                  ENDIF
                ENDIF
              ENDIF
              variable_idx=variable_idx+1
            ENDIF
          ENDDO
        ENDIF
      ELSE
        CALL FlagError("Field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_VARIABLE_GET")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_VARIABLE_GET",ERR,ERROR)
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

    ENTERS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS)) THEN
        CALL FlagError("Dirichlet boundary conditions are already associated for this boundary conditions variable." &
           & ,ERR,ERROR,*999)
      ELSE
        ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate Dirichlet Boundary Conditions",ERR,ERROR,*999)
        BOUNDARY_CONDITIONS_DIRICHLET=>BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS
        NUMBER_OF_DIRICHLET_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
        ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%DIRICHLET_DOF_INDICES(NUMBER_OF_DIRICHLET_CONDITIONS),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate Dirichlet DOF indices array",ERR,ERROR,*999)

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
                  CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO
          ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES(SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS, &
                & MAX_NUMBER_LINEAR_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Dirichlet linear sparsity indices array",ERR,ERROR,*999)
          ALLOCATE(BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES(SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,&
                & MAX_NUMBER_DYNAMIC_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Dirichlet dynamic sparsity indices array",ERR,ERROR,*999)
          DO equations_set_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            DO matrix_idx=1,MAX_NUMBER_LINEAR_MATRICES
              NULLIFY(BOUNDARY_CONDITIONS_DIRICHLET%LINEAR_SPARSITY_INDICES(equations_set_idx,matrix_idx)%PTR)
            ENDDO
            DO matrix_idx=1,MAX_NUMBER_DYNAMIC_MATRICES
              NULLIFY(BOUNDARY_CONDITIONS_DIRICHLET%DYNAMIC_SPARSITY_INDICES(equations_set_idx,matrix_idx)%PTR)
            ENDDO
          ENDDO
        ELSE
          CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE")
    RETURN
!!TODO \todo write BOUNDARY_CONDITIONS_DIRICHLET_FINALISE
999 ERRORSEXITS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialise Sparsity Indices type
  SUBROUTINE BoundaryConditions_SparsityIndicesInitialise(SPARSITY_INDICES,NUMBER_OF_DIRICHLET,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES !<A pointer to the Sparsity Indices type tp initialise
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIRICHLET !<The number of dirichlet conditions this sparsity indices type will hold
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_SparsityIndicesInitialise",ERR,ERROR,*999)

    IF(ASSOCIATED(SPARSITY_INDICES)) THEN
     CALL FlagError("Sparsity Indices are already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(SPARSITY_INDICES,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate sparsity indicies.",ERR,ERROR,*999)
      ALLOCATE(SPARSITY_INDICES%SPARSE_COLUMN_INDICES(NUMBER_OF_DIRICHLET+1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate sparsity column indices array",ERR,ERROR,*999)
    ENDIF

    EXITS("BoundaryConditions_SparsityIndicesInitialise")
    RETURN
!!TODO \todo write BOUNDARY_CONDITIONS_SPARSITY_INDICES_FINALISE
999 ERRORS("BoundaryConditions_SparsityIndicesInitialise",ERR,ERROR)
    EXITS("BoundaryConditions_SparsityIndicesInitialise")
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SparsityIndicesInitialise

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

    ENTERS("BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)) THEN
        CALL FlagError("Pressure incremented boundary conditions are already associated for this boundary conditions variable." &
           & ,ERR,ERROR,*999)
      ELSE
        ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate Pressure incremented Boundary Conditions",ERR,ERROR,*999)
        BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED=>BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS
        NUMBER_OF_PRESSURE_INCREMENTED_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        ALLOCATE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED%PRESSURE_INCREMENTED_DOF_INDICES &
          & (NUMBER_OF_PRESSURE_INCREMENTED_CONDITIONS),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate Pressure incremented DOF indices array",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE")
    RETURN
!!TODO \todo write BOUNDARY_CONDITIONS_DIRICHLET_FINALISE
999 ERRORSEXITS("BOUNDARY_CONDITIONS_DIRICHLET_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE BOUNDARY_CONDITIONS_ROUTINES
